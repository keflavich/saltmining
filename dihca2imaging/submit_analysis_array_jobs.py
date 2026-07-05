#!/usr/bin/env python3
"""Create and submit a SLURM array job with one task per image cube.

This script enumerates image cubes (FITS) in a directory, writes a filelist,
creates a small SBATCH script that maps $SLURM_ARRAY_TASK_ID to a single
filename, substitutes that filename into a user-provided command template,
and submits the array job.

Example:
  # Dry-run, just print the generated SBATCH script
  python submit_array_jobs.py --image-dir /data/images --command "echo processing {file}" --dry-run

  # Submit an array job, limit concurrent tasks to 20
  python submit_array_jobs.py --image-dir /data/images --command "python analyze_cube.py --cube {file} --outdir /scratch/out" --max-concurrent 20

Notes:
- The command template must include the token `{file}` which will be replaced
  by the absolute path to each image for the corresponding array task.
- Default SLURM account/QoS come from the DIHCA README; override with CLI flags.
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import textwrap
from pathlib import Path


DEFAULT_ACCOUNT = "astronomy-dept"
DEFAULT_QOS = "astronomy-dept-b"


def build_sbatch_script(job_name: str, out_dir: str, filelist_path: str, n_tasks: int,
                        command_template: str, account: str, qos: str, time: str,
                        mem: str, cpus: int, max_concurrent: int | None) -> str:
    """Return SBATCH script text that runs `command_template` with {file}.

    The generated script reads the N-th line from `filelist_path` and substitutes
    it for `{file}` in the `command_template` before executing.
    """

    array_spec = f"1-{n_tasks}"
    if max_concurrent and max_concurrent > 0:
        array_spec = f"{array_spec}%{max_concurrent}"

    sbatch = textwrap.dedent(f"""
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --output={os.path.join(out_dir, job_name + '_%A_%a_%j.out')}
    #SBATCH --error={os.path.join(out_dir, job_name + '_%A_%a_%j.err')}
    #SBATCH --account={account}
    #SBATCH --qos={qos}
    #SBATCH --time={time}
    #SBATCH --mem={mem}
    #SBATCH --cpus-per-task={cpus}
    #SBATCH --array={array_spec}

    FILELIST={filelist_path}

    # Map the SLURM array task id to a file (1-based index)
    FILE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" "$FILELIST")

    # The command template is substituted so that occurrences of {{file}} are
    # replaced by the selected filename. Using sed ensures robust quoting.
    COMMAND='{command_template}'
    EXEC_CMD=$(echo "$COMMAND" | sed "s|{{file}}|$FILE|g")

    echo "[array task ${{SLURM_ARRAY_TASK_ID}}] running: $EXEC_CMD"
    eval "$EXEC_CMD"
    """.lstrip("\n"))

    return sbatch


def main():
    parser = argparse.ArgumentParser(description="Submit SLURM array with one task per image cube")
    parser.add_argument("--image-dir", required=True, help="Directory containing image cubes (FITS)")
    parser.add_argument("--pattern", default="*_spw*.image.fits",
                        help="Glob pattern to find image cubes (default: '*_spw*.image.fits')")
    parser.add_argument("--command", required=False,
                        help="Command template to run for each file. Must include the token {file}")
    parser.add_argument("--spectral-analysis", action="store_true",
                        help="Shortcut: run the repository's spectral_analysis.py for each cube (overrides --command)")
    parser.add_argument("--job-name", default="dihca_array", help="SLURM job name")
    parser.add_argument("--out-dir", default="./slurm_out", help="Directory for sbatch stdout/err and filelist")
    parser.add_argument("--account", default=DEFAULT_ACCOUNT)
    parser.add_argument("--qos", default=DEFAULT_QOS)
    parser.add_argument("--time", default="48:00:00", help="SLURM time limit")
    parser.add_argument("--mem", default="16G", help="Memory per job (e.g. 16G or 16000M)")
    parser.add_argument("--cpus", default=1, type=int, help="CPUs per task")
    parser.add_argument("--max-concurrent", type=int, default=0,
                        help="Limit concurrent array tasks (e.g. 20). 0 for no limit")
    parser.add_argument("--dry-run", action="store_true", help="Print SBATCH script and do not submit")
    parser.add_argument("--submit", action="store_true", help="Actually submit the job (overrides dry-run)")

    args = parser.parse_args()

    # If requested, build the command to run the spectral analysis script per-cube
    if args.spectral_analysis:
        script_path = Path(__file__).resolve().parent / "spectral_analysis.py"
        if not script_path.exists():
            raise SystemExit(f"spectral_analysis.py not found at expected location: {script_path}")
        # Override the command template to invoke the script with --process-cube {file}
        args.command = f"python {script_path} --process-cube {{file}}"

    if not args.command:
        raise SystemExit("Either --command or --spectral-analysis must be provided")

    image_dir = Path(args.image_dir).expanduser()
    if not image_dir.exists():
        raise SystemExit(f"Image directory does not exist: {image_dir}")

    pattern = str(image_dir / args.pattern)
    files = sorted(glob.glob(pattern))
    if len(files) == 0:
        raise SystemExit(f"No files matched pattern: {pattern}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filelist_path = out_dir / f"{args.job_name}_filelist.txt"
    with open(filelist_path, "w") as fh:
        for fn in files:
            fh.write(os.path.abspath(fn) + "\n")

    n_tasks = len(files)
    max_concurrent = args.max_concurrent if args.max_concurrent > 0 else None

    sbatch_text = build_sbatch_script(
        job_name=args.job_name,
        out_dir=str(out_dir),
        filelist_path=str(filelist_path),
        n_tasks=n_tasks,
        command_template=args.command,
        account=args.account,
        qos=args.qos,
        time=args.time,
        mem=args.mem,
        cpus=args.cpus,
        max_concurrent=max_concurrent,
    )

    sbatch_path = out_dir / f"{args.job_name}_array.sh"
    with open(sbatch_path, "w") as fh:
        fh.write(sbatch_text)

    print(f"Wrote filelist ({n_tasks} entries) -> {filelist_path}")
    print(f"Wrote SBATCH script -> {sbatch_path}")

    if args.dry_run and not args.submit:
        print("--- SBATCH script (dry-run) ---")
        print(sbatch_text)
        return

    # Submit the job
    submit_cmd = ["sbatch", str(sbatch_path)]
    print(f"Submitting SBATCH: {' '.join(submit_cmd)}")
    proc = subprocess.run(submit_cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        print("Failed to submit job:")
        print(proc.stderr)
        raise SystemExit(1)

    print(proc.stdout.strip())


if __name__ == '__main__':
    main()
