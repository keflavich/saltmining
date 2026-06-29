"""Download ALMA 2017.1.01102.S data for OrionB-Flame.

Per Adam: ~60 AU resolution at d=410 pc -> 0.15" beam. Add to uvdata/.
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astroquery.alma import Alma
import numpy as np

OUTDIR = Path("/orange/adamginsburg/salt/survey_2026/uvdata/2017.1.01102.S/OrionB-Flame")
OUTDIR.mkdir(parents=True, exist_ok=True)

TARGET_RA = 85.42917
TARGET_DEC = -1.84167
SEARCH_RAD_ARCSEC = 60.0


def main():
    print("Querying ALMA archive for 2017.1.01102.S...", flush=True)
    q = Alma.query(payload={"project_code": "2017.1.01102.S"},
                     public=None, science=True)
    print(f"  {len(q)} obs entries", flush=True)
    if len(q) == 0:
        return
    ra = np.asarray(q["s_ra"], dtype=float)
    dec = np.asarray(q["s_dec"], dtype=float)
    sep = np.sqrt(((ra - TARGET_RA) * np.cos(np.radians(TARGET_DEC))) ** 2
                   + (dec - TARGET_DEC) ** 2) * 3600
    mask = sep < SEARCH_RAD_ARCSEC
    print(f"  {mask.sum()} obs within {SEARCH_RAD_ARCSEC}\" of target",
          flush=True)
    if mask.sum() == 0:
        return
    sub = q[mask]
    uids = sorted(set(sub["member_ous_uid"]))
    print(f"  unique MOUSes: {len(uids)}", flush=True)
    for uid in uids:
        print(f"    {uid}", flush=True)

    # Download the data
    try:
        urls = Alma.get_data_info(uids, expand_tarfiles=True)
    except (ConnectionError, OSError, ValueError) as e:
        print(f"  url fetch err: {e}", flush=True)
        return
    print(f"  {len(urls)} files to fetch", flush=True)
    # Filter to pbcor + mfs / cube fits
    keep_urls = [u for u in urls["access_url"]
                 if (".pbcor.fits" in u and ("cube" in u or "cont" in u))]
    print(f"  filtered: {len(keep_urls)} pbcor fits files", flush=True)
    Alma.cache_location = str(OUTDIR)
    # Use download_files
    if keep_urls:
        files = Alma.download_files(keep_urls, savedir=str(OUTDIR),
                                       continuation=True)
        print(f"  downloaded: {len(files) if files else 0}")


if __name__ == "__main__":
    main()
