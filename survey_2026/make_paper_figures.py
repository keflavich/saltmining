"""Build paper figures for demography_2026 from the current sample +
analysis_products. Outputs go to demography_2026/figures/.

Figures
-------
detect_vs_lbol_{species}.pdf    -- logistic regression P(detect) vs log10 L_bol
                                    for species in {NaCl, KCl, H2O, RRL, SiO, SiS, SO}
detect_corner_logistic.pdf      -- 6x6 logistic regression corner over species
intensity_corner.pdf            -- pairwise intensity / UL scatter for all
                                    species pairs
cooccurrence_heatmap.pdf        -- P(B|A) matrix
detect_vs_distance.pdf          -- logistic regression P(detect) vs distance
                                    (per species, faceted)
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.optimize import minimize

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
FIG_DIR = Path("/orange/adamginsburg/salt/demography_2026/figures")
FIG_DIR.mkdir(exist_ok=True)

SPECIES_GROUPS = {
    "NaCl": ("NaCl", lambda L: L.startswith("NaCl_")),
    "KCl":  ("KCl",  lambda L: L.startswith("KCl_")),
    "H2O":  ("H2O",  lambda L: L.startswith("H2O")),
    "RRL":  ("RRL",  lambda L: bool(re.match(r"^H\d+(alpha|beta|gamma|delta)$", L))),
    "SiO":  ("SiO",  lambda L: L.startswith("SiO") or L.startswith("29SiO")),
    "SiS":  ("SiS",  lambda L: L.startswith("SiS") or L.startswith("29SiS")),
    "SO":   ("SO",   lambda L: L.startswith("SO_") or L.startswith("SO2_") or L == "SO"),
}


def brightest_id(proposal_dir):
    cont = proposal_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def load_per_target_detections():
    """Return DataFrame: one row per target with columns
    {species: True/False detected at brightest source, species_peak_K,
    species_int, species_upper_K, species_in_band}
    plus dist_kpc, lbol_lsun, name.
    """
    src = pd.read_csv(SRC_CSV)
    rows = []
    for _, sr in src.iterrows():
        name = sr["name"]
        dist = float(sr["dist_kpc"])
        lbol = float(sr["lbol_lsun"])
        adirs = sorted((ANALYSIS / name).glob("2*")) if (ANALYSIS / name).is_dir() else []
        # Aggregate across all proposals at the brightest source: max(snr) per species.
        agg = {sp: dict(snr=np.nan, peak=np.nan, integ=np.nan, sigma=np.nan,
                         in_band=False) for sp in SPECIES_GROUPS}
        for ad in adirs:
            bid = brightest_id(ad)
            if bid is None:
                continue
            meas = ad / "line_measurements.csv"
            if not meas.exists():
                continue
            try:
                m = pd.read_csv(meas)
            except pd.errors.EmptyDataError:
                continue
            if m.empty or not {"source", "snr", "line"} <= set(m.columns):
                continue
            sub = m[m["source"] == bid]
            for sp, (_, pred) in SPECIES_GROUPS.items():
                hits = sub[sub["line"].astype(str).apply(pred)]
                if hits.empty:
                    continue
                agg[sp]["in_band"] = True
                best = hits.nlargest(1, "snr")
                snr = float(best.iloc[0]["snr"])
                if not np.isfinite(agg[sp]["snr"]) or snr > agg[sp]["snr"]:
                    agg[sp].update(
                        snr=snr,
                        peak=float(best.iloc[0].get("peak_Kkms_or_unit", np.nan)),
                        integ=float(best.iloc[0].get("integ", np.nan)),
                        sigma=float(best.iloc[0].get("sigma", np.nan)),
                    )
        row = dict(name=name, dist_kpc=dist, lbol_lsun=lbol,
                    analyzed=any(a[1].is_dir() for a in [(0, ANALYSIS / name)]))
        for sp, info in agg.items():
            row[f"{sp}_in_band"] = info["in_band"]
            row[f"{sp}_snr"] = info["snr"]
            row[f"{sp}_peak"] = info["peak"]
            row[f"{sp}_integ"] = info["integ"]
            row[f"{sp}_sigma"] = info["sigma"]
            row[f"{sp}_det"] = bool(np.isfinite(info["snr"]) and info["snr"] >= 5.0)
            row[f"{sp}_intensity_or_UL"] = (
                info["integ"] if info["snr"] >= 5.0 else (3.0 * info["sigma"])
            )
        rows.append(row)
    df = pd.DataFrame(rows)
    # Merge audit warning flags
    audit_csv = ROOT / "data/evidence_audit.csv"
    if audit_csv.exists():
        a = pd.read_csv(audit_csv)
        if not a.empty and "flags" in a.columns:
            flags_per_target = (a.groupby("target")["flags"]
                                  .apply(lambda s: "|".join(sorted(set(
                                      f for v in s.dropna() for f in str(v).split("|") if f
                                  ))))
                                  .reset_index().rename(columns={"flags": "warning_flags"}))
            df = df.merge(flags_per_target, left_on="name", right_on="target", how="left")
            df["warning_flags"] = df["warning_flags"].fillna("")
            # Per-species suspicious flag: TRUE if det is single-line OR vexc-only
            for sp in SPECIES_GROUPS:
                col = f"{sp}_suspect"
                df[col] = df["warning_flags"].apply(
                    lambda fs: bool(fs) and any(
                        tag in fs for tag in (
                            f"{sp}_VEXC_ONLY",
                            f"{sp}_SINGLE_LINE_DET",
                            "H2O_NONDET_WITH_NACL_DET" if sp == "NaCl" else "",
                            "KCL_DET_NACL_NONDET" if sp == "KCl" else "",
                            "HIGH_CONFUSION",
                        ) if tag
                    ))
    else:
        df["warning_flags"] = ""
        for sp in SPECIES_GROUPS:
            df[f"{sp}_suspect"] = False
    return df


def logistic_fit(x, y):
    """Fit P(y=1) = expit(a + b*x). Returns (a, b)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]
    if x.size < 4 or np.unique(y).size < 2:
        return None, None
    def nll(params):
        a, b = params
        p = expit(a + b * x)
        p = np.clip(p, 1e-9, 1 - 1e-9)
        return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
    res = minimize(nll, [0.0, 0.0], method="L-BFGS-B")
    if not res.success:
        return None, None
    return float(res.x[0]), float(res.x[1])


def plot_logistic(ax, x, y, title="", xlabel="$\\log_{10} L_\\mathrm{bol}/L_\\odot$",
                   show_pts=True):
    a, b = logistic_fit(x, y)
    ax.set_title(title)
    if show_pts:
        # jitter the 0/1 dots a little
        rng = np.random.default_rng(0)
        jit = rng.uniform(-0.02, 0.02, size=len(y))
        ax.scatter(x, np.asarray(y) + jit, s=22, c="black", alpha=0.7)
    if a is not None:
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 200)
        ax.plot(xx, expit(a + b * xx), "r-", lw=1.5,
                 label=f"a={a:.2f}, b={b:.2f}")
        ax.legend(fontsize=7, loc="center right")
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["non-det", "det"])
    ax.set_xlabel(xlabel)


def fig_detect_vs_lbol(df, species_list, robust=False, name_suffix=""):
    fig, axes = plt.subplots(1, len(species_list), figsize=(3.4 * len(species_list), 3.4),
                              sharey=True)
    if len(species_list) == 1:
        axes = [axes]
    use = df.copy()
    use["lbol_log"] = np.log10(use["lbol_lsun"])
    for ax, sp in zip(axes, species_list):
        sub = use[use[f"{sp}_in_band"]].copy()
        det_col = sub[f"{sp}_det"].astype(bool)
        suspect_col = sub.get(f"{sp}_suspect", pd.Series([False]*len(sub),
                                                            index=sub.index)).astype(bool)
        if robust:
            # Demote suspect detections to non-detection
            det_col = det_col & ~suspect_col
        y = det_col.astype(int).values
        x = sub["lbol_log"].values
        plot_logistic(ax, x, y, title=f"{sp} ({len(sub)} obs)")
        # Mark suspect points with red ring
        suspect = sub[sub.get(f"{sp}_suspect", False) & sub[f"{sp}_det"]]
        if not suspect.empty:
            ax.scatter(np.log10(suspect["lbol_lsun"]),
                        [1.0] * len(suspect), s=60, marker="o",
                        facecolors="none", edgecolors="red", linewidths=1.5,
                        label="suspect")
            ax.legend(fontsize=6, loc="lower right")
    axes[0].set_ylabel("Detection")
    fig.tight_layout()
    out = FIG_DIR / f"detect_vs_lbol{name_suffix}.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_detect_corner(df, species_list):
    n = len(species_list)
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n, 2.4 * n))
    for i, sp_y in enumerate(species_list):
        for j, sp_x in enumerate(species_list):
            ax = axes[i, j]
            if i == j:
                # diagonal: histogram of det/nondet
                vals = df[f"{sp_y}_det"].astype(int)
                ax.hist(vals, bins=[-0.5, 0.5, 1.5], rwidth=0.7, color="gray")
                ax.set_xticks([0, 1])
                ax.set_xticklabels(["n", "y"])
                ax.set_title(sp_y, fontsize=9)
                continue
            if j > i:
                ax.axis("off")
                continue
            # off-diagonal: P(sp_y det | sp_x det) logistic-ish scatter on
            # x = sp_x intensity, y = sp_y detection
            sub = df[df[f"{sp_x}_in_band"] & df[f"{sp_y}_in_band"]]
            x = np.log10(np.maximum(sub[f"{sp_x}_intensity_or_UL"].abs(), 1e-6))
            y = sub[f"{sp_y}_det"].astype(int).values
            plot_logistic(ax, x, y, title="", xlabel=f"log {sp_x} int", show_pts=True)
            ax.set_ylabel(f"{sp_y} det")
    fig.tight_layout()
    out = FIG_DIR / "detect_corner_logistic.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_intensity_corner(df, species_list):
    n = len(species_list)
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n, 2.4 * n))
    for i, sp_y in enumerate(species_list):
        for j, sp_x in enumerate(species_list):
            ax = axes[i, j]
            if i == j:
                ax.axis("off")
                ax.text(0.5, 0.5, sp_y, ha="center", va="center", fontsize=12,
                         transform=ax.transAxes)
                continue
            if j > i:
                ax.axis("off")
                continue
            sub = df[df[f"{sp_x}_in_band"] & df[f"{sp_y}_in_band"]]
            for _, r in sub.iterrows():
                x_val = abs(r[f"{sp_x}_intensity_or_UL"])
                y_val = abs(r[f"{sp_y}_intensity_or_UL"])
                xd = bool(r[f"{sp_x}_det"])
                yd = bool(r[f"{sp_y}_det"])
                marker = "o" if xd and yd else ("^" if xd or yd else "v")
                color = "k" if xd and yd else ("blue" if xd else "red")
                ax.loglog(x_val, y_val, marker=marker, color=color, ms=4,
                          markeredgewidth=0.4)
            ax.set_xlabel(f"{sp_x}")
            ax.set_ylabel(f"{sp_y}")
            ax.tick_params(labelsize=6)
    fig.tight_layout()
    out = FIG_DIR / "intensity_corner.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_cooccurrence(df, species_list):
    n = len(species_list)
    mat = np.full((n, n), np.nan)
    for i, sp_a in enumerate(species_list):
        det_a = df[df[f"{sp_a}_in_band"]][f"{sp_a}_det"]
        if det_a.sum() == 0:
            continue
        for j, sp_b in enumerate(species_list):
            sub = df[df[f"{sp_a}_in_band"] & df[f"{sp_b}_in_band"] & df[f"{sp_a}_det"]]
            if sub.empty:
                continue
            mat[i, j] = sub[f"{sp_b}_det"].mean()
    fig, ax = plt.subplots(figsize=(0.7 * n + 2, 0.7 * n + 1.5))
    im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(species_list); ax.set_yticklabels(species_list)
    ax.set_xlabel("Detected B")
    ax.set_ylabel("Conditioned on A detected")
    for i in range(n):
        for j in range(n):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                         color="white" if mat[i, j] < 0.5 else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="$P(B \\mid A)$")
    fig.tight_layout()
    out = FIG_DIR / "cooccurrence_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_detect_vs_distance(df, species_list):
    fig, axes = plt.subplots(1, len(species_list), figsize=(3.4 * len(species_list), 3.4),
                              sharey=True)
    if len(species_list) == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_list):
        sub = df[df[f"{sp}_in_band"]]
        y = sub[f"{sp}_det"].astype(int).values
        x = sub["dist_kpc"].values
        plot_logistic(ax, x, y, title=f"{sp} ({len(sub)} obs)",
                       xlabel="distance (kpc)")
    axes[0].set_ylabel("Detection")
    fig.tight_layout()
    out = FIG_DIR / "detect_vs_distance.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def load_continuum_peak(name):
    """Max peak_Jybeam across all proposals of this target."""
    p_dir = ANALYSIS / name
    if not p_dir.is_dir():
        return np.nan
    peaks = []
    for prop_dir in p_dir.glob("2*"):
        cont = prop_dir / "continuum_sources.csv"
        if not cont.exists():
            continue
        try:
            d = pd.read_csv(cont)
        except pd.errors.EmptyDataError:
            continue
        if d.empty or "peak_Jybeam" not in d.columns:
            continue
        peaks.append(float(d["peak_Jybeam"].max()))
    return float(np.nanmax(peaks)) if peaks else np.nan


def load_f5sigma(name):
    aux = ROOT / "data/data_summary_aux.csv"
    if not aux.exists():
        return np.nan
    d = pd.read_csv(aux)
    sub = d[d["target"] == name]
    if sub.empty:
        return np.nan
    return float(sub["f5sigma"].max())


def load_effective_beam_au(name):
    """Smallest synthesized beam major axis (in AU) across proposals at the
    brightest source of this target."""
    p_dir = ANALYSIS / name
    if not p_dir.is_dir():
        return np.nan
    beams_au = []
    src = pd.read_csv(SRC_CSV)
    dist = float(src.loc[src["name"] == name, "dist_kpc"].iloc[0]) \
        if (src["name"] == name).any() else np.nan
    if not np.isfinite(dist):
        return np.nan
    for prop_dir in p_dir.glob("2*"):
        meas = prop_dir / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            m = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if m.empty or "cube" not in m.columns:
            continue
        from astropy.io import fits
        for cube_name in m["cube"].dropna().unique():
            cp = ROOT / "uvdata" / prop_dir.name / name / str(cube_name)
            if not cp.exists():
                continue
            h = fits.getheader(cp, ext=0)
            if "BMAJ" in h:
                beams_au.append(float(h["BMAJ"]) * 3600.0 * dist * 1000.0)
            break  # one cube is enough
    return float(np.nanmin(beams_au)) if beams_au else np.nan


def fig_detect_vs_continuum(df, species_list):
    """Logistic regression P(detect) vs log10 peak continuum flux (mJy/beam)."""
    df = df.copy()
    df["cont_peak_Jybeam"] = [load_continuum_peak(n) for n in df["name"]]
    df["cont_peak_log"] = np.log10(df["cont_peak_Jybeam"] * 1000.0)
    fig, axes = plt.subplots(1, len(species_list), figsize=(3.4 * len(species_list), 3.4),
                              sharey=True)
    if len(species_list) == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_list):
        sub = df[df[f"{sp}_in_band"] & np.isfinite(df["cont_peak_log"])]
        y = sub[f"{sp}_det"].astype(int).values
        x = sub["cont_peak_log"].values
        plot_logistic(ax, x, y, title=f"{sp} ({len(sub)} obs)",
                       xlabel=r"$\log_{10}\,S_\nu^\mathrm{peak}$ (mJy/beam)")
    axes[0].set_ylabel("Detection")
    fig.tight_layout()
    out = FIG_DIR / "detect_vs_continuum.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_detect_vs_f5sigma(df, species_list):
    df = df.copy()
    df["f5sigma"] = [load_f5sigma(n) for n in df["name"]]
    fig, axes = plt.subplots(1, len(species_list), figsize=(3.4 * len(species_list), 3.4),
                              sharey=True)
    if len(species_list) == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_list):
        sub = df[df[f"{sp}_in_band"] & np.isfinite(df["f5sigma"])]
        y = sub[f"{sp}_det"].astype(int).values
        x = sub["f5sigma"].values
        plot_logistic(ax, x, y, title=f"{sp} ({len(sub)} obs)",
                       xlabel=r"$f(>5\sigma)$ line-confusion fraction")
    axes[0].set_ylabel("Detection")
    fig.tight_layout()
    out = FIG_DIR / "detect_vs_f5sigma.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_detect_vs_beam_au(df, species_list):
    df = df.copy()
    df["beam_au"] = [load_effective_beam_au(n) for n in df["name"]]
    df["beam_au_log"] = np.log10(df["beam_au"])
    fig, axes = plt.subplots(1, len(species_list), figsize=(3.4 * len(species_list), 3.4),
                              sharey=True)
    if len(species_list) == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_list):
        sub = df[df[f"{sp}_in_band"] & np.isfinite(df["beam_au_log"])]
        y = sub[f"{sp}_det"].astype(int).values
        x = sub["beam_au_log"].values
        plot_logistic(ax, x, y, title=f"{sp} ({len(sub)} obs)",
                       xlabel=r"$\log_{10}\,\theta_\mathrm{beam}$ (AU)")
    axes[0].set_ylabel("Detection")
    fig.tight_layout()
    out = FIG_DIR / "detect_vs_beam_au.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_phi_correlation(df, species_list):
    """Phi correlation bar chart: matthews correlation between species detection flags.
    phi = (n11*n00 - n10*n01) / sqrt((n1.)(n0.)(n.1)(n.0))."""
    pairs = []
    for i, sp_a in enumerate(species_list):
        for sp_b in species_list[i + 1:]:
            sub = df[df[f"{sp_a}_in_band"] & df[f"{sp_b}_in_band"]]
            if len(sub) < 4:
                continue
            a = sub[f"{sp_a}_det"].astype(int).values
            b = sub[f"{sp_b}_det"].astype(int).values
            n11 = int(((a == 1) & (b == 1)).sum())
            n10 = int(((a == 1) & (b == 0)).sum())
            n01 = int(((a == 0) & (b == 1)).sum())
            n00 = int(((a == 0) & (b == 0)).sum())
            denom = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)
            if denom == 0:
                continue
            phi = (n11 * n00 - n10 * n01) / np.sqrt(denom)
            pairs.append((f"{sp_a}-{sp_b}", phi, len(sub)))
    pairs.sort(key=lambda r: r[1], reverse=True)
    labels = [p[0] for p in pairs]
    vals = [p[1] for p in pairs]
    ns = [p[2] for p in pairs]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(labels))))
    colors = ["tab:blue" if v >= 0 else "tab:red" for v in vals]
    ax.barh(range(len(labels)), vals, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    for k, (v, n) in enumerate(zip(vals, ns)):
        ax.text(v + (0.01 if v >= 0 else -0.01), k, f" n={n}",
                 va="center", ha="left" if v >= 0 else "right", fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"$\phi$ correlation coefficient")
    ax.set_title("Pairwise species-detection $\\phi$ correlation")
    ax.invert_yaxis()
    fig.tight_layout()
    out = FIG_DIR / "phi_correlation.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    df = load_per_target_detections()
    df.to_csv(ROOT / "data/paper_figure_inputs.csv", index=False)
    print(f"loaded {len(df)} targets")
    species_list = list(SPECIES_GROUPS.keys())
    fig_detect_vs_lbol(df, ["NaCl", "H2O", "RRL"])
    fig_detect_vs_lbol(df, ["NaCl", "H2O", "RRL"], robust=True,
                         name_suffix="_robust")
    fig_detect_vs_lbol(df, species_list)  # extended version
    fig_detect_corner(df, species_list)
    fig_intensity_corner(df, species_list)
    fig_cooccurrence(df, species_list)
    fig_detect_vs_distance(df, species_list)
    fig_detect_vs_continuum(df, ["NaCl", "H2O", "RRL"])
    fig_detect_vs_f5sigma(df, ["NaCl", "H2O", "RRL"])
    fig_detect_vs_beam_au(df, ["NaCl", "H2O", "RRL"])
    fig_phi_correlation(df, species_list)
    print(f"\nfigures in {FIG_DIR}")


if __name__ == "__main__":
    main()
