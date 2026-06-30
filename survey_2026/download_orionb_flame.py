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

TARGET_RA = 85.4188
TARGET_DEC = -1.9045
SEARCH_RAD_ARCSEC = 600.0


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

    # Download the data — retry on 502 / proxy errors
    import time as _t
    for attempt in range(5):
        try:
            urls = Alma.get_data_info(uids, expand_tarfiles=True)
            break
        except Exception as e:
            print(f"  attempt {attempt+1}: url fetch err: {e}", flush=True)
            if attempt < 4:
                _t.sleep(30)
    else:
        print(f"  giving up after 5 attempts", flush=True)
        return
    print(f"  {len(urls)} files to fetch", flush=True)
    # Print a sample of URLs so we know what's there
    for u in list(urls["access_url"])[:5]:
        print(f"    sample: {u}", flush=True)
    # Filter: keep any URL referencing pbcor or _sci science data
    keep_urls = [u for u in urls["access_url"]
                 if "pbcor" in u or ("_sci" in u and ".fits" in u)]
    print(f"  filtered: {len(keep_urls)} pbcor/sci-fits files", flush=True)
    Alma.cache_location = str(OUTDIR)
    if keep_urls:
        files = Alma.download_files(keep_urls, savedir=str(OUTDIR),
                                       continuation=True)
        print(f"  downloaded: {len(files) if files else 0}")


if __name__ == "__main__":
    main()
