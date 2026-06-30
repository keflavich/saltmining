"""Download G081.6802+00.5405A (DR21(OH) A) from 2019.1.00263.S (task #47).

Single ALMA program at 0.26" resolution. Field also contains DR21(OH) B
(separate target row in sources_L4_d2.csv); the cubes serve both since
the targets are <2 arcsec apart in the same MOUS pointing.
"""
import re
import warnings
from pathlib import Path

from astroquery.alma import Alma

warnings.filterwarnings("ignore")

UVDIR = Path("/orange/adamginsburg/salt/survey_2026/uvdata")
PROPS = ["2019.1.00263.S"]
DEST_NAME = "G081.6802+00.5405A"
NAMES = ["DR21(OH)", "DR21OH", "DR21_OH", "G081.682+0.541",
          "G081.682+00.541", "G81.683+00.541", "DR21OH", "DR_21_OH",
          "DR21(OH) North", "DR21", "DR21_DF1"]

SCI = re.compile(r"_sci\..*\.(cube|mfs|cont)(\.[a-z]+)?\.I\.pbcor\.fits$",
                  re.IGNORECASE)


def keep_url(url):
    return bool(SCI.search(url.rsplit("/", 1)[-1]))


def main():
    for prop in PROPS:
        print(f"\n=== {prop} ===", flush=True)
        q = Alma.query(payload={"project_code": prop}, public=True)
        if q is None or len(q) == 0:
            print("  no rows")
            continue
        df = q.to_pandas()
        mask = df["target_name"].astype(str).str.contains(
            "|".join(re.escape(n) for n in NAMES), case=False, regex=True)
        sub = df[mask]
        if len(sub) == 0:
            print(f"  target not found in {prop}; targets:")
            for t in sorted(set(df["target_name"])):
                print(f"    {t!r}")
            continue
        mous_list = sorted(set(sub["member_ous_uid"]))
        print(f"  {len(mous_list)} MOUS: {mous_list}", flush=True)
        dest = UVDIR / prop / DEST_NAME
        dest.mkdir(parents=True, exist_ok=True)
        for i, mous in enumerate(mous_list, 1):
            print(f"\n  [{i}/{len(mous_list)}] {mous}", flush=True)
            try:
                info = Alma.get_data_info(mous, expand_tarfiles=True)
            except (OSError, ValueError) as e:
                print(f"    ! get_data_info fail: {e}")
                continue
            if info is None or len(info) == 0:
                continue
            urls = list(info["access_url"])
            keep = [u for u in urls if keep_url(u)]
            to_get = [u for u in keep
                      if not (dest / u.rsplit("/", 1)[-1]).exists()]
            print(f"    {len(urls)} files, {len(keep)} sci pbcor, "
                  f"{len(to_get)} to fetch")
            if not to_get:
                continue
            Alma.cache_location = str(dest)
            files = Alma.download_files(to_get, savedir=str(dest),
                                          continuation=True)
            print(f"    downloaded {len(files) if files else 0}")


if __name__ == "__main__":
    main()
