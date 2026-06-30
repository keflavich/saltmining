"""Download IRAS17233-3606 (task #39).

7 proposals available with <500 AU resolution at d=1.3 kpc:
  2013.1.00260.S  0.067"
  2015.1.00496.S  0.022"  (best resolution, likely largest)
  2016.1.01137.S  0.110"
  2017.1.00237.S  0.049"
  2021.1.00095.S  0.274"
  2021.1.00311.S  0.269"
  2023.1.01382.S  0.160"

Choose a moderate-resolution-and-size proposal to start; rerun with
expanded PROPS list once disk usage is verified.
"""
import re
import warnings
from pathlib import Path

from astroquery.alma import Alma

warnings.filterwarnings("ignore")

UVDIR = Path("/orange/adamginsburg/salt/survey_2026/uvdata")
PROPS = ["2017.1.00237.S", "2016.1.01137.S"]
DEST_NAME = "IRAS17233-3606"
NAMES = ["IRAS17233-3606", "IRAS_17233-3606", "IRAS_17233_3606",
         "I17233", "IRAS17233", "G351.7745-00.5377", "G351.77-0.54"]

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
