# Salt-survey ALMA imaged-cube coverage

Per-source summary of which pipeline-imaged FITS cubes have been
downloaded under `uvdata/<proposal>/<source>/` and which salt-relevant
lines are covered by their spw frequency footprints.

Salt lines checked (rest GHz):

- NaCl_J7: 91.1697
- NaCl_J8: 104.1896
- NaCl_J10: 130.2234
- NaCl_J14: 182.2663
- NaCl_J17: 221.2727
- NaCl_J18: 234.2519
- NaCl_J19: 247.2278
- NaCl_J20: 260.2000
- NaCl_J22: 286.1304
- NaCl_J25: 324.9886
- NaCl_J27: 350.8682
- NaCl_v1_18-17: 232.5100
- NaCl_v1_17-16: 219.6149
- NaCl_v2_17-16: 217.9800
- KCl_J12: 92.1008
- KCl_J13: 99.7738
- KCl_J20: 153.4294
- KCl_J28: 215.0083
- KCl_J30: 230.3206
- KCl_J32: 245.6255
- KCl_J34: 260.9162
- KCl_J38: 291.4603
- KCl_J45: 344.7888
- KClv3_29: 218.5797
- KClv3_31: 233.6057
- H2O_232: 232.6867
- H2O_325: 325.1530
- H2O_183: 183.3101
- H2O_336: 336.2280
- H2O_321: 321.2256
- H2O_437: 437.3467
- SiO_2-1: 86.8470
- SiO_5-4: 217.1050
- SiO_8-7: 347.3306

## Per-source coverage

| Source | Proposal | n_cubes | freq range (GHz) | covered lines | needs reimaging | notes |
|---|---|---:|---|---|:-:|---|
| G018.3412+01.7681 | 2016.1.01036.S | 24 | 216.86-235.44 | NaCl_J18, NaCl_v1_18-17, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 24 cubes; 9 unique spws |
| G336.4917-01.4741A | 2016.1.01182.S | 4 | 218.00-235.91 | NaCl_J18, NaCl_v1_18-17, NaCl_v1_17-16, KClv3_29, KClv3_31, H2O_232 | no | 4 cubes; 4 unique spws |
| G336.4917-01.4741B | 2016.1.01182.S | 4 | 218.00-235.91 | NaCl_J18, NaCl_v1_18-17, NaCl_v1_17-16, KClv3_29, KClv3_31, H2O_232 | no | 4 cubes; 4 unique spws |
| G035.1979-00.7427 | 2017.1.00181.S | 9 | 216.50-234.83 | NaCl_J18, KClv3_31, SiO_5-4 | no | 9 cubes; 9 unique spws |
| G014.9958-00.6732 | 2019.1.00195.L | 600 | 216.96-220.91 | NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, SiO_5-4 | no | 600 cubes; 8 unique spws |
| G015.1288-00.6717 | 2019.1.00195.L | 600 | 216.96-220.91 | NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, SiO_5-4 | no | 600 cubes; 8 unique spws |
| G339.6221-00.1209 | 2019.1.00195.L | 844 | 217.02-220.97 | NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, SiO_5-4 | no | 844 cubes; 8 unique spws |
| G081.6802+00.5405B | 2019.1.00263.S | 12 | 216.76-233.94 | NaCl_v1_18-17, NaCl_v1_17-16, NaCl_v2_17-16, KCl_J30, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 12 cubes; 8 unique spws |
| G192.6005-00.0479 | 2019.1.00315.S | 4 | 334.62-349.53 | H2O_336, SiO_8-7 | no | 4 cubes; 4 unique spws |
| G339.8838-01.2588 | 2019.1.00517.S | 1 | 259.86-260.80 | NaCl_J20 | no | 1 cubes; 1 unique spws |
| G026.3819+01.4057A | 2021.1.00095.S | 108 | 216.97-234.43 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 108 cubes; 7 unique spws |
| G344.2207-00.5953 | 2021.1.00095.S | 76 | 216.97-234.43 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 76 cubes; 5 unique spws |
| G015.0357-00.6795 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G316.8112-00.0566 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G328.2523-00.5320A | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G328.2523-00.5320B | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G328.8074+00.6324 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G343.1261-00.0623 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G345.4881+00.3148 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G345.4938+01.4677 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G348.5312-00.9714 | 2022.1.00700.S | 1 | 231.00-232.87 | NaCl_v1_18-17, H2O_232 | no | 1 cubes; 1 unique spws |
| G326.6618+00.5207 | 2022.1.01344.S | 21 | 217.11-234.48 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 21 cubes; 7 unique spws |
| G345.0034-00.2240A | 2022.1.01344.S | 14 | 217.08-234.45 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 14 cubes; 7 unique spws |
| G345.5043+00.3480 | 2022.1.01344.S | 14 | 217.08-234.45 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 14 cubes; 7 unique spws |
| G348.7250-01.0435 | 2022.1.01344.S | 14 | 217.07-234.44 | NaCl_J18, NaCl_v1_17-16, NaCl_v2_17-16, KClv3_29, KClv3_31, H2O_232, SiO_5-4 | no | 14 cubes; 7 unique spws |
| G012.8062-00.1987 | 2023.1.01346.S | 9 | 214.74-260.54 | NaCl_J20, NaCl_v1_18-17, NaCl_v2_17-16, KCl_J28, KCl_J32, H2O_232, SiO_5-4 | no | 9 cubes; 9 unique spws |
| G018.3029-00.3910 | 2013.1.01193.S | 0 | - | - | yes | no cubes on disk |
| G017.6380+00.1566 | 2017.1.00098.S | 0 | - | - | yes | no cubes on disk |
| G012.9090-00.2607 | 2018.1.00458.S | 0 | - | - | yes | no cubes on disk |
| G076.3829-00.6210 | 2022.1.00367.S | 0 | - | - | yes | no cubes on disk |
| G010.8411-02.5919 | 2022.1.00386.S | 0 | - | - | yes | no cubes on disk |
| G026.4207+01.6858 | 2022.1.01160.S | 4 | 238.76-257.71 | - | yes | 4 cubes; 4 unique spws; no salt lines in covered spws |
| G189.0307+00.7821 | 2022.1.01160.S | 0 | - | - | yes | no cubes on disk |
| G232.6207+00.9959 | 2023.1.00394.S | 0 | - | - | yes | no cubes on disk |
| G345.0061+01.7944C | 2023.1.01346.S | 0 | - | - | yes | no cubes on disk |
| G348.6972-01.0263 | 2023.1.01346.S | 0 | - | - | yes | no cubes on disk |
| G080.8645+00.4197 | 2024.1.01642.S | 0 | - | - | yes | no cubes on disk |
| G345.2244+01.0304 | 2024.1.01709.S | 0 | - | - | yes | no cubes on disk |

**Summary:** 26/38 sources have at least one salt line covered by an existing pipeline cube; the remaining 12/38 need reimaging from raw uv data (use `download_alma_data.py --include-asdm` then `imaging/run_imaging.py`).
