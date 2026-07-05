"""Run myXCLASS over 215-260 GHz with the comprehensive molfit."""
import os, sys
from pathlib import Path

XCLASSRootDir = os.environ.get("XCLASSRootDir",
                                "/orange/adamginsburg/software/XCLASS-Interface")
sys.path.insert(0, os.path.join(XCLASSRootDir, "build_tasks"))

from task_myXCLASS import myXCLASS

HERE = Path(__file__).parent
MOLFIT = HERE / "comprehensive_band6.molfit"

res = myXCLASS(
    FreqMin=215000.0, FreqMax=260000.0, FreqStep=2.0,
    TelescopeSize=1.5,
    BMIN=0.5, BMAJ=0.5, BPA=0.0,
    Inter_Flag=True, Redshift=0.0,
    t_back_flag=True, tBack=2.7, tslope=0.0,
    BackgroundFileName="",
    N_H=0.0, beta_dust=0.0, kappa_1300=0.0, DustFileName="",
    Te_ff=None, EM_ff=None,
    kappa_sync=None, B_sync=None, p_sync=None, l_sync=None,
    ContPhenFuncID=None,
    ContPhenFuncParam1=None, ContPhenFuncParam2=None,
    ContPhenFuncParam3=None, ContPhenFuncParam4=None,
    ContPhenFuncParam5=None,
    MolfitsFileName=str(MOLFIT),
    iso_flag=False, IsoTableFileName="",
    CollisionFileName="",
    NumModelPixelXX=100, NumModelPixelYY=100,
    LocalOverlapFlag=False, NoSubBeamFlag=True,
    dbFilename="",
    RestFreq=0.0, vLSR=0.0,
)
print(f"\nReturn objs: {len(res)}")
# Find JobDir by inspecting the most recent job
import glob
jobs = sorted(glob.glob("/orange/adamginsburg/software/XCLASS-Interface/run/myXCLASS/job__30-06-2026*"),
              key=os.path.getmtime)
if jobs:
    jd = jobs[-1]
    print(f"Most recent job dir: {jd}")
    tt = os.path.join(jd, "transition_energies.dat")
    if os.path.exists(tt):
        print(f"  transition table: {tt}  size={os.path.getsize(tt)}")
