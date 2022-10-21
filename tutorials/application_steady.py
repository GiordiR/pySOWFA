import os

from pySOWFA import *

def steady_flex():
    turbName = 'turbine0'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_rigid'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')
    expDir = os.path.abspath(os.path.join(path, os.pardir, 'UNAFLOW', 'WAKE'))

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readWakeExperiment(expDir=expDir)

    pURef.readSets(postProcDir=postProcDir, var='Ux')
