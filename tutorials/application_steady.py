import os

from pySOWFA import *

def steady_flex():
    turbName = 'turbine0'
    turbineDir = '../../WTMresults/WTM_TURBULENT_fine_Rated_rigid'
    turbineFileName = 'DTU10MW_POLIMI_WTM'

    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readWakeExperiment()

