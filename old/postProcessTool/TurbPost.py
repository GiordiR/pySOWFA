
import numpy as np
import matplotlib.pyplot as plt
from PostProcessTool import *

####################################################################################
# TURBINE 0
turbName = 'turbine0'
turbineDir = './Steady_coarse_noNac'
turbineFileName = 'DTU10MW_POLIMI_WTM'

postProcDir = './Steady_coarse/postProcessing'


# probe reference wind streamwise velocity
probeName = 'probe_URef'
pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
pURef.readProbes(postProcDir)


# horizontal probes
for i in range(1, 9):
    probeName = 'probe_wake_h'+str(i)+'D'
    plotDir = './plot/' + probeName
    p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_h.readProbes(postProcDir)
    p_h.plotWakeProfile(i, plotDir, var='UMeanx')
    endPlot()


# vertical probes
for i in range(1, 9):
    probeName = 'probe_wake_v'+str(i)+'D'
    plotDir = './plot/' + probeName
    p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_v.readProbes(postProcDir)
    p_v.plotWakeProfile(i, plotDir, var='UMeanx')
    endPlot()

# experimental cross
for i in range(1, 3):
    probeName = 'probe_exp_cross'+str(i)
    plotDir = './plot/'+probeName
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir)
    p_c.plotWakeProfile(i, plotDir, var='UMeanx')
    p_ec = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_ec.readWakeExperiment(probeSet='cross')
    p_ec.plotWakeExperiment(i, plotDir, expProbe='probe_exp_cross'+str(i))
    endPlot()

# experimental along
for i in range(1, 3):
    probeName = 'probe_exp_along'+str(i)
    plotDir = './plot/'+probeName
    p_a = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_a.readProbes(postProcDir)
    p_a.plotWakeProfile(i, plotDir, var='UMeanx')
    p_ea = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_ea.readWakeExperiment(probeSet='along')
    p_ea.plotWakeExperiment(i, plotDir, expProbe='probe_exp_along'+str(i))
    endPlot()


# SOWFA turbine
turbineOutDir = './Steady_coarse_noNac/postProcessing/turbineOutput/0/'
plotDir = './plot/bladeSOWFA/'
bladeSowfa = SOWFA(turbName, turbineDir, turbineFileName)
bladeSowfa.readTurbineOutput(turbineOutDir, turbineNumber=1, nTurbines=1, FAST='on')
bladeSowfa.plotTurbine(1, 1, plotDir)
endPlot()

# FAST
outDir = './Steady_coarse_noNac/'
fastDir = outDir+'WTM'
outFileName = 'WTM.out'
plotDir = './plot/fastOut/'
fast = FAST(turbName, turbineDir, turbineFileName)
fast.readOutput(outDir, outFileName)
fast.readBladeProp(fastDir)

fast.plotBladeTipDeflections(1, plotDir, xlim=[13, 14])
endPlot()
fast.plotBladeRoot(2, plotDir, xlim=[13, 14])
endPlot()
fast.plotBladeOverTime(3, plotDir, ylim=True, xlim=[13, 14])
endPlot()
fast.plotBladeOverSpan(4, plotDir)
endPlot()
fast.plotRotor(5, plotDir, ylim=True)
endPlot()
fast.plotTowerBaseLoads(7, plotDir, ylim=True, xlim=[13, 14])
endPlot()
fast.plotTowerTopDisplacements(8, plotDir, xlim=[13, 14])
endPlot()
fast.plotLowSpeedShaft(9, plotDir, ylim=True)
endPlot()
fast.plotHighSpeedShaft(10, plotDir, ylim=True)
endPlot()

fast.plotBladeTipDeflectionsPSD(1, plotDir)
endPlot()
fast.plotBladeRootPSD(2, plotDir)
endPlot()
fast.plotRotorPSD(2, plotDir)
