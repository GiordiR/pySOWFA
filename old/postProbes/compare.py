
import numpy as np
import matplotlib.pyplot as plt
from PostProcessTool import *

####################################################################################
# STEADY COARSE
turbName = 'Rigid'
turbineDir = './STEADY_medium_rigid2'
turbineFileName = 'DTU10MW_POLIMI_WTM'
postProcDir = './STEADY_medium_rigid2/postProcessing'

'''
# horizontal probes
for i in range(1, 9):
    probeName = 'probe_wake_h'+str(i)+'D'
    plotDir = './plot_compare/' + probeName
    p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_h.readProbes(postProcDir)
    p_h.plotWakeProfile(i, plotDir, var='UMeanx')

# vertical probes
for i in range(1, 9):
    probeName = 'probe_wake_v'+str(i)+'D'
    plotDir = './plot_compare/' + probeName
    p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_v.readProbes(postProcDir)
    p_v.plotWakeProfile(i+8, plotDir, var='UMeanx')
'''
# experimental cross
for i in range(1, 3):
    probeName = 'probe_exp_cross'+str(i)
    plotDir = './plot_compare/'+probeName
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir)
    p_c.plotWakeProfile(i+16, plotDir, var='UMeanx')
    p_ec = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_ec.readWakeExperiment(probeSet='cross')
    p_ec.plotWakeExperiment(i+16, plotDir, expProbe='probe_exp_cross'+str(i))
    p_ec.plotWakeError(i+18, plotDir, probeRef=[p_ec, 'xVarWakeExp', 'yVarWakeExp'], probeCompare=[p_c, 'xVarWake'], labelCompare='Rigid')
'''
# SOWFA turbine
turbineOutDir = './STEADY_medium_rigid2/postProcessing/turbineOutput/0/'
plotDir = './plot/SOWFA/'
bladeSowfa_c = SOWFA(turbName, turbineDir, turbineFileName)
bladeSowfa_c.readTurbineOutput(turbineOutDir, turbineNumber=1, nTurbines=1, FAST='on')
bladeSowfa_c.plotTurbine(1, 1, plotDir, var='rotorAxialForce', ylim=[36, 40])
bladeSowfa_c.plotTurbine(2, 1, plotDir, var='rotorHorizontalForce')
bladeSowfa_c.plotTurbine(3, 1, plotDir, var='rotorPower', ylim=[80, 100])
bladeSowfa_c.plotTurbine(4, 1, plotDir, var='rotorSpeed')
bladeSowfa_c.plotTurbine(5, 1, plotDir, var='rotorTorque', ylim=[3,4])
bladeSowfa_c.plotTurbine(6, 1, plotDir, var='rotorVerticalForce')


# FAST
outDir = './STEADY_medium_rigid2/'
fastDir = outDir+'WTM'
outFileName = 'WTM.out'
plotDir = './plot/fastOut/'
fast = FAST(turbName, turbineDir, turbineFileName)
fast.readOutput(outDir, outFileName)
fast.readBladeProp(fastDir)
#fast.generateStatistics()
'''
#####################################################################################
# STEADY MEDIUM
turbName = 'Elastic'
turbineDir = './STEADY_medium_flex2'
turbineFileName = 'DTU10MW_POLIMI_WTM'
postProcDir = './STEADY_medium_flex2/postProcessing'

'''
# horizontal probes
for i in range(1, 9):
    probeName = 'probe_wake_h'+str(i)+'D'
    plotDir = './plot_compare/' + probeName
    p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_h.readProbes(postProcDir)
    p_h.plotWakeProfile(i, plotDir, var='UMeanx', compareID=i)
    
# vertical probes
for i in range(1, 9):
    probeName = 'probe_wake_v'+str(i)+'D'
    plotDir = './plot_compare/' + probeName
    p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_v.readProbes(postProcDir)
    p_v.plotWakeProfile(i+8, plotDir, var='UMeanx', compareID=i+8)
'''
# experimental cross
for i in range(1, 3):
    probeName = 'probe_exp_cross'+str(i)
    plotDir = './plot_compare/'+probeName
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir)
    p_c.plotWakeProfile(i+16, plotDir, var='UMeanx', compareID=i+16)
    p_ec.readWakeExperiment(probeSet='cross')
    p_ec.plotWakeExperiment(i + 16, plotDir, expProbe='probe_exp_cross' + str(i))
    p_c.plotWakeError(i + 18, plotDir, probeRef=[p_ec, 'xVarWakeExp', 'yVarWakeExp'], probeCompare=[p_c, 'xVarWake'], labelCompare='Elastic', compareID=i+18)
'''
# SOWFA turbine
turbineOutDir = './STEADY_medium_flex2/postProcessing/turbineOutput/0/'
plotDir = './plot/SOWFA/'
bladeSowfa_f = SOWFA(turbName, turbineDir, turbineFileName)
bladeSowfa_f.readTurbineOutput(turbineOutDir, turbineNumber=1, nTurbines=1, FAST='on')
bladeSowfa_f.plotTurbine(1, 1, plotDir, var='rotorAxialForce', ylim=[36, 40], compareID=1)
bladeSowfa_f.plotTurbine(2, 1, plotDir, var='rotorHorizontalForce', compareID=2)
bladeSowfa_f.plotTurbine(3, 1, plotDir, var='rotorPower', ylim=[80, 100], compareID=3)
bladeSowfa_f.plotTurbine(4, 1, plotDir, var='rotorSpeed', compareID=4)
bladeSowfa_f.plotTurbine(5, 1, plotDir, var='rotorTorque', ylim=[3,4], compareID=5)
bladeSowfa_f.plotTurbine(6, 1, plotDir, var='rotorVerticalForce', compareID=6)

# FAST
outDir = './STEADY_medium_flex2/'
fastDir = outDir+'WTM'
outFileName = 'WTM.out'
plotDir = './plot/fastOut/'
fast = FAST(turbName, turbineDir, turbineFileName)
fast.readOutput(outDir, outFileName)
fast.readBladeProp(fastDir)
fast.generateStatistics()
'''
################################################################################################à
# MEDIUM SOWFA
turbName = 'SOWFA'
turbineDir = './STEADY_medium_SOWFA'
turbineFileName = 'DTU10MW_POLIMI_WTM'
postProcDir = './STEADY_medium_SOWFA/postProcessing'

# experimental cross
for i in range(1, 3):
    probeName = 'probe_exp_cross'+str(i)
    plotDir = './plot_compare/'+probeName
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir)
    p_c.plotWakeProfile(i+16, plotDir, var='UMeanx', compareID=i+16)
    p_ec.readWakeExperiment(probeSet='cross')
    p_ec.plotWakeExperiment(i + 16, plotDir, expProbe='probe_exp_cross' + str(i))
    p_c.plotWakeError(i + 18, plotDir, probeRef=[p_ec, 'xVarWakeExp', 'yVarWakeExp'], probeCompare=[p_c, 'xVarWake'], labelCompare='SOWFA', compareID=i+18)

endPlot()