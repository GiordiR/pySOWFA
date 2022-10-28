import os

from pySOWFA import *

def turbulent_turbine():
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
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in range(1, 9):
        if not i == 2:
            probeName = 'probe_wake_h' + str(i) + 'D'
            plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
            p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
            p_h.readSets(postProcDir=postProcDir, var='UMeanx')
            #p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'])
            p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', xLim=[0.5, 4.5], filter=None)
            endPlot()

    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets')
        endPlot()

    # experimental cross
    probeName = 'probe_wake_cross1'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readSets(postProcDir=postProcDir, var='UMeanx')
    p_c.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets')
    p_ec = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_ec.readWakeExperiment(expDir=expDir, probeSet='cross')
    p_ec.plotWakeExperiment(i, plotDir, expProbe='probe_exp_cross' + str(i))
    endPlot()
    """
    # SOWFA turbine
    turbineOutDir = './Steady_coarse_noNac/postProcessing/turbineOutput/0/'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'SOWFA'))
    bladeSowfa = SOWFA(turbName, turbineDir, turbineFileName)
    bladeSowfa.readTurbineOutput(turbineOutDir, turbineNumber=1, nTurbines=1, FAST='on')
    bladeSowfa.plotTurbine(1, 1, plotDir)
    endPlot()

    # FAST
    outDir = './Steady_coarse_noNac/'
    fastDir = outDir + 'WTM'
    outFileName = 'WTM.out'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFAST'))
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
    """
