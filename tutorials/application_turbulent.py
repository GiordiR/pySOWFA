from pySOWFA import *

def turbulent_turbine(turbName):

    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_' + turbName.lower()

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
    for i in [x for x in range(1,9) if x!=2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.readSets(postProcDir=postProcDir, var='UPrime2Meanxx')
        p_h.getTurbulenceIntensity()
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', xLim=[0.5, 4.5], filter=None)
        p_h.plotWakeProfile(i+7, plotDir, var='TIx', sampleType='sets')

    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.readSets(postProcDir=postProcDir, var='UPrime2Meanxx')
        p_v.getTurbulenceIntensity()
        p_v.plotWakeProfile(i+14, plotDir, var='UMeanx', sampleType='sets')
        p_v.plotWakeProfile(i+21, plotDir, var='TIx', sampleType='sets')
    endPlot()

    # SOWFA turbine
    turbineOutPath = os.path.join(postProcDir, 'turbineOutput')
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'SOWFA'))
    bladeSowfa = SOWFA(turbName, turbineDir, turbineFileName)
    bladeSowfa.readTurbineOutput(turbineOutPath, turbineNumber=1, nTurbines=1, fast='on')
    bladeSowfa.plotTurbine(1, 1, plotDir)
    endPlot()

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    outFileName = 'WTM.out'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFAST'))
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
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

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics.txt'))


def turbulent_turbine_comparison():
    turbName = 'Rigid'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_rigid'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1,9) if x!=2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic')
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2])

    ##################################################################
    turbName = 'Flex'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)

    ##################################################################
    turbName = 'Flex1'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex1'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)

    ##################################################################
    turbName = 'Flex2'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)

    ##################################################################
    turbName = 'Flex3'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex3'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)

    ##################################################################
    turbName = 'Flex4'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex4'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)


    ##################################################################
    turbName = 'Flex5'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex5'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2, 1.2],
                            interp='quadratic', compareID=i)
    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)


def steady_mesh_comparison():
    filter=[15,3]
    xLim=[0.2, 1.2]

    turbName = 'Fine mesh'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_STEADY_fine_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')
    expDir = os.path.abspath(os.path.join(path, os.pardir, 'UNAFLOW', 'WAKE'))

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # experimental cross
    probeName = 'probe_exp_cross1'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir=postProcDir)
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='probe', normVar=[pURef, 'UMeanx'], interp='linear', filter=filter, xLim=xLim)
    p_ec = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_ec.readWakeExperiment(expDir=expDir, probeSet='cross')
    p_ec.plotWakeExperiment(1, plotDir, expProbe='probe_exp_cross', norm=True, interp='linear', filter=[15,3])
    p_ec.plotWakeError(2, plotDir, probeRef=[p_ec, 'UcrossMeanx', 'yCross'], probeCompare=[p_c, 'UMeanx'], labelCompare=turbName)

    ##################################################################
    turbName = 'Medium mesh'
    turbineDirName = 'WTM_STEADY_medium_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readProbes(postProcDir=postProcDir)

    # experimental cross
    probeName = 'probe_exp_cross1'
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir=postProcDir)
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='probe', normVar=[pURef, 'UMeanx'], compareID=1, interp='linear', filter=filter, xLim=xLim)
    p_ec.plotWakeError(2, plotDir, probeRef=[p_ec, 'UcrossMeanx', 'yCross'], probeCompare=[p_c, 'UMeanx'],labelCompare=turbName, compareID=2)

    ##################################################################
    turbName = 'Coarse mesh'
    turbineDirName = 'WTM_STEADY_coarse_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readProbes(postProcDir=postProcDir)

    # experimental cross
    probeName = 'probe_exp_cross1'
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir=postProcDir)
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='probe', normVar=[pURef, 'UMeanx'], compareID=1, interp='linear', filter=filter, xLim=xLim)
    p_ec.plotWakeError(2, plotDir, probeRef=[p_ec, 'UcrossMeanx', 'yCross'], probeCompare=[p_c, 'UMeanx'],labelCompare=turbName, compareID=2)


def steady_epsilon_comparison():
    filter = None
    xLim = [0.2, 1.2]

    turbName = r'\epsilon = 2 \Delta_g'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_STEADY_fine_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')
    expDir = os.path.abspath(os.path.join(path, os.pardir, 'UNAFLOW', 'WAKE'))

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # experimental cross
    probeName = 'probe_wake_cross1'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readSets(postProcDir=postProcDir, var='UMeanx')
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], filter=filter,
                        xLim=xLim)

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    outFileName = 'WTM.out'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFAST'))
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatisticsE2.txt'))

    ##################################################################
    turbName = r'\epsilon = 1.25 \Delta_g'
    turbineDirName = 'WTM_STEADY_fine_flex1.25'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # experimental cross
    probeName = 'probe_wake_cross1'
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readSets(postProcDir=postProcDir, var='UMeanx')
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], compareID=1,
                        filter=filter, xLim=xLim)

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    outFileName = 'WTM.out'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFAST'))
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatisticsE125.txt'))


def test_interpolate_wakes(turbType):
    turbName = 'Normal'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_' + turbType

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # horizontal probes
    for i in [x for x in range(1,9) if x!=2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2,1.2])

    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', normVar=[pURef, 'UMeanx'], sampleType='sets', xLim=[0.2,1.2])

    ##################################################################
    turbName = 'Interp'

    # horizontal probes
    for i in [x for x in range(1, 9) if x != 2]:
        probeName = 'probe_wake_h' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_h = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_h.readSets(postProcDir=postProcDir, var='UMeanx')
        p_h.plotWakeProfile(i, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], xLim=[0.2,1.2], interp='quadratic', compareID=i)

    # vertical probes
    for i in range(1, 9):
        probeName = 'probe_wake_v' + str(i) + 'D'
        plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
        p_v = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
        p_v.readSets(postProcDir=postProcDir, var='UMeanx')
        p_v.plotWakeProfile(i+10, plotDir, var='UMeanx', sampleType='sets', normVar=[pURef, 'UMeanx'], interp='linear', filter=[15,5], xLim=[0.2,1.2], compareID=i+10)

    ##################################################################


def test_interpolate_exp(turbType):
    filter = [11, 5]
    xLim = [0.2, 1.2]

    turbName = turbType
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_STEADY_' + turbType + '_flex2'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # probe reference wind streamwise velocity
    probeName = 'probe_URef'
    pURef = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    pURef.readSets(postProcDir=postProcDir, var='UMeanx')

    # experimental cross
    probeName = 'probe_exp_cross1'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir=postProcDir)
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='probe', normVar=[pURef, 'UMeanx'], xLim=xLim)

    ##################################################################
    turbName = 'Interp'

    # experimental cross
    probeName = 'probe_exp_cross1'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFOAM'))
    p_c = OpenFOAM(turbName, probeName, turbineDir, turbineFileName)
    p_c.readProbes(postProcDir=postProcDir)
    p_c.plotWakeProfile(1, plotDir, var='UMeanx', sampleType='probe', normVar=[pURef, 'UMeanx'], interp='linear',
                        filter=filter, xLim=xLim, compareID=1)
    ##################################################################


def generate_statistics():
    turbName = 'Rigid'
    turbineFileName = 'DTU10MW_POLIMI_WTM'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_rigid'

    # Set directories paths
    path = os.getcwd()
    turbineDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    postProcDir = os.path.join(turbineDir, 'postProcessing')

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    outFileName = 'WTM.out'
    plotDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', 'plots', 'OpenFAST'))
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex1'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex1'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex2'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex2'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex3'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex3'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex4'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex4'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

    ##################################################################
    turbName = 'Flex5'
    turbineDirName = 'WTM_TURBULENT_fine_Rated_flex5'

    # FAST
    fastOutDir = os.path.abspath(os.path.join(path, os.pardir, 'WTMresults', turbineDirName))
    fastDir = os.path.join(fastOutDir, 'WTM')
    fast = FAST(turbName, turbineDir, turbineFileName)
    fast.readOutput(fastOutDir, outFileName)
    fast.readBladeProp(fastDir)

    fast.generateStatistics(statFile=os.path.join(plotDir, 'FASTstatistics'+turbName+'.txt'))

