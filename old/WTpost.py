from PostProcessTool import *
import gc
####################################################################################

# WIND TUNNEL
turbName='windTunnel - Above Rated'
turbineDir = './WTA/'
postProcDir = turbineDir + 'postProcessing/'
'''
# Longitudinal probe
probeName = 'probe_long'
plotDir = './plot_windTunnel/' + probeName
plong = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
plong.readProbes(postProcDir)
plong.plotSpaceCorrelations(plotDir)

'''
# Spectrum 0D probe
probeName = 'probe_spectrum0D'
plotDir = './plot_windTunnel/' + probeName
pspec0 = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
pspec0.readProbes(postProcDir)
pspec0.plotEnergySpectrum(plotDir)
pspec0.plotTimeCorrelations(plotDir)

'''
# Spectrum 1D front probe
probeName = 'probe_spectrum1DFront'
plotDir = './plot_windTunnel/' + probeName
pspec1F = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
pspec1F.readProbes(postProcDir)
pspec1F.plotEnergySpectrum(plotDir)
pspec1F.plotTimeCorrelations(plotDir)


# Spectrum 1D back probe
probeName = 'probe_spectrum1DBack'
plotDir = './plot_windTunnel/' + probeName
pspec1B = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
pspec1B.readProbes(postProcDir)
pspec1B.plotEnergySpectrum(plotDir)
pspec1B.plotTimeCorrelations(plotDir)


# Transversal probe
probeName = 'probe_trasv'
plotDir = './plot_windTunnel/' + probeName
ptrasv = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
ptrasv.readProbes(postProcDir)
ptrasv.plotSpaceCorrelations(plotDir)


'''
# probe reference wind streamwise velocity
probeName = 'probe_URef'
pURef = OpenFOAM(turbName, probeName, turbineDir)
pURef.readProbes(postProcDir)

# Vertical 0R probe
probeName = 'probe_vert0R'
plotDir = './plot_windTunnel/' + probeName
pvert0 = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
pvert0.readProbes(postProcDir)
pvert0.plotWindProfile(1, plotDir, var='UMeanx', label='y = 0R')
pvert0.getTurbulenceIntesity()
pvert0.plotWindProfile(2, plotDir, var='TIx', label='y = 0R')
pvert0.plotWindProfile(3, plotDir, var='UPrime2Meanxz')
pvert0.plotWindProfile(4, plotDir, var='UPrime2Meanxx')
pvert0.plotWindProfile(5, plotDir, var='UPrime2Meanyy')
pvert0.plotWindProfile(6, plotDir, var='UPrime2Meanzz')
pvert0.plotResiduals(plotDir, var='UMeanx')
pvert0.plotResiduals(plotDir, var='TIx')
pvert0.plotResiduals(plotDir, var='UPrime2Meanxz')


# Vertical 1R probe
probeName = 'probe_vert1R'
plotDir = './plot_windTunnel/' + probeName
pvert1 = OpenFOAM(turbName, probeName, turbineDir=turbineDir)
pvert1.readProbes(postProcDir)
pvert1.plotWindProfile(1, plotDir, var='UMeanx', label='y = 1R', compareID=1)
pvert1.getTurbulenceIntesity()
pvert1.plotWindProfile(2, plotDir, var='TIx', label='y = 1R', compareID=2)
endPlot()

