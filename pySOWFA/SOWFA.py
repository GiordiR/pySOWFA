import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as scsig
import scipy as sc

import pySOWFA.Turbine as Turbine
from pySOWFA.plotUtils import plot, plotUtils, plotCompare, plotLog, plotLogCompare, plotLogUtils, plotLogLog, \
                      plotLogLogCompare, plot3D, endPlot, getAxes, getTitle


class SOWFA(Turbine):
    """
    Class to postprocess SOWFA related output files (e.g. Turbine output)

    :param str turbineName: turbine name [turbine name, windTunnel, precursor, noTurbine]
    :param str turbineDir: turbine directory path
    :param str turbineFileName: turbine file name
    """

    def __init__(self, turbineName, turbineDir=None, turbineFileName=None):
        Turbine.__init__(self, turbineName, turbineDir, turbineFileName)

    def readTurbineOutput(self, turbineOutDir=None, turbineNumber=1, nTurbines=1, fast=False):
        """
        Read SOWFA turbine output files

        :param str turbineOutDir: Wind Turbine output directory path
        :param int turbineNumber: Wind Turbine ID number
        :param int nTurbines: number of Wind Turbines considered
        :param boolean fast: change output selection if SOWFA coupled with OpenFAST
        """
        if turbineOutDir is None:
            turbineOutDir = "./postProcessing/turbineOutput/0/"

        # Find files in the directory
        files = os.listdir(turbineOutDir)
        self.SOWFAturbine = []
        self.SOWFAnacelle = []
        self.SOWFAtower = []
        self.SOWFAblade = []

        # Only rotor files for SOWFA-FAST
        if fast:
            for a in files[:]:
                if a.startswith(('blade', 'nacelle', 'tower')):
                    files.remove(a)

        # Read turbine output files
        for file in files:
            turb_db = pd.read_csv(turbineOutDir + file, sep=' ', skiprows=1, header=None)
            turb_db.dropna(how='all', axis=1, inplace=True)
            values = turb_db.values
            for i in range(0, nTurbines):
                if file.startswith(('rotor', 'generator', 'bladePitch')):
                    vars(self)[file + str(turbineNumber)] = values[i::nTurbines, -1]
                    self.turbineTime = values[i::nTurbines, 1]
                    self.SOWFAturbine.append(file)
                elif file.startswith('nacelle'):
                    # Collect global nacelle values and first point nacelle values
                    vars(self)[file + str(turbineNumber)] = values[i::nTurbines, 3]
                    self.SOWFAnacelle.append(file)
                elif file.startswith('tower'):
                    # Collect global tower values and tower top values
                    vars(self)[file + str(turbineNumber)] = values[i::nTurbines, -1]
                    self.SOWFAtower.append(file)
                elif file.startswith('blade') and file != 'bladePitch':
                    vars(self)[file + 'Time'] = values[0::3 * nTurbines, 2]
                    vars(self)[file + 'dt'] = values[0, 3]
                    vars(self)[file + str(turbineNumber)] = values[i::3 * nTurbines, 4:]
                    vars(self)[file + str(turbineNumber) + 'Root'] = values[i::3 * nTurbines, 4]
                    vars(self)[file + str(turbineNumber) + 'Tip'] = values[i::3 * nTurbines, -1]
                    vars(self)[file + str(turbineNumber) + 'Ft'] = values[-nTurbines + i, 4:]
                    blade_span = np.linspace(self.blade_r[0], self.blade_r[-1], num=self.numBladePoints, dtype=float)
                    self.blade_span_norm = blade_span / (self.blade_r[-1])
                    self.SOWFAblade.append(file)

    def plotTurbine(self, figID, turbineNumber=1, plotDir=None, var='all', compareID=None):
        """
        Plot Wind Turbine performance

        :param int figID: figure identification number
        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAturbine:
                xVar = self.turbineTime
                yVar = vars(self)[file + str(turbineNumber)]
                xlabel = 'Time [s]'
                ylabel = getAxes(file)
                label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
                title = getTitle(file)
                figName = file + '_turbine'
                if compareID is None:
                    plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
                else:
                    plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                figID += 1
        else:
            xVar = self.turbineTime
            yVar = vars(self)[var + str(turbineNumber)]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
            title = getTitle(var)
            figName = var + '_turbine'
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotNacelle(self, figID, turbineNumber=1, plotDir=None, var='all', compareID=None):
        """
        Plot Nacelle performance

        :param int figID: figure identification number
        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAnacelle:
                xVar = self.turbineTime
                yVar = vars(self)[file + str(turbineNumber)]
                xlabel = 'Time [s]'
                ylabel = getAxes(file)
                label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.nacelleepsilon0)
                title = getTitle(file)
                figName = file + '_nacelle'
                if compareID is None:
                    plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
                else:
                    plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                figID += 1
        else:
            xVar = self.turbineTime
            yVar = vars(self)[var + str(turbineNumber)]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.nacelleepsilon0)
            title = getTitle(var)
            figName = var + '_nacelle'
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTower(self, figID, turbineNumber=1, plotDir=None, var='all', compareID=None):
        """
        Plot Tower performance

        :param int figID: figure identification number
        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAtower:
                xVar = self.turbineTime
                yVar = vars(self)[file + str(turbineNumber)]
                xlabel = 'Time [s]'
                ylabel = getAxes(file)
                label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
                title = getTitle(file)
                figName = file + '_tower'
                if compareID is None:
                    plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
                else:
                    plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                figID += 1
        else:
            xVar = self.turbineTime
            yVar = vars(self)[var + str(turbineNumber)]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
            title = getTitle(var)
            figName = var + '_tower'
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeOverSpan(self, figID, turbineNumber=1, plotDir=None, var='all', compareID=None):
        """
        Plot Wind Turbine blade performance over blade span

        :param int figID: figure identification number
        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAblade:
                if 'blade' in file:
                    xVar = self.blade_span_norm
                    yVar = vars(self)[file + str(turbineNumber) + 'Ft']
                    xlabel = 'r/R'
                    ylabel = getAxes(file)
                    label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
                    title = getTitle(file) + ' Final Time'
                    figName = file + '_bladeFt'
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.blade_span_norm
            yVar = vars(self)[var + str(turbineNumber) + 'Ft']
            xlabel = 'r/R'
            ylabel = getAxes(var)
            label = self.bladeForceProjectionType + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
            title = getTitle(var) + ' Final Time'
            figName = var + '_bladeFt'
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeOverTime(self, figID, turbineNumber=1, plotDir=None, var='all', Root_Tip='Root', compareID=None):
        """
        Plot Wind Turbine blade root and tip performance over time

        :param int figID: figure identification number
        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param str Root_Tip: blade station to be plotted ['Root', 'Tip']
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAblade:
                if 'blade' in file:
                    xVar = vars(self)[file + 'Time']
                    yVar = vars(self)[file + str(turbineNumber) + Root_Tip]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(file)
                    label = Root_Tip + self.bladeForceProjectionType \
                            + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
                    title = getTitle(file)
                    figName = file + '_blade_' + Root_Tip
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = vars(self)[var + 'Time']
            yVar = vars(self)[var + str(turbineNumber) + Root_Tip]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = Root_Tip + self.bladeForceProjectionType \
                    + ' - ' + chr(949) + '=' + str(self.bladeepsilon0)
            title = getTitle(var)
            figName = var + '_blade_' + Root_Tip
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeFrequency(self, turbineNumber=1, plotDir=None, var='all', point='Tip'):
        """
        Plot Wind Turbine blade performance in frequency domain

        :param int turbineNumber: Wind Turbine ID number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param str point: blade station to be plotted ['Root', 'Tip']
        """
        if plotDir is None:
            plotDir = './postProcessing/turbineOutput/plots/'
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for file in self.SOWFAblade:
                if 'blade' in file:
                    fSampl = 1.0 / vars(self)[file + 'dt']
                    totTime = vars(self)[file + 'Time'][-1] - vars(self)[file + 'Time'][0]
                    dimSampl = np.shape(vars(self)[file + str(turbineNumber) + point])[0]
                    window = sc.signal.hann(dimSampl)
                    valWindow = vars(self)[file + str(turbineNumber) + point] * window
                    df = 1.0 / (dimSampl * 1 / fSampl)
                    freq = np.arange(0, dimSampl / 2 * df, df)
                    specPos = freq * 2.0
                    specTot = sc.fft(valWindow)
                    specPos[0] = specTot[0] / dimSampl
                    if dimSampl % 2 == 0:
                        specPos[1:round(dimSampl / 2) - 1] = specTot[1:round(dimSampl / 2) - 1] / [dimSampl / 2]
                        specPos[round(dimSampl / 2) - 1] = specTot[round(dimSampl / 2)] / dimSampl
                    else:
                        specPos[1:round(dimSampl / 2)] = specTot[1:round(dimSampl / 2)] / [dimSampl / 2]
                    fig, (ax1, ax2) = plt.subplots(1, 2)
                    fig.suptitle(file)
                    ax1.plot(vars(self)[file + 'Time'], vars(self)[file + str(turbineNumber) + point])
                    ax1.set_xlabel('Time [s]')
                    ax1.set_ylabel(file)
                    ax2.plot(freq, abs(specPos))
                    ax2.set_xlabel("Frequency [Hz]")
                    ax2.set_ylabel("Power")
                    plt.savefig(plotDir + file + '_Frequency.eps')
                    plt.close()
        else:
            fSampl = 1.0 / vars(self)[var + 'dt']
            totTime = vars(self)[var + 'Time'][-1] - vars(self)[var + 'Time'][0]
            dimSampl = np.shape(vars(self)[var + point])[0]
            window = scsig.hann(dimSampl)
            valWindow = vars(self)[var + point] * window
            df = 1.0 / (dimSampl * 1 / fSampl)
            freq = np.arange(0, dimSampl / 2 * df, df)
            specPos = freq * 2.0
            specTot = sc.fft(valWindow)
            specPos[0] = specTot[0] / dimSampl
            if dimSampl % 2 == 0:
                specPos[1:round(dimSampl / 2) - 1] = specTot[1:round(dimSampl / 2) - 1] / [dimSampl / 2]
                specPos[round(dimSampl / 2) - 1] = specTot[round(dimSampl / 2)] / dimSampl
            else:
                specPos[1:round(dimSampl / 2)] = specTot[1:round(dimSampl / 2)] / [dimSampl / 2]
            plt.figure()
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.suptitle(var)
            ax1.plot(vars(self)[var + 'Time'], vars(self)[var + point])
            ax1.set_xlabel('Time [s]')
            ax1.set_ylabel(var)
            ax2.plot(freq, abs(specPos))
            ax2.set_xlabel("Frequency [Hz]")
            ax2.set_ylabel("Power")
            plt.savefig(plotDir + var + '_Frequency.eps')
            plt.close()

