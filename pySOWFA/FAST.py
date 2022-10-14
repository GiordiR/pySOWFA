import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as scsig
import scipy as sc
import string
from scipy.io import loadmat
from scipy.signal import savgol_filter

import Turbine
import utils
from mathUtils import FFT, xcorr_fft
from plotUtils import plot, plotUtils, plotCompare, plotLog, plotLogCompare, plotLogUtils, plotLogLog, \
                      plotLogLogCompare, plot3D, endPlot, getAxes, getTitle


class FAST(Turbine):
    """
    Class to postprocess OpenFAST related output files (e.g. Turbine output, ElastoDyn outputs)

    :param str turbineName: turbine name [turbine name, windTunnel, precursor, noTurbine]
    :param str turbineDir: turbine directory path
    :param str turbineFileName: turbine file name
    """

    def __init__(self, turbineName, turbineDir=None, turbineFileName=None):
        Turbine.__init__(self, turbineName, turbineDir, turbineFileName)

    def readBladeProp(self, fastDir=None):
        """
        Read Wind Turbine blade properties from OpenFAST input files

        :param str fastDir: OpenFAST case directory path
        """
        if fastDir is None:
            self.fastDir = './WTM'
        else:
            self.fastDir = fastDir

        with open(self.fastDir + '/Blade_AeroDyn.dat', "r") as AB:
            dbAB = pd.read_csv(AB, skiprows=6, sep='\s+', header=None, error_bad_lines=False)
            dataAB = dbAB.to_numpy(dtype=float)
            bldSpn = dataAB[:, 0]
            self.normBldSpn = bldSpn / bldSpn[-1]

    def readOutput(self, outDir=None, outFileName=None):
        """
        Read OpenFAST output

        :param str outDir: output directory path
        :param str outFileName: output file name
        """
        if outDir is None:
            self.outDir = './FASTout/'
        else:
            self.outDir = outDir

        if outFileName is None:
            outFileName = 'WTM.out'

        with open(self.outDir + outFileName, "r") as file:
            db_FAST = pd.read_csv(file, sep='\t', skiprows=6, header=None, low_memory=False)
            db_FAST = db_FAST.rename(columns=db_FAST.iloc[0])
            db_FAST = db_FAST.drop(db_FAST.index[0:2])
            self.dataNames = db_FAST.columns
            data_FAST = db_FAST.to_numpy(dtype=float)
            for i in range(0, db_FAST.shape[1]):
                vars(self)[db_FAST.columns[i]] = data_FAST[:, i]

        self.dt = self.Time[1] - self.Time[0]
        self.frq = 1 / self.dt

        # Find the number of blade nodes
        self.bladeNodes = 0
        for var in self.dataNames:
            if var.startswith('AB1'):
                if var.endswith('Re'):
                    self.bladeNodes += 1

        # Create spanwise values from node values
        self.bladeRe = []
        self.bladeM = []
        self.bladeAlpha = []
        self.bladeCl = []
        self.bladeCd = []
        self.bladeCm = []
        for var in self.dataNames:
            if var.startswith('AB1'):
                if var.endswith('Re'):
                    self.bladeRe.append(vars(self)[var][-1])
                elif var.endswith('M'):
                    self.bladeM.append(vars(self)[var][-1])
                elif var.endswith('Alpha'):
                    self.bladeAlpha.append(vars(self)[var][-1])
                elif var.endswith('Cl'):
                    self.bladeCl.append(vars(self)[var][-1])
                elif var.endswith('Cd'):
                    self.bladeCd.append(vars(self)[var][-1])
                elif var.endswith('Cm'):
                    self.bladeCm.append(vars(self)[var][-1])
        self.bladeRe = np.array(self.bladeRe, dtype=float)
        self.bladeM = np.array(self.bladeM, dtype=float)
        self.bladeAlpha = np.array(self.bladeAlpha, dtype=float)
        self.bladeCl = np.array(self.bladeCl, dtype=float)
        self.bladeCd = np.array(self.bladeCd, dtype=float)
        self.bladeCm = np.array(self.bladeCm, dtype=float)

    def readForceExperiment(self, expDir=None):
        """
        Read force output from experimental data

        :param str expDir: experimental data directory path
        """
        if expDir is None:
            expDir = '/home/giordi/Desktop/File_galleria/POLIMI_UNAFLOW_DATA'

        force = loadmat(expDir + '/' + 'FORCE/Surge_V4A0F0.mat')
        self.time = force['t']
        self.TBFx = force['FTBxyz'][:, 0]
        self.TBFy = force['FTBxyz'][:, 1]
        self.TBFz = force['FTBxyz'][:, 2]
        self.TBMx = force['FTBxyz'][:, 3]
        self.TBMy = force['FTBxyz'][:, 4]
        self.TBMz = force['FTBxyz'][:, 5]

        self.TTFx = force['FTTxyz'][:, 0]
        self.TTFy = force['FTTxyz'][:, 1]
        self.TTFz = force['FTTxyz'][:, 2]
        self.TTMx = force['FTTxyz'][:, 3]
        self.TTMy = force['FTTxyz'][:, 4]
        self.TTMz = force['FTTxyz'][:, 5]

    def plotBladeTipDeflections(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot blade tip deflections

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Tip'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = 'Tip'
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = 'Tip'
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeTipDeflectionsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot blade tip deflections power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Tip'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if var.startswith('TipDxb') or var.startswith('TipDyb'):
                        for i in range(1, 4):
                            rotFrq = i * 241.0 / 60.0
                            xUtils = [rotFrq, rotFrq]
                            yUtils = [min(yVar), max(yVar)]
                            plotLogUtils(figID, xUtils, yUtils)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if var.startswith('TipDxb') or var.startswith('TipDyb'):
                for i in range(1, 4):
                    rotFrq = i * 241.0 / 60.0
                    xUtils = [rotFrq, rotFrq]
                    yUtils = [min(yVar), max(yVar)]
                    plotLogUtils(figID, xUtils, yUtils)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeRoot(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot blade root output variables

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Root'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = 'Root'
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = 'Root'
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeRootPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot blade root output variables power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Root'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeOverTime(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot blade output variables over time

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Bl'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)
            figID += 1

    def plotBladeOverSpan(self, figID, plotDir=None, window=None, compareID=None):
        """
        Plot blade output variables over blade span

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param list(float) window: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        # Blade span
        '''
        blade_R = self.rotor_D/2.0
        bldSpn = np.linspace(0, blade_R, self.bladeNodes)
        normBldSpn = bldSpn/blade_R
        '''

        yVariables = ['bladeRe', 'bladeM', 'bladeAlpha', 'bladeCl', 'bladeCd', 'bladeCm']

        for i in range(0, len(yVariables)):
            xVar = self.normBldSpn
            yVar = vars(self)[yVariables[i]]
            xlabel = 'r/R'
            ylabel = getAxes(yVariables[i])
            label = yVariables[i]
            title = getTitle(yVariables[i])
            figName = yVariables[i] + '_' + str(figID)
            if window is not None:
                xlim = [window[0], window[1]]
            else:
                xlim = None
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)
            figID += 1

    def plotRotor(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot rotor output variables

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Rt'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_Time_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_Time_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotRotorPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot rotor output variables power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('Rt'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerBaseLoads(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot tower base loads

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('TwrBs'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerBaseLoadsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot tower base loads power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('TwrBs'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopDisplacements(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot tower top displacements

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('TT'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopDisplacementsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot tower top displacements power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('TT'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopLoads(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot tower top loads

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('YawBr'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopLoadsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):
        """
        Plot tower top loads power spectral density

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('YawBr'):
                    xVar, yVar = FFT(vars(self)[var], self.frq)
                    title = getTitle(var) + " - Frequency Spectrum"
                    xlabel = 'frequency [Hz]'
                    ylabel = 'PSD'
                    label = var
                    figName = var + '_PSD_' + str(figID)
                    if compareID is None:
                        plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar, yVar = FFT(vars(self)[var], self.frq)
            title = getTitle(var) + " - Frequency Spectrum"
            xlabel = 'frequency [Hz]'
            ylabel = 'PSD'
            label = var
            figName = var + '_PSD_' + str(figID)
            if compareID is None:
                plotLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotLowSpeedShaft(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot low speed shaft variables

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('LSS') or var.startswith('RotPwr'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotHighSpeedShaft(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):
        """
        Plot high speed shaft variables

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: variable to be plotted
        :param float ylim: plot y-axis limits
        :param float xlim: plot x-axis limits
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = self.outDir + 'plot/'

        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if var == 'all':
            for var in self.dataNames:
                if var.startswith('HSS') or var.startswith('GenSpeed'):
                    xVar = self.Time
                    yVar = vars(self)[var]
                    xlabel = 'Time [s]'
                    ylabel = getAxes(var)
                    label = var
                    title = getTitle(var)
                    figName = var + '_' + str(figID)
                    if ylim is True:
                        if min(yVar) < 0:
                            low = 1.2 * min(yVar)
                        else:
                            low = 0.3 * min(yVar)
                        if max(yVar) < 0:
                            up = 0.3 * max(yVar)
                        else:
                            up = 1.2 * max(yVar)
                        ylim = [low, up]
                    if compareID is None:
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                    else:
                        plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                    figID += 1
        else:
            xVar = self.Time
            yVar = vars(self)[var]
            xlabel = 'Time [s]'
            ylabel = getAxes(var)
            label = var
            title = getTitle(var)
            figName = var + '_' + str(figID)
            if ylim is True:
                if min(yVar) < 0:
                    low = 1.2 * min(yVar)
                else:
                    low = 0.3 * min(yVar)
                if max(yVar) < 0:
                    up = 0.3 * max(yVar)
                else:
                    up = 1.2 * max(yVar)
                ylim = [low, up]
            if compareID is None:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def generateStatistics(self, statFile=None):
        """
        Generate rotor statistics file

        :param str statFile: rotor statistics file name
        """
        if statFile is None:
            statFile = './FASTstatistics.txt'

        ## BLADE
        # TIP DISPLACEMENTS
        # Max tip displacements
        maxTipDxc1 = max(self.TipDxc1)
        maxTipDyc1 = max(self.TipDyc1)
        maxTipDzc1 = max(self.TipDzc1)
        maxTipDxb1 = max(self.TipDxb1)
        maxTipDyb1 = max(self.TipDyb1)

        # Min tip displacements
        minTipDxc1 = min(self.TipDxc1)
        minTipDyc1 = min(self.TipDyc1)
        minTipDzc1 = min(self.TipDzc1)
        minTipDxb1 = min(self.TipDxb1)
        minTipDyb1 = min(self.TipDyb1)

        # Mean tip displacements
        meanTipDxc1 = np.mean(self.TipDxc1)
        meanTipDyc1 = np.mean(self.TipDyc1)
        meanTipDzc1 = np.mean(self.TipDzc1)
        meanTipDxb1 = np.mean(self.TipDxb1)
        meanTipDyb1 = np.mean(self.TipDyb1)

        # ROOT FORCES AND MOMENTS
        # Max root forces and moments
        maxRootFxb1 = max(self.RootFxb1)
        maxRootFyb1 = max(self.RootFyb1)
        maxRootMxb1 = max(self.RootMxb1)
        maxRootMyb1 = max(self.RootMyb1)
        maxRootMzb1 = max(self.RootMzb1)

        # Min root forces and moments
        minRootFxb1 = min(self.RootFxb1)
        minRootFyb1 = min(self.RootFyb1)
        minRootMxb1 = min(self.RootMxb1)
        minRootMyb1 = min(self.RootMyb1)
        minRootMzb1 = min(self.RootMzb1)

        # Mean root forces and moments
        meanRootFxb1 = np.mean(self.RootFxb1)
        meanRootFyb1 = np.mean(self.RootFyb1)
        meanRootMxb1 = np.mean(self.RootMxb1)
        meanRootMyb1 = np.mean(self.RootMyb1)
        meanRootMzb1 = np.mean(self.RootMzb1)

        ## TOWER
        # TOWER TOP DISPLACEMENTS
        # Max tower top displacement
        maxTTDspFA = max(self.TTDspFA)
        maxTTDspSS = max(self.TTDspSS)

        # Min tower top displacement
        minTTDspFA = min(self.TTDspFA)
        minTTDspSS = min(self.TTDspSS)

        # Mean tower top displacement
        meanTTDspFA = np.mean(self.TTDspFA)
        meanTTDspSS = np.mean(self.TTDspSS)

        # TOWER TOP FORCES AND MOMENTS
        # Max tower top displacement
        maxYawBrFxp = max(self.YawBrFxp)
        maxYawBrFyp = max(self.YawBrFyp)

        # Min tower top displacement
        minYawBrFxp = min(self.YawBrFxp)
        minYawBrFyp = min(self.YawBrFyp)

        # Mean tower top displacement
        meanYawBrFxp = np.mean(self.YawBrFxp)
        meanYawBrFyp = np.mean(self.YawBrFyp)

        ## ROTOR
        # Mean rotor parameters
        meanRotPwr = np.mean(self.RotPwr, dtype=float)
        meanRtAeroPwr = np.mean(self.RtAeroPwr, dtype=float)
        meanRtAeroFxh = np.mean(self.RtAeroFxh, dtype=float)
        meanRtAeroMxh = np.mean(self.RtAeroMxh, dtype=float)

        with open(statFile, "w+") as file:
            file.write("\n{0:#^50s}".format('FAST STATISTICS'))
            file.write("\n")
            file.write("\n")
            file.write("\n{0:*^50s}".format('BLADE'))
            file.write("\n{0:-^50s}".format('Tip Displacements'))
            file.write("\n")
            file.write("\n{0:<40s}".format('Max Out-Of-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(maxTipDxc1))
            file.write("\n{0:<40s}".format('Max In-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(maxTipDyc1))
            file.write("\n{0:<40s}".format('Max Axial Tip Deflection [m]'))
            file.write("{0:>50f}".format(maxTipDzc1))
            file.write("\n{0:<40s}".format('Max Flapwise Tip Displacement [m]'))
            file.write("{0:>50f}".format(maxTipDxb1))
            file.write("\n{0:<40s}".format('Max Edgewise Tip Displacement [m]'))
            file.write("{0:>50f}".format(maxTipDyb1))
            file.write("\n")
            file.write("\n{0:<40s}".format('Min Out-Of-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(minTipDxc1))
            file.write("\n{0:<40s}".format('Min In-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(minTipDyc1))
            file.write("\n{0:<40s}".format('Min Axial Tip Deflection [m]'))
            file.write("{0:>50f}".format(minTipDzc1))
            file.write("\n{0:<40s}".format('Min Flapwise Tip Displacement [m]'))
            file.write("{0:>50f}".format(minTipDxb1))
            file.write("\n{0:<40s}".format('Min Edgewise Tip Displacement [m]'))
            file.write("{0:>50f}".format(minTipDyb1))
            file.write("\n")
            file.write("\n{0:<40s}".format('Mean Out-Of-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(meanTipDxc1))
            file.write("\n{0:<40s}".format('Mean In-Plane Tip Displacement [m]'))
            file.write("{0:>50f}".format(meanTipDyc1))
            file.write("\n{0:<40s}".format('Mean Axial Tip Deflection [m]'))
            file.write("{0:>50f}".format(meanTipDzc1))
            file.write("\n{0:<40s}".format('Mean Flapwise Tip Displacement [m]'))
            file.write("{0:>50f}".format(meanTipDxb1))
            file.write("\n{0:<40s}".format('Mean Edgewise Tip Displacement [m]'))
            file.write("{0:>50f}".format(meanTipDyb1))
            file.write("\n")
            file.write("\n{0:-^50s}".format('Root Forces and Moments'))
            file.write("\n")
            file.write("\n{0:<40s}".format('Max flapwise blade root shear force [N]'))
            file.write("{0:>50f}".format(maxRootFxb1))
            file.write("\n{0:<40s}".format('Max edgewise blade root shear force [N]'))
            file.write("{0:>50f}".format(maxRootFyb1))
            file.write("\n{0:<40s}".format('Max edgewise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(maxRootMxb1))
            file.write("\n{0:<40s}".format('Max flapwise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(maxRootMyb1))
            file.write("\n{0:<40s}".format('Max pitching moment [kN-m]'))
            file.write("{0:>50f}".format(maxRootMzb1))
            file.write("\n")
            file.write("\n{0:<40s}".format('Min flapwise blade root shear force [N]'))
            file.write("{0:>50f}".format(minRootFxb1))
            file.write("\n{0:<40s}".format('Min edgewise blade root shear force [N]'))
            file.write("{0:>50f}".format(minRootFyb1))
            file.write("\n{0:<40s}".format('Min edgewise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(minRootMxb1))
            file.write("\n{0:<40s}".format('Min flapwise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(minRootMyb1))
            file.write("\n{0:<40s}".format('Min pitching moment [kN-m]'))
            file.write("{0:>50f}".format(minRootMzb1))
            file.write("\n")
            file.write("\n{0:<40s}".format('Mean flapwise blade root shear force [N]'))
            file.write("{0:>50f}".format(meanRootFxb1))
            file.write("\n{0:<40s}".format('Mean edgewise blade root shear force [N]'))
            file.write("{0:>50f}".format(meanRootFyb1))
            file.write("\n{0:<40s}".format('Mean edgewise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(meanRootMxb1))
            file.write("\n{0:<40s}".format('Mean flapwise blade root moment [kN-m]'))
            file.write("{0:>50f}".format(meanRootMyb1))
            file.write("\n{0:<40s}".format('Mean pitching moment [kN-m]'))
            file.write("{0:>50f}".format(meanRootMzb1))
            file.write("\n")
            file.write("\n")
            file.write("\n{0:*^50s}".format('TOWER'))
            file.write("\n{0:-^50s}".format('Tower Top Displacements'))
            file.write("\n")
            file.write("\n{0:<40s}".format('Max Tower Top Fore-Aft Deflections [m]'))
            file.write("{0:>50f}".format(maxTTDspFA))
            file.write("\n{0:<40s}".format('Max Tower Top Side-To-Side Deflections [m]'))
            file.write("{0:>50f}".format(maxTTDspSS))
            file.write("\n")
            file.write("\n{0:<40s}".format('Min Tower Top Fore-Aft Deflections [m]'))
            file.write("{0:>50f}".format(minTTDspFA))
            file.write("\n{0:<40s}".format('Min Tower Top Side-To-Side Deflections [m]'))
            file.write("{0:>50f}".format(minTTDspSS))
            file.write("\n")
            file.write("\n{0:<40s}".format('Mean Tower Top Fore-Aft Deflections [m]'))
            file.write("{0:>50f}".format(meanTTDspFA))
            file.write("\n{0:<40s}".format('Mean Tower Top Side-To-Side Deflections [m]'))
            file.write("{0:>50f}".format(meanTTDspSS))
            file.write("\n")
            file.write("\n{0:-^50s}".format('Tower Top Forces'))
            file.write("\n")
            file.write("\n{0:<40s}".format('Max Tower Top Fore-Aft shear force [N]'))
            file.write("{0:>50f}".format(maxYawBrFxp))
            file.write("\n{0:<40s}".format('Max Tower Top Side-To-Side shear force [N]'))
            file.write("{0:>50f}".format(maxYawBrFyp))
            file.write("\n")
            file.write("\n{0:<40s}".format('Min Tower Top Fore-Aft shear force [N]'))
            file.write("{0:>50f}".format(minYawBrFxp))
            file.write("\n{0:<40s}".format('Min Tower Top Side-To-Side shear force [N]'))
            file.write("{0:>50f}".format(minYawBrFyp))
            file.write("\n")
            file.write("\n{0:<40s}".format('Mean Tower Top Fore-Aft shear force [N]'))
            file.write("{0:>50f}".format(meanYawBrFxp))
            file.write("\n{0:<40s}".format('Mean Tower Top Side-To-Side shear force [N]'))
            file.write("{0:>50f}".format(meanYawBrFyp))
            file.write("\n")
            file.write("\n")
            file.write("\n{0:*^50s}".format('ROTOR'))
            file.write("\n")
            file.write("\n{0:<40s}".format('Rotor Power [kW]'))
            file.write("{0:>50f}".format(meanRotPwr))
            file.write("\n{0:<40s}".format('Rotor Aerodynamic Power [W]'))
            file.write("{0:>50f}".format(meanRtAeroPwr))
            file.write("\n{0:<40s}".format('Rotor Thrust [N]'))
            file.write("{0:>50f}".format(meanRtAeroFxh))
            file.write("\n{0:<40s}".format('Rotor Torque [N-m]'))
            file.write("{0:>50f}".format(meanRtAeroMxh))
