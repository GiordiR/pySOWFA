import os
import re
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import math

import pySOWFA.Turbine as Turbine
import pySOWFA.utils as utils
from pySOWFA.mathUtils import FFT, xcorr_fft
from pySOWFA.plotUtils import plot, plotUtils, plotCompare, plotLog, plotLogCompare, plotLogUtils, plotLogLog, \
                      plotLogLogCompare, plot3D, endPlot, getAxes, getTitle


class OpenFOAM(Turbine):
    """
    Class to postprocess OpenFOAM related output files (e.g. probes, sample)

    :param str turbineName: turbine name [turbine name, windTunnel, precursor, noTurbine]
    :param str probeName: name of the probe set to be post-processed or created
    :param str turbineDir: turbine directory path
    :param str turbineFileName: turbine file name
    """

    def __init__(self, turbineName, probeName, turbineDir=None, turbineFileName=None):
        self.probeName = probeName
        Turbine.__init__(self, turbineName, turbineDir, turbineFileName)

    def makeProbes(self, probeDir='./', probeSets=1, outCtrl='timeStep', outInt=1, tStart=None, tEnd=None, fields=None,
                   x=np.array([(1, 1)]), y=np.array([(1, 1)]), z=np.array([(1, 1)]), nProbes=None, stepProbes=None):
        """
        Creates probes file to sample fields from OpenFOAM simulation

        :param str probeDir: probe's output file directory
        :param int probeSets: number of probe sets
        :param str outCtrl: probe's output control
        :param int outInt: probe's data output interval
        :param int tStart: probe's sampling start time
        :param int tEnd: probe's sampling end time
        :param list(str) fields: probe's output fields
        :param float x: probe's initial x coordinate
        :param float y: probe's initial y coordinate
        :param float z: probe's initial z coordinate
        :param int nProbes: number of probes in a set
        :param float stepProbes: spatial steps between consecutive probe's points
        """

        # Fill default parameters
        if nProbes is None:
            nProbes = np.ones((probeSets, 3))
        elif stepProbes is None:
            stepProbes = np.zeros((probeSets, 3))

        # Change parameter type
        x = x.astype(np.float)
        y = y.astype(np.float)
        z = z.astype(np.float)
        nProbes = nProbes.astype(np.int)
        stepProbes = stepProbes.astype(np.float)

        # Check for directory or make it
        if not os.path.isdir(probeDir):
            print('Making directory: ' + probeDir)
            os.mkdir(probeDir)

        # Check for field ( fields=['u', 'p', ....] )
        if fields is None:
            raise NameError('FieldListError')

        # Open the file and write
        with open(probeDir + str(self.probeName), 'w') as file:
            file.write("{}".format(self.probeName))
            file.write("\n{")
            file.write("\n{:4}{:<30}{}".format('', 'type', 'probes;'))
            file.write("\n{:4}{:<30}{}".format('', 'functionObjectLibs', '("libsampling.so");'))
            file.write("\n{:4}{:<30}{}".format('', 'enabled', 'true;'))
            file.write("\n{:4}{:<30}{}{}".format('', 'probeName', self.probeName, ';'))
            file.write("\n{:4}{:<30}{}{}".format('', 'outputControl', outCtrl, ';'))
            file.write("\n{:4}{:<30}{}{}".format('', 'outputInterval', int(outInt), ';'))
            if tStart and tEnd is not None:
                file.write("\n{:4}{:<30}{}{}".format('', 'timeStart', int(tStart), ';'))
                file.write("\n{:4}{:<30}{}{}".format('', 'timeEnd', int(tEnd), ';'))
            file.write("\n\n{:4}{}".format('', 'fields'))
            file.write("\n{:4}{}".format('', '('))
            for field in fields:
                file.write("\n{:8}{}".format('', field))
            file.write("\n{:4}{}".format('', ');'))
            file.write("\n\n{:4}{}".format('', 'probeLocations'))
            file.write("\n{:4}{}".format('', '('))
            # Write Probe Locations
            minimum = np.zeros((probeSets, 3))
            maximum = np.zeros((probeSets, 3))
            iterValue = np.zeros((probeSets, 3))
            for i in range(0, probeSets):
                count_n = 0
                minimum[i, 0] = x[i, 0]
                minimum[i, 1] = y[i, 0]
                minimum[i, 2] = z[i, 0]
                maximum[i, 0] = x[i, 1]
                maximum[i, 1] = y[i, 1]
                maximum[i, 2] = z[i, 1]
                iterValue[i, 0] = x[i, 0]
                iterValue[i, 1] = y[i, 0]
                iterValue[i, 2] = z[i, 0]
                if not (nProbes[i, :] == np.ones((1, 3))).all():
                    step = np.zeros(3)
                    step[0] = (maximum[i, 0] - minimum[i, 0]) / nProbes[i, 0]
                    step[1] = (maximum[i, 1] - minimum[i, 1]) / nProbes[i, 1]
                    step[2] = (maximum[i, 2] - minimum[i, 2]) / nProbes[i, 2]
                    while (iterValue[i, :] <= maximum[i, :]).all() and (count_n <= max(nProbes[i, :])):
                        file.write("\n{:8}{} {:f} {:f} {:f} {}".format('', '(', iterValue[i, 0], iterValue[i, 1],
                                                                       iterValue[i, 2], ')'))
                        iterValue[i, :] += step
                        count_n += 1
                    print("{} probes in set {} ".format(count_n, i + 1))
                if not (stepProbes[i, :] == np.zeros((1, 3))).all():
                    while (iterValue[i, :] <= maximum[i, :]).all():
                        file.write("\n{:8}{} {:f} {:f} {:f} {}".format('', '(', iterValue[i, 0], iterValue[i, 1],
                                                                       iterValue[i, 2], ')'))
                        iterValue[i, :] += stepProbes[i, :]
                        count_n += 1
                    print("{} probes in set {} ".format(count_n, i + 1))
            file.write("\n{:4}{}".format('', ');'))
            file.write("\n}")

    def readProbes(self, postProcDir=None):
        """
        Read probes data output from OpenFOAM output file

        :param str postProcDir: post-processing OpenFOAM folder path
        """

        # Set post-processing folder path or default one
        if postProcDir is None:
            postProcDir = './postProcessing/'
        probePath = os.path.join(postProcDir, self.probeName, '')

        # Read probes files in post-processing time folders
        timeDir = os.listdir(probePath)
        for t in range(0, len(timeDir)):
            probeDir = os.path.join(probePath, timeDir[t], '')
            files = os.listdir(probeDir)

            # Read probes locations
            with open(probeDir + files[0], 'r') as loc:
                numericPattern = '([+-]?\d+\.?\d* [+-]?\d+\.?\d* [+-]?\d+\.?\d*)'
                probeLoc = []
                n_line = 1
                for line in loc:
                    probeLoc.append(re.findall(numericPattern, line))
                    if 'Probe' not in line:
                        break
            probeLoc = probeLoc[:-2]
            # split sublists
            for i in range(0, len(probeLoc)):
                dummy = [el.split() for el in probeLoc[i]]
                probeLoc[i] = dummy
            # convert to numpy array
            self.probeLoc = np.array(probeLoc, dtype=float).reshape(len(probeLoc), 3)

            # Number of probes
            self.nProbes = len(self.probeLoc)

            # Find number of probe sets and their dimensions
            nProbeSets = 0
            setsDim = []
            old = 0
            for i in range(2, self.nProbes):
                helpSets = 0
                for j in range(0, 3):
                    if (self.probeLoc[i, j] != self.probeLoc[i - 1, j]) or (
                            self.probeLoc[i, j] != self.probeLoc[i - 2, j]):
                        helpSets += 1
                if helpSets == 2 and old != i - 1:
                    nProbeSets += 1
                    setsDim.append(int(i - old))
                    old = i
                if i == self.nProbes - 1:
                    setsDim.append(int(i - old + 1))
            self.nProbeSets = nProbeSets + 1
            self.setsDim = np.array(setsDim)

            # Read field values from file
            scalar_fields, vector_fields, tensor_fields = utils.fieldsOF()
            for file in files:
                if file in scalar_fields:
                    scalar_db = pd.read_csv(probeDir + file, sep='\s+', skiprows=self.nProbes + 2, header=None)
                    values = scalar_db.to_numpy(dtype=float)
                    self.time = values[:, 0]
                    if self.time.shape[0] > 1:
                        self.dt = self.time[1] - self.time[0]
                        self.frq = 1.0 / self.dt
                    vars(self)[file] = values[:, 1:]
                elif file in vector_fields:
                    vector_pattern = '\s+|\(|\)'
                    vector_db = pd.read_csv(probeDir + file, sep=vector_pattern, skiprows=self.nProbes + 2,
                                            header=None, engine='python', keep_default_na=True)
                    vector_db.dropna(how='all', axis=1, inplace=True)
                    values = vector_db.to_numpy(dtype=float)
                    self.time = values[:, 0]
                    if self.time.shape[0] > 1:
                        self.dt = self.time[1] - self.time[0]
                        self.frq = 1.0 / self.dt
                    vars(self)[file + 'x'] = values[:, 1::3]
                    vars(self)[file + 'y'] = values[:, 2::3]
                    vars(self)[file + 'z'] = values[:, 3::3]
                elif file in tensor_fields:
                    tensor_pattern = '\s+|\(|\)'
                    tensor_db = pd.read_csv(probeDir + file, sep=tensor_pattern, skiprows=self.nProbes + 2,
                                            header=None, engine='python', keep_default_na=True)
                    tensor_db.dropna(how='all', axis=1, inplace=True)
                    values = tensor_db.to_numpy(dtype=float)
                    self.time = values[:, 0]
                    if self.time.shape[0] > 1:
                        self.dt = self.time[1] - self.time[0]
                        self.frq = 1.0 / self.dt
                    vars(self)[file + 'xx'] = values[:, 1::6]
                    vars(self)[file + 'xy'] = values[:, 2::6]
                    vars(self)[file + 'xz'] = values[:, 3::6]
                    vars(self)[file + 'yy'] = values[:, 4::6]
                    vars(self)[file + 'yz'] = values[:, 5::6]
                    vars(self)[file + 'zz'] = values[:, 6::6]
                else:
                    raise NameError('FieldTypeError')

        # Calculate fluctuating velocities
        #self.ux = self.Ux - self.UMeanx
        #self.uy = self.Uy - self.UMeany
        #self.uz = self.Uz - self.UMeanz

    def makeSampleDict(self, fields=['U', 'UMean', 'UPrime2Mean'], sampleDictDir=None, setFormat='raw', surfaceFormat='raw', interpolationScheme='cell'):
        """
        Create sampleDict file for OpenFOAM simulation postprocess
        TODO:write sets
        :param list(str) fields: List of fields to be sampled
        :param str sampleDictDir: sampleDict file directory path
        :param str setFormat: output set format [raw, csv, gnuplot, ...]
        :param str interpolationScheme: interpolation scheme [cell, cellPoint, cellPointFace, ...]
        """
        # Check for directory or make it
        if not os.path.isdir(sampleDictDir):
            print('Making directory: ' + sampleDictDir)
            os.mkdir(sampleDictDir)

        # Open the file and write
        with open(sampleDictDir + '/sampleDict', 'w') as file:
            file.write("\n{:<15}{:<30}{}".format('setFormat', setFormat, ';'))
            file.write("\n{:<15}{:<30}{}".format('surfaceFormat', surfaceFormat, ';'))
            file.write("\n{:<15}{:<30}{}".format('interpolationScheme', interpolationScheme, ';'))
            file.write("\n\n}{}".format('fields'))
            file.write("\n{}".format('('))
            for field in fields:
                file.write("\n{:4}{}".format('', field))
            file.write("\n{}".format(');'))
            file.write("\n\n{}".format( 'sets'))
            file.write("\n{}".format('('))
            # Write sets
            '''
            minimum = np.zeros((probeSets, 3))
            maximum = np.zeros((probeSets, 3))
            iterValue = np.zeros((probeSets, 3))
            for i in range(0, probeSets):
                count_n = 0
                minimum[i, 0] = x[i, 0]
                minimum[i, 1] = y[i, 0]
                minimum[i, 2] = z[i, 0]
                maximum[i, 0] = x[i, 1]
                maximum[i, 1] = y[i, 1]
                maximum[i, 2] = z[i, 1]
                iterValue[i, 0] = x[i, 0]
                iterValue[i, 1] = y[i, 0]
                iterValue[i, 2] = z[i, 0]
                if not (nProbes[i, :] == np.ones((1, 3))).all():
                    step = np.zeros(3)
                    step[0] = (maximum[i, 0] - minimum[i, 0]) / nProbes[i, 0]
                    step[1] = (maximum[i, 1] - minimum[i, 1]) / nProbes[i, 1]
                    step[2] = (maximum[i, 2] - minimum[i, 2]) / nProbes[i, 2]
                    while (iterValue[i, :] <= maximum[i, :]).all() and (count_n <= max(nProbes[i, :])):
                        file.write("\n{:8}{} {:f} {:f} {:f} {}".format('', '(', iterValue[i, 0], iterValue[i, 1],
                                                                       iterValue[i, 2], ')'))
                        iterValue[i, :] += step
                        count_n += 1
                    print("{} probes in set {} ".format(count_n, i + 1))
                if not (stepProbes[i, :] == np.zeros((1, 3))).all():
                    while (iterValue[i, :] <= maximum[i, :]).all():
                        file.write("\n{:8}{} {:f} {:f} {:f} {}".format('', '(', iterValue[i, 0], iterValue[i, 1],
                                                                       iterValue[i, 2], ')'))
                        iterValue[i, :] += stepProbes[i, :]
                        count_n += 1
                    print("{} probes in set {} ".format(count_n, i + 1))
            file.write("\n{:4}{}".format('', ');'))
            file.write("\n}")
            '''

    def readSets(self, postProcDir=None, var='Ux'):
        """
        Read sample set output data

        :param str postProcDir: path to the postProcessing directory
        :param str var: variable to be read from file
        """

        # Set post-processing folder path or default one
        if postProcDir is None:
            postProcDir = './postProcessing/'
        setsPath = os.path.join(postProcDir, 'sets')

        # Read probes files in post-processing time folders
        timeDir = os.listdir(setsPath)
        for t in range(0, len(timeDir)):
            setsDir = os.path.join(setsPath, timeDir[t])

        vector_components, tensor_components = utils.fieldsComponentsOF()
        if var in vector_components and not var.startswith('UPrime2Mean'):
            with open(os.path.join(setsDir, self.probeName + '_U_UMean.xy'), "r") as file:
                db = pd.read_csv(file, sep='\s+', skiprows=0, header=None)
                values = db.to_numpy(dtype=float)
                self.xSets = values[:, 0]
                self.ySets = values[:, 1]
                self.zSets = values[:, 2]
                self.Ux = values[:, 3]
                self.Uy = values[:, 4]
                self.Uz = values[:, 5]
                self.UMeanx = values[:, 6]
                self.UMeany = values[:, 7]
                self.UMeanz = values[:, 8]
        elif var.startswith('UPrime2Mean'):
            with open(os.path.join(setsDir, self.probeName + '_UPrime2Mean.xy'), "r") as file:
                db2 = pd.read_csv(file, sep='\s+', skiprows=0, header=None)
                values2 = db2.to_numpy(dtype=float)
                self.xSets = values2[:, 0]
                self.ySets = values2[:, 1]
                self.zSets = values2[:, 2]
                self.UPrime2Meanxx = values2[:, 3]
                self.UPrime2Meanxy = values2[:, 4]
                self.UPrime2Meanxz = values2[:, 5]
                self.UPrime2Meanyy = values2[:, 6]
                self.UPrime2Meanyz = values2[:, 7]
                self.UPrime2Meanzz = values2[:, 8]
        self.setsDim = len(self.xSets)
        self.nProbeSets = 1

    def readWakeExperiment(self, expDir=None, probeSet='cross'):
        """
        Read velocity data from experimental data

        :param str expDir: experimental data file directory path
        :param str probeSet: probe's name
        """

        if expDir is None:
            expDir = '/home/giordi/Scrivania/UNAFLOW/WAKE'

        # WAKE CROSS
        if probeSet == 'cross':
            self.xCross = np.array([5.48, 5.48])
            self.yCross = np.arange(-1.47, 1.740, 0.1)
            self.zCross = np.array([2.085, 2.085])

            wake_cross = pd.read_csv(os.path.join(expDir, 'CW_Steady_V4.dat'), delimiter='\t', header=None).to_numpy()

            self.time = wake_cross[:, 0]
            self.Ucrossx = wake_cross[:-1, 2::3]
            self.Ucrossy = wake_cross[:-1, 3::3]
            self.Ucrossz = wake_cross[:-1, 4::3]

            self.UcrossMeanx = np.sum(self.Ucrossx, axis=0) / len(self.Ucrossx)
            self.UcrossMeany = np.sum(self.Ucrossy, axis=0) / len(self.Ucrossy)
            self.UcrossMeanz = np.sum(self.Ucrossz, axis=0) / len(self.Ucrossz)

        elif probeSet == 'along':
            self.xAlong = np.arange(2.18, 5.81, 0.33)
            self.yAlong = np.array([0.7, 0.9])
            self.zAlong = np.array([2.1, 2.1])

            wake_along = pd.read_csv(os.path.join(expDir, 'AW_Steady_V4.dat'), delimiter='\t', header=None).to_numpy()

            self.time = wake_along[:, 0]
            self.Ucrossx = wake_along[:-1, 2::3]
            self.Ucrossy = wake_along[:-1, 3::3]
            self.Ucrossz = wake_along[:-1, 4::3]

            self.UcrossMeanx = np.sum(self.Ucrossx, axis=0) / len(self.Ucrossx)
            self.UcrossMeany = np.sum(self.Ucrossy, axis=0) / len(self.Ucrossy)
            self.UcrossMeanz = np.sum(self.Ucrossz, axis=0) / len(self.Ucrossz)

    def getTurbulenceIntensity(self, rms='UPrime2Mean'):
        """
        Calculate turbulence intensity index from RMean data or from UPrime2Mean data

        :param str rms: type of rms data [RMean, UPrime2Mean]
        """
        if rms == 'RMean':
            self.TIx = np.sqrt(self.RMeanxx) / self.UMeanx
            self.TIy = np.sqrt(self.RMeanyy) / self.UMeanx
            self.TIz = np.sqrt(self.RMeanzz) / self.UMeanx
            self.TI = np.sqrt(self.RMeanxx + self.RMeanyy + self.RMeanzz) / np.sqrt(
                self.UMeanx ** 2 + self.UMeany ** 2 + self.UMeanz ** 2)
        elif rms == 'UPrime2Mean':
            self.TIx = np.sqrt(self.UPrime2Meanxx) / self.UMeanx
            self.TIy = np.sqrt(self.UPrime2Meanyy) / self.UMeanx
            self.TIz = np.sqrt(self.UPrime2Meanzz) / self.UMeanx
            self.TI = np.sqrt(self.UPrime2Meanxx + self.UPrime2Meanyy + self.UPrime2Meanzz) / np.sqrt(
                self.UMeanx ** 2 + self.UMeany ** 2 + self.UMeanz ** 2)

    def getAddedTurbulenceIntensity(self, TIref=None):
        """
        Calculate added turbulence intensity

        :param TIref: reference turbulence intensity
        """
        self.addTIx = np.sqrt(vars(TIref)['TIx'] ** 2 - self.TIx ** 2)
        self.addTIy = np.sqrt(vars(TIref)['TIy'] ** 2 - self.TIy ** 2)
        self.addTIz = np.sqrt(vars(TIref)['TIz'] ** 2 - self.TIz ** 2)
        self.addTI = np.sqrt(vars(TIref)['TI'] ** 2 - self.TI ** 2)

    def getLESRatio(self):
        pass

    def getPowerProfile(self, z, alpha, normVar):
        """
        Calculate reference velocity power profile

        :param float z: z coordinate
        :param float alpha: power law exponent
        :param list normVar: normalizing variable
        :return Upow: power profile velocity field
        """
        hub_H = 2.09
        # Upow = getattr(normVar[0], normVar[1])[-1] * (z/hub_H)**alpha
        Upow = 4.0 * (z / hub_H) ** alpha
        return Upow

    def plotEnergySpectrum(self, plotDir=None):
        """
        Compute and plot energy spectrum

        :param str plotDir: plot saving directory path
        """

        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        ux = self.ux[:, 0]
        uy = self.uy[:, 0]
        uz = self.uz[:, 0]

        # velocity in frequency domain
        uxfrq, uxamp = FFT(ux, self.frq)
        uyfrq, uyamp = FFT(uy, self.frq)
        uzfrq, uzamp = FFT(uz, self.frq)

        # Generate autocorrelation
        R11, tauR11 = xcorr_fft(ux, norm='biased')
        R22, tauR22 = xcorr_fft(uy, norm='biased')
        R33, tauR33 = xcorr_fft(uz, norm='biased')

        # Time energy spectrum
        S11frq, S11 = FFT(R11, self.frq)
        S22frq, S22 = FFT(R22, self.frq)
        S33frq, S33 = FFT(R33, self.frq)

        plotLogLog(1, S11frq, S11, xlabel=r'$k_x$ [Hz]', ylabel='E(k)', label=r"$S_{11}$", plotDir=plotDir, figName='EnergySpectrumS', title='Energy Spectrum')
        plotLogLogCompare(1, S22frq, S22, label=r"$S_{22}$", plotDir=plotDir, figName='EnergySpectrumS')
        plotLogLogCompare(1, S33frq, S33, label=r"$S_{33}$", plotDir=plotDir, figName='EnergySpectrumS')

        plotLogLog(2, uxfrq, uxamp, xlabel=r'$k_x$ [Hz]', ylabel='E(k)', label=r"$u_{x}$", plotDir=plotDir,
                   figName='EnergySpectrumU', title='Energy Spectrum')
        plotLogLogCompare(2, uxfrq, uyamp, label=r"$u_{y}$", plotDir=plotDir, figName='EnergySpectrumU')
        plotLogLogCompare(2, uxfrq, uzamp, label=r"$u_{z}$", plotDir=plotDir, figName='EnergySpectrumU')
        plotLogLogCompare(2, uxfrq, 0.001 * uxfrq ** (-5 / 3), style='--', label='5/3 law', plotDir=plotDir, figName='EnergySpectrumU')
        endPlot()

    def plotSpaceCorrelations(self, plotDir=None):
        """
        Compute and plot space correlations of the turbulent quantities

        :param str plotDir: plot saving directory path
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        ux = self.ux[-1, :]
        uy = self.uy[-1, :]
        uz = self.uz[-1, :]

        self.r11, self.taur11 = xcorr_fft(ux)
        self.r22, self.taur22 = xcorr_fft(uy)
        self.r33, self.taur33 = xcorr_fft(uz)
        self.r12, self.taur12 = xcorr_fft(ux, uy)
        self.r13, self.taur13 = xcorr_fft(ux, uz)
        self.r23, self.taur23 = xcorr_fft(uy, uz)

        probeStart = 0
        for i in range(0, self.nProbeSets):
            probeEnd = probeStart + self.setsDim[i]
            x = self.probeLoc[probeStart:probeEnd, 0]
            y = self.probeLoc[probeStart:probeEnd, 1]
            z = self.probeLoc[probeStart:probeEnd, 2]
            # Find x axis for plot
            if x[1] != x[2]:
                dx = x[1] - x[0]
                X = x[-1] - x[0]
                xVar11 = self.taur11 * dx / X
                xVar22 = self.taur22 * dx / X
                xVar33 = self.taur33 * dx / X
                xlabel = 'x/Lx'
            elif y[1] != y[2]:
                dy = y[1] - y[0]
                Y = y[-1] - y[0]
                xVar11 = self.taur11 * dy / Y
                xVar22 = self.taur22 * dy / Y
                xVar33 = self.taur33 * dy / Y
                xlabel = 'y/Ly'
            elif z[1] != z[2]:
                dz = z[1] - z[0]
                Z = z[-1] - z[0]
                xVar11 = self.taur11 * dz / Z
                xVar22 = self.taur22 * dz / Z
                xVar33 = self.taur33 * dz / Z
                xlabel = 'z/Lz'

        plot(1, xVar11, self.r11, xLabel=xlabel, yLabel=r'$\rho$', label=r'$\rho_{11}$', plotDir=plotDir,
             figName='space_correlations', title='Space correlations')
        plotCompare(1, xVar22, self.r22, label=r'$\rho_{22}$', plotDir=plotDir, figName='space_correlations')
        plotCompare(1, xVar33, self.r33, label=r'$\rho_{33}$', plotDir=plotDir, figName='space_correlations')
        endPlot()

    def plotTimeCorrelations(self, plotDir=None):
        """
        Compute and plot time correlations of the turbulent quantities

        :param str plotDir: plot saving directory path
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        ux = self.ux[:, 0]
        uy = self.uy[:, 0]
        uz = self.uz[:, 0]

        self.r11, self.taur11 = xcorr_fft(ux)
        self.r22, self.taur22 = xcorr_fft(uy)
        self.r33, self.taur33 = xcorr_fft(uz)
        self.r12, self.taur12 = xcorr_fft(ux, uy)
        self.r13, self.taur13 = xcorr_fft(ux, uz)
        self.r23, self.taur23 = xcorr_fft(uy, uz)

        plot(1, self.taur11, self.r11, xLabel=r'$\tau$', yLabel=r'$\rho$', label=r'$\rho_{11}$', plotDir=plotDir,
             figName='time_correlations', title='Time correlations')
        plotCompare(1, self.taur22, self.r22, label=r'$\rho_{22}$', plotDir=plotDir, figName='time_correlations')
        plotCompare(1, self.taur33, self.r33, label=r'$\rho_{33}$', plotDir=plotDir, figName='time_correlations')
        endPlot()

    def plotWindProfile(self, figID, plotDir=None, var='p', normVar=None, normalize=False, compareID=None):
        """
        Compute and plot wind profile

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: field to be plotted
        :param str normVar: normalizing variable
        :param bool normalize: boolean to use normalization
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = os.path.join('./postProcessing', self.probeName, 'plots')
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        # Plot every set of probes
        probeStart = 0
        for i in range(0, self.nProbeSets):
            probeEnd = probeStart + self.setsDim[i]
            x = self.probeLoc[probeStart:probeEnd, 0]
            y = self.probeLoc[probeStart:probeEnd, 1]
            z = self.probeLoc[probeStart:probeEnd, 2]
            # Find x-axis for plot
            if x[1] != x[2]:
                yVar = x
                ylabel = 'x [m]'
                ylim = None
            elif y[1] != y[2]:
                yVar = y
                ylabel = 'y [m]'
                ylim = None
            elif z[1] != z[2]:
                yVar = z
                ylabel = 'z [m]'
                ylim = [0, 3]
            if normalize is False:
                xVar = vars(self)[var][-1, probeStart:probeEnd]
                xlabel = getAxes(var)
                figName = var + str(figID)
            else:
                xVar = np.divide(vars(self)[var][-1, probeStart:probeEnd], getattr(normVar[0], normVar[1])[-1])
                xlabel = getAxes(var) + '/' + normVar[1]
                figName = var + str(figID) + '_norm'
            if compareID is not None:
                plotCompare(compareID, xVar, yVar, var, plotDir, figName)
            else:
                if var == 'UMeanx':
                    plotUtils(figID, self.getPowerProfile(z, 0.2, normVar), yVar, 'Power law')
                plot(figID, xVar, yVar, xlabel, ylabel, var, plotDir, figName, ylim=ylim)
            probeStart = probeEnd

    def plotWakeProfile(self, figID, plotDir=None, sampleType='probe', var='p', normVar=None, xLim=None, filter=None, interp=None, compareID=None):
        """
        Compute and plot wake profile

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str sampleType: type of sample method used [probe, sets]
        :param str var: field to be plotted
        :param str normVar: normalizing variable
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = os.path.join('./postProcessing', self.probeName, 'plots')
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        # Plot every set of probes
        probeStart = 0
        for i in range(0, self.nProbeSets):
            probeEnd = probeStart + (self.setsDim if sampleType=='sets' else self.setsDim[i])
            if sampleType == 'probe':
                x = self.probeLoc[probeStart:probeEnd, 0]
                y = self.probeLoc[probeStart:probeEnd, 1]
                z = self.probeLoc[probeStart:probeEnd, 2]
            elif sampleType == 'sets':
                x = self.xSets
                y = self.ySets
                z = self.zSets
                vars(self)[var] = np.tile(vars(self)[var], (2, 1))

            if normVar is None:
                xVar = vars(self)[var][-1, probeStart:probeEnd]
                xlabel = getAxes(var)
                figName = '/' + var + str(figID)
            else:
                xVar = np.divide(vars(self)[var][-1, probeStart:probeEnd], getattr(normVar[0], normVar[1])[-1])
                xlabel = getAxes(var) + '/' + getAxes(normVar[1])
                figName = var + str(figID) + '_norm'

            xVar, x, y, z = utils.cleanDoubleData(xVar, x, y, z)

            # Specify x-axis limits manually
            if xLim is not None:
                xlim = [xLim[0], xLim[1]]

            # Find x-axis for plot
            if x[1] != x[2]:
                yVar = x / self.rotor_D
                ylabel = 'x/D'
                xlim = None
                ylim = None
            elif y[1] != y[2]:
                yVar = y / self.rotor_D
                ylabel = 'y/D'
                ylim = [-1, 1]
                if xLim is None:
                    xlim = [min(xVar) - 0.5, max(xVar) + 0.5]
                plotUtils(figID, [xlim[0], xlim[1]], [0.5, 0.5])
                plotUtils(figID, [xlim[0], xlim[1]], [0, 0])
                plotUtils(figID, [xlim[0], xlim[1]], [-0.5, -0.5])
                figName = figName + '_H'
            elif z[1] != z[2]:
                yVar = z
                ylabel = 'z [m]'
                ylim = [0, 3.5]
                if xLim is None:
                    xlim = [min(xVar) - 0.5, max(xVar) + 0.5]
                plotUtils(figID, [xlim[0], xlim[1]],
                          [self.nacelle_H - self.rotor_D / 2, self.nacelle_H - self.rotor_D / 2])
                plotUtils(figID, [xlim[0], xlim[1]], [self.nacelle_H, self.nacelle_H])
                plotUtils(figID, [xlim[0], xlim[1]],
                          [self.nacelle_H + self.rotor_D / 2, self.nacelle_H + self.rotor_D / 2])
                figName = figName + '_V'
            else:
                raise ValueError("Error! Probe points are equal!")

            # Smooth data
            if interp is not None:
                yInterp = np.linspace(min(yVar), max(yVar), num=200, endpoint=True)
                f = interp1d(yVar, xVar, kind=interp)
                xVar = f(yInterp)
                yVar = yInterp
            if filter is not None:
                xVar = savgol_filter(xVar, filter[0], filter[1])

            # distance = round((x[0] / self.rotor_D) / 0.5) * 0.5
            distance = math.floor((x[0] / self.rotor_D) * 10) / 10
            title = 'Distance: ' + str(distance) + 'D'
            label = self.turbineName
            if compareID is not None:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)
            else:
                if var.startswith('TI'):
                    xlim = [0.0, 0.1]
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            probeStart = probeEnd

        self.xVarWake = xVar
        self.yVarWake = yVar
        self.ylabelWake = ylabel

    def plotWakeExperiment(self, compareID, plotDir=None, expProbe='probe_exp_cross', norm=False, interp=None, filter=None):
        """
        Plot wake profile from experimental data

        :param int compareID: figure identification number for comparison
        :param str plotDir: plot saving directory path
        :param str expProbe: experimental probe name
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        if expProbe == 'probe_exp_cross':
            if norm:
                xVar = self.UcrossMeanx / 4.0
            else:
                xVar = self.UcrossMeanx
            yVar = self.yCross / self.rotor_D

            # Smooth data
            if interp is not None:
                yInterp = np.linspace(min(yVar), max(yVar), num=200, endpoint=True)
                f = interp1d(yVar, xVar, kind=interp)
                xVar = f(yInterp)
                yVar = yInterp
            if filter is not None:
                xVar = savgol_filter(xVar, filter[0], filter[1])


            label = 'Experimental Data'
            figName = '/expComparison_cross_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)
        elif expProbe == 'probe_exp_along':
            xVar = self.UcrossMeanx
            yVar = self.xAlong / self.rotor_D
            # Smooth data
            if filter is not None:
                xVar = savgol_filter(xVar, filter[0], filter[1])
            label = 'Experimental Data'
            figName = '/expComparison_along_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotWakeDeficitProfile(self, figID, plotDir=None, var='p', normVar=None, compareID=None):
        """
        Compute and plot wake deficit profile

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str var: field to be plotted
        :param str normVar: normalizing variable
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        # Plot every set of probes
        probeStart = 0
        for i in range(0, self.nProbeSets):
            probeEnd = probeStart + self.setsDim[i]
            x = self.probeLoc[probeStart:probeEnd, 0]
            y = self.probeLoc[probeStart:probeEnd, 1]
            z = self.probeLoc[probeStart:probeEnd, 2]
            # Find x axis for plot
            if x[1] != x[2]:
                yVar = x / self.rotor_D
                ylabel = 'x/D'
                xlim = None
                ylim = None
            elif y[1] != y[2]:
                yVar = y / self.rotor_D
                ylabel = 'y/D'
                ylim = [-1, 1]
                xlim = [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                        max(vars(self)[var][-1, probeStart:probeEnd]) + 1]
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1], [0.5, 0.5])
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1], [0, 0])
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1], [-0.5, -0.5])

            elif z[1] != z[2]:
                yVar = z
                ylabel = 'z [m]'
                ylim = [0, 3.5]
                xlim = [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                        max(vars(self)[var][-1, probeStart:probeEnd]) + 1]
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1],
                          [self.nacelle_H - self.rotor_D / 2, self.nacelle_H - self.rotor_D / 2])
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1], [self.nacelle_H, self.nacelle_H])
                plotUtils(figID, [min(vars(self)[var][-1, probeStart:probeEnd]) - 1,
                                  max(vars(self)[var][-1, probeStart:probeEnd]) + 1],
                          [self.nacelle_H + self.rotor_D / 2, self.nacelle_H + self.rotor_D / 2])
            xVar = np.divide(getattr(normVar[0], normVar[1])[-1] - vars(self)[var][-1, probeStart:probeEnd],
                             getattr(normVar[0], normVar[1])[-1])
            xlabel = r'$\Delta$' + getAxes(var) + '/' + normVar[1]
            title = 'Distance: ' + self.probeName[-2:]
            figName = var + str(figID) + '_deficit'
            label = self.turbineName
            if compareID is not None:
                if normVar is None:
                    xVar = vars(self)[var][-1, probeStart:probeEnd]
                    plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                else:
                    xVar = np.divide(vars(self)[var][-1, probeStart:probeEnd], getattr(normVar[0], normVar[1])[-1])
                    plotCompare(compareID, xVar, yVar, var, plotDir, figName)
            else:
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            probeStart = probeEnd

    def plotWakeError(self, figID, plotDir=None, probeRef=None, probeCompare=None, labelCompare=None, compareID=None):
        """
        Compute and plot wake error estimate between the calculated values and the reference ones

        :param int figID: figure identification number
        :param str plotDir: plot saving directory path
        :param str probeRef: reference object direction ['x variable', 'y variable']
        :param str probeCompare: comparative object variable
        :param int compareID: figure identification number for comparison
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        Uref = getattr(probeCompare[2], probeCompare[3])
        if len(Uref) == 1:
            Uref = getattr(probeCompare[2], probeCompare[3])
        else:
            Uref = getattr(probeCompare[2], probeCompare[3])[-1]

        Ucomp = np.divide(getattr(probeCompare[0], probeCompare[1])[-1, 1:],
                          Uref)
        Umes = np.divide(getattr(probeRef[0], probeRef[1])[1:], 4.0)

        xVar = (Umes - Ucomp) / Umes
        yVar = getattr(probeRef[0], probeRef[2]) / self.rotor_D
        yVar = yVar[1:]
        xlabel = r'$\epsilon$ '
        ylabel = getattr(probeCompare[0], 'ylabelWake')
        #label = labelCompare + ' - ' + r'$\epsilon_{mean} = $' + str(xMean.round(decimals=3))
        title = 'Relative error distribution'
        figName = 'wakeError_' + str(figID)
        if compareID is None:
            plot(figID, xVar, yVar, xlabel, ylabel, labelCompare, plotDir, figName, title=title)
        else:
            plotCompare(compareID, xVar, yVar, labelCompare, plotDir, figName)

        xMean = np.mean(xVar)
        N = len(getattr(probeRef[0], probeRef[1]))
        eRMS = np.sqrt(np.sum((Ucomp - Umes) ** 2)) / N

        with open(os.path.join(plotDir, 'ErrorStatistics.txt'),'a+') as file:
            file.write("\n{0:*^50s}".format('CASE ' + labelCompare))
            file.write("\n")
            file.write("\n{0:<50s}".format('Averaged Relative Velocity Error'))
            file.write("{0:>60f}".format(xMean))
            file.write("\n{0:<50s}".format('Root Mean Squared Error'))
            file.write("{0:>60f}".format(eRMS))



    def plotResiduals(self, plotDir=None, var='UMeanx'):
        """
        Compute and plot residuals of the simulations

        :param str plotDir: plot saving directory path
        :param str var: field to be plotted
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        xVar = self.probeLoc[:, 2]
        yVar = self.time
        zVar = vars(self)[var]

        # Plot residuals
        xLabel = 'Height [m]'
        yLabel = 'Time [s]'
        zLabel = getAxes(var)

        plot3D(1, xVar, yVar, zVar, xLabel, yLabel, zLabel, plotDir=plotDir, figName='Residual_'+var)
        endPlot()

