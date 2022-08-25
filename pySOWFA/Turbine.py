import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as scsig
import scipy as sc
import string
from scipy.io import loadmat

import utils
from mathUtils import FFT, xcorr_fft
from plotUtils import plot, plotUtils, plotCompare, plotLogLog, plotLogLogCompare, endPlot, getAxes, getTitle

class Turbine(object):
    """
    Tubine class initialization: read the turbine parameters given the turbine name and directory

    :param str turbineName: turbine name [turbine name, windTunnel, precursor, noTurbine]
    :param str turbineDir: turbine directory path
    :param str turbineFileName: turbine file name
    """

    def __init__(self, turbineName, turbineDir, turbineFileName):
        self.turbineName = turbineName

        if turbineName == 'windTunnel' or turbineName == 'precursor' or turbineName == 'noTurbine':
            self.turbineFileName = None
            self.turbineDir = turbineDir
        else:
            if turbineDir is None:
                self.turbineDir = "./"
            else:
                self.turbineDir = turbineDir
            self.turbinePropDir = self.turbineDir + '/constant/turbineProperties/'
            self.turbineArrayProp = self.turbineDir + '/constant/turbineArrayProperties'
            if turbineFileName is None:
                self.turbineFileName = 'DTU10MW_POLIMI_WTM'
            else:
                self.turbineFileName = turbineFileName

            # Read turbine properties
            with open(self.turbinePropDir + self.turbineFileName, "r") as turbProp:
                flag = False
                blade_data = []
                numeric_pattern = '[+-]?\d*\.\d+ [+-]?\d*\.\d+ [+-]?\d*\.\d+ [+-]?\d*\.\d+ [+-]?\d*\.\d+ [+-]?\d+ '
                for line in turbProp:
                    if 'TipRad' in line:
                        self.rotor_R = float(re.findall('\d+\.?\d*', line)[0])
                        self.rotor_D = 2 * self.rotor_R
                    elif 'HubRad' in line:
                        self.hub_R = float(re.findall('\d+\.?\d*', line)[0])
                    elif 'TowerHt' in line:
                        self.tower_H = float(re.findall('\d+\.?\d*', line)[0])
                    elif 'BladeData' in line:
                        flag = True
                    elif 'TowerData' in line:
                        flag = False
                    if flag:
                        blade_data.append(re.findall(numeric_pattern, line))
            # remove empty elements from blade_data
            blade_data = list(filter(None, blade_data))
            # split sublists
            for i in range(0, len(blade_data)):
                dummy = [el.split() for el in blade_data[i]]
                blade_data[i] = dummy
            # convert to numpy array
            blade_data = np.array(blade_data, dtype=float).reshape(len(blade_data), 6)
            self.blade_r = blade_data[:, 0]
            self.blade_c = blade_data[:, 1]
            self.blade_twist = blade_data[:, 2]
            # self.nacelle_H = self.tower_H+self.offsetNacelle

            with open(self.turbineArrayProp, "r") as turbArr:
                for line in turbArr:
                    if 'baseLocation' in line:
                        self.turbX0 = float(re.findall('\d+\.?\d*', line)[0])
                        self.turbY0 = float(re.findall('\d+\.?\d*', line)[1])
                        self.turbZ0 = float(re.findall('\d+\.?\d*', line)[2])
                        self.nacelle_H = self.tower_H + self.turbZ0
                    elif 'numBladePoints' in line:
                        self.numBladePoints = int(re.findall('\d+\.?\d*', line)[0])
                    elif 'numNacellePoints' in line:
                        self.numNacellePoints = int(re.findall('\d+\.?\d*', line)[0])
                    elif 'numTowerPoints' in line:
                        self.numTowerPoints = int(re.findall('\d+\.?\d*', line)[0])
                    elif 'bladeForceProjectionType' in line:
                        self.bladeForceProjectionType = line.split()[1].translate(
                            str.maketrans('', '', string.punctuation))
                    elif 'bladeEpsilon' in line:
                        self.bladeepsilon0 = float(re.findall('\d+\.?\d*', line)[0])
                        self.bladeepsilon1 = float(re.findall('\d+\.?\d*', line)[1])
                        self.bladeepsilon2 = float(re.findall('\d+\.?\d*', line)[2])
                    elif 'nacelleEpsilon' in line:
                        self.nacelleepsilon0 = float(re.findall('\d+\.?\d*', line)[0])
                        self.nacelleepsilon1 = float(re.findall('\d+\.?\d*', line)[1])
                        self.nacelleepsilon2 = float(re.findall('\d+\.?\d*', line)[2])
                    elif 'towerEpsilon' in line:
                        self.towerepsilon0 = float(re.findall('\d+\.?\d*', line)[0])
                        self.towerepsilon1 = float(re.findall('\d+\.?\d*', line)[1])
                        self.towerepsilon2 = float(re.findall('\d+\.?\d*', line)[2])

        with open(self.turbineDir + '/setUp', "r") as setUp:
            for line in setUp:
                if 'xMin' in line:
                    self.xMin = float(re.findall('\d+\.?\d*', line)[0])
                elif 'yMin' in line:
                    self.yMin = float(re.findall('\d+\.?\d*', line)[0])
                elif 'zMin' in line:
                    self.zMin = float(re.findall('\d+\.?\d*', line)[0])
                elif 'xMax' in line:
                    self.xMax = float(re.findall('\d+\.?\d*', line)[0])
                elif 'yMax' in line:
                    self.yMax = float(re.findall('\d+\.?\d*', line)[0])
                elif 'zMax' in line:
                    self.zMax = float(re.findall('\d+\.?\d*', line)[0])
                elif 'U0Mag' in line:
                    self.U0Mag = float(re.findall('\d+\.?\d*', line)[0])
        self.XDomain = self.xMax - self.xMin
        self.YDomain = self.yMax - self.yMin
        self.ZDomain = self.zMax - self.zMin


class OpenFOAM(Turbine):
    """
    Class to postprocess OpenFOAM related output files (e.g. probes)

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

        TODO: change self attributes to return ones

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
        self.ux = self.Ux - self.UMeanx
        self.uy = self.Uy - self.UMeany
        self.uz = self.Uz - self.UMeanz

    def readWakeExperiment(self, expDir=None, probeSet='cross'):
        """
        Read velocity data from experimental data

        :param str expDir: experimental data file directory path
        :param str probeSet: probe's name
        """

        if expDir is None:
            expDir = '/home/giordi/Desktop/File_galleria/POLIMI_UNAFLOW_DATA'

        # WAKE CROSS
        if probeSet == 'cross':
            self.xCross = np.array([5.48, 5.48])
            self.yCross = np.arange(-1.47, 1.740, 0.1)
            self.zCross = np.array([2.085, 1.735])

            wake_cross = loadmat(expDir + '/' + 'WAKE_CROSS/TN02.mat')
            self.time = wake_cross['t']
            self.U1x = wake_cross['u'][:, :, 0]  # probe 1
            self.U1y = wake_cross['v'][:, :, 0]  # probe 1
            self.U1z = wake_cross['w'][:, :, 0]  # probe 1
            self.U2x = wake_cross['u'][:, :, 1]  # probe 2
            self.U2y = wake_cross['v'][:, :, 1]  # probe 2
            self.U2z = wake_cross['w'][:, :, 1]  # probe 2

            self.U1Meanx = np.sum(self.U1x, axis=0) / len(self.U1x)
            self.U1Meany = np.sum(self.U1y, axis=0) / len(self.U1y)
            self.U1Meanz = np.sum(self.U1z, axis=0) / len(self.U1z)
            self.U2Meanx = np.sum(self.U2x, axis=0) / len(self.U2x)
            self.U2Meany = np.sum(self.U2y, axis=0) / len(self.U2y)
            self.U2Meanz = np.sum(self.U2z, axis=0) / len(self.U2z)
        elif probeSet == 'along':
            self.xAlong = np.arange(2.18, 5.81, 0.33)
            self.yAlong = np.array([0.7, 0.9])
            self.zAlong = np.array([2.1, 2.1])

            wake_along = loadmat(expDir + '/' + 'WAKE_ALONG/TN16.mat')

            self.time = wake_along['t']
            self.U1x = wake_along['u'][:, :, 0]  # probe 1
            self.U1y = wake_along['v'][:, :, 0]  # probe 1
            self.U1z = wake_along['w'][:, :, 0]  # probe 1
            self.U2x = wake_along['u'][:, :, 1]  # probe 2
            self.U2y = wake_along['v'][:, :, 1]  # probe 2
            self.U2z = wake_along['w'][:, :, 1]  # probe 2

            self.U1Meanx = np.sum(self.U1x, axis=0) / len(self.U1x)
            self.U1Meany = np.sum(self.U1y, axis=0) / len(self.U1y)
            self.U1Meanz = np.sum(self.U1z, axis=0) / len(self.U1z)
            self.U2Meanx = np.sum(self.U2x, axis=0) / len(self.U2x)
            self.U2Meany = np.sum(self.U2y, axis=0) / len(self.U2y)
            self.U2Meanz = np.sum(self.U2z, axis=0) / len(self.U2z)
            '''
            for i in range(0, wake_along['u'].size[1]):
                self.U1x = wake_along['u'][:, i, 0]  # probe 1
                self.U1y = wake_along['v'][:, i, 0]  # probe 1
                self.U1z = wake_along['w'][:, i, 0]  # probe 1
                self.U2x = wake_along['u'][:, i, 1]  # probe 2
                self.U2y = wake_along['v'][:, i, 1]  # probe 2
                self.U2z = wake_along['w'][:, i, 1]  # probe 2

                vars(self)['U1Meanx'+str(i)] = np.sum(self.U1x, axis=0) / len(self.U1x)
                vars(self)['U1Meany'+str(i)] = np.sum(self.U1y, axis=0) / len(self.U1y)
                vars(self)['U1Meanz'+str(i)] = np.sum(self.U1z, axis=0) / len(self.U1z)
                vars(self)['U2Meanx'+str(i)] = np.sum(self.U2x, axis=0) / len(self.U2x)
                vars(self)['U2Meany'+str(i)] = np.sum(self.U2y, axis=0) / len(self.U2y)
                vars(self)['U2Meanz'+str(i)] = np.sum(self.U2z, axis=0) / len(self.U2z)
            '''

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

    def plotWakeProfile(self, figID, plotDir=None, var='p', normVar=None, compareID=None):
        """
        Compute and plot wake profile

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
            if normVar is None:
                xVar = vars(self)[var][-1, probeStart:probeEnd]
                xlabel = getAxes(var)
                figName = '/' + var + str(figID)
                distance = round((x[0] / self.rotor_D) / 0.5) * 0.5
                title = 'Distance: ' + str(distance) + 'D'
            else:
                xVar = np.divide(vars(self)[var][-1, probeStart:probeEnd], getattr(normVar[0], normVar[1])[-1])
                xlabel = getAxes(var) + '/' + normVar[1]
                title = 'Distance: ' + str(distance) + 'D'
                figName = var + str(figID) + '_norm'
            label = self.turbineName
            if compareID is not None:
                if normVar is None:
                    xVar = vars(self)[var][-1, probeStart:probeEnd]
                    plotCompare(compareID, xVar, yVar, label, plotDir, figName)
                else:
                    xVar = np.divide(vars(self)[var][-1, probeStart:probeEnd], getattr(normVar[0], normVar[1])[-1])
                    plotCompare(compareID, xVar, yVar, var, plotDir, figName)
            else:
                if var.startswith('TI'):
                    xlim = [0.0, 0.1]
                if normVar is None:
                    plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
                else:
                    plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim, xlim, title)
            probeStart = probeEnd

            self.xVarWake = xVar
            self.yVarWake = yVar
            self.ylabelWake = ylabel

    def plotWakeExperiment(self, compareID, plotDir=None, expProbe='probe_exp_cross1'):
        """
        Plot wwke profile from experimental data

        :param int compareID: figure identification number for comparison
        :param str plotDir: plot saving directory path
        :param str expProbe: experimental probe name
        """
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        if expProbe == 'probe_exp_cross1':
            xVar = self.U1Meanx
            yVar = self.yCross / self.rotor_D
            label = 'Experimental Data'
            figName = '/expComparison_cross1_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)
        elif expProbe == 'probe_exp_cross2':
            xVar = self.U2Meanx
            yVar = self.yCross / self.rotor_D
            label = 'Experimental Data'
            figName = '/expComparison_cross2_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)
        elif expProbe == 'probe_exp_along1':
            xVar = self.U1Meanx
            yVar = self.xAlong / self.rotor_D
            label = 'Experimental Data'
            figName = '/expComparison_along1_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)
        elif expProbe == 'probe_exp_along2':
            xVar = self.U2Meanx
            yVar = self.xAlong / self.rotor_D
            label = 'Experimental Data'
            figName = '/expComparison_along2_' + str(compareID)
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)

        self.xVarWakeExp = xVar
        self.yVarWakeExp = yVar

    def plotWakeDeficitProfile(self, figID, plotDir=None, var='p', normVar=None, compareID=None):

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
        '''
        :param probeRef: [reference object, 'x variable', 'y variable']
        :param probeCompare: [comparative object, 'variable']
        '''
        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        xVar = (getattr(probeRef[0], probeRef[1]) - getattr(probeCompare[0], probeCompare[1])) / getattr(probeRef[0],
                                                                                                         probeRef[1])
        yVar = getattr(probeRef[0], probeRef[2])
        xlabel = r'$\epsilon$ '
        ylabel = getattr(probeCompare[0], 'ylabelWake')
        xMean = np.mean(xVar)
        label = labelCompare + ' - ' + r'$\epsilon_{mean} = $' + str(xMean.round(decimals=2))
        title = 'Relative error distribution'
        figName = 'wakeError_' + str(figID)
        if compareID is None:
            plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, title=title)
        else:
            plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotResiduals(self, plotDir=None, var='UMeanx'):

        if plotDir is None:
            plotDir = './postProcessing/' + self.probeName + '/plots/'
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        xVar = self.probeLoc[:, 2]
        yVar = self.time
        zVar = vars(self)[var]

        fig = plt.figure()
        # ax = Axes3D(fig)
        ax = plt.axes(projection='3d')
        col_map = plt.get_cmap('viridis')
        X, Y = np.meshgrid(xVar, yVar, sparse=True)
        ax.plot_surface(X, Y, zVar, cmap=col_map)
        ax.set_xlabel('Height [m]')
        ax.set_ylabel('Time [s]')
        ax.set_zlabel(getAxes(var))
        ax.view_init(20, 130)
        fig.savefig(plotDir + 'Residual_' + var + '.eps', format='eps', dpi=1200)
        plt.close()


class SOWFA(Turbine):

    def __init__(self, turbName, turbineDir=None, turbineFileName=None):
        Turbine.__init__(self, turbName, turbineDir, turbineFileName)

    def readTurbineOutput(self, turbineOutDir=None, turbineNumber=1, nTurbines=1, FAST='off'):

        if turbineOutDir is None:
            turbineOutDir = "./postProcessing/turbineOutput/0/"

        # Find files in the directory
        files = os.listdir(turbineOutDir)
        self.SOWFAturbine = []
        self.SOWFAnacelle = []
        self.SOWFAtower = []
        self.SOWFAblade = []

        # Only rotor files for SOWFA-FAST
        if FAST == 'on':
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
        Plot root and tip values over time
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


class FAST(Turbine):

    def __init__(self, turbName, turbineDir=None, turbineFileName=None):
        Turbine.__init__(self, turbName, turbineDir, turbineFileName)

    def readBladeProp(self, fastDir=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeTipDeflectionsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeRoot(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeRootPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotBladeOverTime(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)
            figID += 1

    def plotBladeOverSpan(self, figID, plotDir=None, window=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotRotorPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerBaseLoads(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerBaseLoadsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopDisplacements(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopDisplacementsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopLoads(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotTowerTopLoadsPSD(self, figID, plotDir=None, var='all', ylim=None, xlim=None, compareID=None):

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
                        pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                               title=title)
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
                pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotLogCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotLowSpeedShaft(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def plotHighSpeedShaft(self, figID, plotDir=None, var='all', ylim=True, xlim=None, compareID=None):

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
                        plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim,
                             title=title)
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
                plot(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, ylim=ylim, xlim=xlim, title=title)
            else:
                plotCompare(compareID, xVar, yVar, label, plotDir, figName)

    def generateStatistics(self, statFile=None):

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