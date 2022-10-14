import re
import numpy as np
import string


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

