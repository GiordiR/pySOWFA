import os
from scipy.io import loadmat
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

class post(object):

    def __init__(self, pathSets, probeName, case):
        self.pathSets = pathSets
        self.probeName = probeName
        self.case = case
        self.hubH = 2.09
        self.rotD = 2.37
        self.rotR = self.rotD/2

    def readSets(self, var):
        vec = ['Ux','Uy','Uz', 'UMeanx', 'UMeany', 'UMeanz']
        if var in vec:
            with open(self.pathSets+'/'+self.probeName+'_U_UMean.xy', "r") as file:
                db = pd.read_csv(file, sep='\s+', skiprows=0, header=None)
                values = db.to_numpy(dtype=float)
                self.x = values[:, 0]
                self.y = values[:, 1]
                self.z = values[:, 2]
                self.Ux = values[:, 3]
                self.Uy = values[:, 4]
                self.Uz = values[:, 5]
                self.UMeanx = values[:, 6]
                self.UMeany = values[:, 7]
                self.UMeanz = values[:, 8]
        elif var.startswith('UPrime2Mean'):
            with open(self.pathSets+'/'+self.probeName+'_UPrime2Mean.xy', "r") as file:
                db2 = pd.read_csv(file, sep='\s+',skiprows=0, header=None)
                values2 = db2.to_numpy(dtype=float)
                self.x = values2[:, 0]
                self.y = values2[:, 1]
                self.z = values2[:, 2]
                self.UPrime2Meanxx = values2[:, 3]
                self.UPrime2Meanxy = values2[:, 4]
                self.UPrime2Meanxz = values2[:, 5]
                self.UPrime2Meanyy = values2[:, 6]
                self.UPrime2Meanyz = values2[:, 7]
                self.UPrime2Meanzz = values2[:, 8]

    def readExp(self):
        self.expDir = '/home/giordi/Desktop/File_galleria/POLIMI_UNAFLOW_DATA'

        # WAKE CROSS
        self.xCross = np.array([5.48, 5.48])
        self.yCross = np.arange(-1.47, 1.740, 0.1)
        self.zCross = np.array([2.085, 1.735])

        wake_cross = loadmat(self.expDir+'/'+'WAKE_CROSS/TN02.mat')
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

    def getTurbulenceIntensity(self):
        self.TIx = np.sqrt(self.UPrime2Meanxx) / self.UMeanx

    def plotWind(self, figID, var, normVar=None):
        self.plotDir = './plot/' + self.probeName
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)

        label = self.case
        if normVar is None:
            xVar = vars(self)[var]
        else:
            xVar = np.divide(vars(self)[var], getattr(normVar[0], normVar[1]))
        yVar = self.z
        xVar = savgol_filter(xVar, 31, 2)
        plt.figure(figID)
        plt.title('Distance: -1D')
        plt.plot(xVar, yVar, label=label)
        plt.ylabel('z')
        if var.startswith('TI'):
            plt.xlabel('$TI_x$')
        else:
            plt.xlabel('$U_x/U_\infty$')
        plt.savefig(self.plotDir + '/' + var + '_vProfile.png', format='png')

    def plotWakeH(self, figID, var, normVar=None, xlim=None, compareID=None):
        self.plotDir = './plot/' + self.probeName
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)

        if normVar == 'Exp':
            xVar = np.divide(vars(self)[var], 4.0)
            yVar = np.divide(self.yCross,self.rotD)
            label = 'Exp. Data'
        else:
            label = self.case
            if normVar is None:
                xVar = vars(self)[var]
            else:
                xVar = np.divide(vars(self)[var], getattr(normVar[0], normVar[1]))
            yVar = np.divide(self.y, self.rotD)
        if compareID is None:
            plt.figure(figID)
            if normVar != 'Exp':
                xVar = savgol_filter(xVar, 31, 2)
            plt.plot(xVar, yVar, label=label)
            if xlim is None:
                #xlim = [min(vars(self)[var])-1, max(vars(self)[var])+1]
                if var.startswith('TI'):
                    xlim = [0.0, 0.5]
                else:
                    xlim = [0.2, 1.2]
            ylim = [-1, 1]
            plt.xlim(xlim[0], xlim[1])
            plt.ylim(ylim[0], ylim[1])
            plt.plot(xlim, [0.5, 0.5], '--g', linewidth=0.5)
            plt.plot(xlim, [0, 0], '--g', linewidth=0.5)
            plt.plot(xlim, [-0.5, -0.5], '--g', linewidth=0.5)
            distance = round(self.x[0]/self.rotD)
            if self.probeName.startswith('probe_wake_cross'):
                plt.title('Distance: 2.3D')
            else:
                # solo per multiple turbine
                # if distance >= 5.0:
                #    distance = distance - 5.0
                #    plt.title('Distance: ' + str(distance) + 'D - Turbine 2')
                # else:
                #    plt.title('Distance: ' + str(distance) + 'D - Turbine 1')
                plt.title('Distance: ' + str(distance) + 'D')
            plt.savefig(self.plotDir+'/'+var+'h_'+self.probeName[-2:]+'.png', format='png')
        else:
            plt.figure(compareID)
            if normVar != 'Exp':
                xVar = savgol_filter(xVar, 31, 2)
            plt.plot(xVar, yVar, label=label)
            plt.ylabel('y/D')
            if var.startswith('TI'):
                plt.xlabel('$TI_x$')
            else:
                plt.xlabel('$U_x/U_\infty$')
            plt.legend(loc='best')
            plt.savefig(self.plotDir + '/' + var + 'h_' + self.probeName[-2:] + '.png', format='png')

    def plotWakeV(self, figID, var, normVar=None, xlim=None, compareID=None):
        self.plotDir = './plot/' + self.probeName
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)
        if normVar is None:
            xVar = vars(self)[var]
        else:
            xVar = np.divide(vars(self)[var], getattr(normVar[0], normVar[1]))
        yVar = self.z
        if compareID is None:
            plt.figure(figID)
            xVar = savgol_filter(xVar, 31, 2)
            plt.plot(xVar, yVar, label=self.case)
            if xlim is None:
                if var.startswith('TI'):
                    xlim = [0.0, 0.5]
                else:
                    xlim = [0.2, 1.2]
            ylim = [0.0, 3.6]
            plt.xlim(xlim[0], xlim[1])
            plt.ylim(ylim[0], ylim[1])
            plt.plot(xlim, [0.905, 0.905], '--g', linewidth=0.5)
            plt.plot(xlim, [2.09, 2.09], '--g', linewidth=0.5)
            plt.plot(xlim, [3.275, 3.275], '--g', linewidth=0.5)
            distance = round(self.x[0]/self.rotD)
            if self.probeName.startswith('probe_wake_cross'):
                plt.title('Distance: 2.3D')
            else:
                # solo per multiple turbine
                # if distance >= 5.0:
                #    distance = distance - 5.0
                #    plt.title('Distance: ' + str(distance) + 'D - Turbine 2')
                # else:
                #    plt.title('Distance: ' + str(distance) + 'D - Turbine 1')
                plt.title('Distance: ' + str(distance) + 'D')
        else:
            plt.figure(compareID)
            xVar = savgol_filter(xVar, 31, 2)
            plt.plot(xVar, yVar, label=self.case)
            plt.legend(loc='best')
        plt.ylabel('z [m]')
        if var.startswith('TI'):
            plt.xlabel('$TI_x$')
        else:
            plt.xlabel('$U_x/U_\infty$')
        plt.savefig(self.plotDir + '/' + var + 'v_' + self.probeName[-2:] + '.png', format='png')