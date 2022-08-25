import matplotlib.pyplot as plt


def plot(figID, xVar, yVar, xLabel, yLabel, label, plotDir, figName, ylim=None, xlim=None, title=None):
    plt.figure(figID)
    plt.plot(xVar, yVar, label=label)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])
    if title is not None:
        plt.title(title, fontweight='bold')
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def plotUtils(figID, xVar, yVar, label=None):
    plt.figure(figID)
    plt.plot(xVar, yVar, '--g', linewidth=0.5, label=label)


def plotCompare(figID, xVar, yVar, label, plotDir, figName):
    plt.figure(figID)
    plt.plot(xVar, yVar, label=label)
    plt.legend(loc='center right', bbox_to_anchor=(1.1, 0.5))
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def pltLogLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, yLim=None, xLim=None, title=None):
    plt.figure(figID)
    plt.loglog(xVar, yVar, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xLim is not None:
        plt.xlim(xLim[0], xLim[1])
    if yLim is not None:
        plt.ylim(yLim[0], yLim[1])
    if title is not None:
        plt.title(title, fontweight='bold')
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def plotLogLogCompare(figID, xVar, yVar, label, plotDir, figName):
    plt.figure(figID)
    plt.loglog(xVar, yVar, label=label)
    plt.legend()
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def pltLog(figID, xVar, yVar, xlabel, ylabel, label, plotDir, figName, yLim=None, xLim=None, title=None):
    plt.figure(figID)
    plt.semilogy(xVar, yVar, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xLim is not None:
        plt.xlim(xLim[0], xLim[1])
    if yLim is not None:
        plt.ylim(yLim[0], yLim[1])
    if title is not None:
        plt.title(title, fontweight='bold')
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def plotLogCompare(figID, xVar, yVar, label, plotDir, figName):
    plt.figure(figID)
    plt.semilogy(xVar, yVar, label=label)
    plt.legend()
    plt.savefig(plotDir + '/' + figName + '.eps', format='eps', dpi=1200)


def plotLogUtils(figID, xVar, yVar, label=None):
    plt.figure(figID)
    plt.semilogy(xVar, yVar, '--g', linewidth=0.5, label=label)


def endPlot():
    plt.close('all')


def getTitle(varName):
    titleDict = {
        # OpenFOAM
        "p": "Pressure",
        "pMean": "Mean pressure",
        "ppMean": "Pressure correlation",
        "k": "Turbulen Kinetic Energy",
        "kMean": "Modeled Turbulent Kinetic Energy",
        "nuSgs": "Sub-grid viscosity",
        "Ux": "Instantaneous streamwise velocity",
        "Uy": "Instantaneous spanwise velocity",
        "Uz": "Instantaneous vertical velocity",
        "TI": "Turbulence Intensity",
        "TIx": "Streamwise turbulence intensity",
        "TIy": "Spanwise turbulence intensity",
        "TIz": "Vertical turbulence intensity",
        "UMeanx": "Mean streamwise velocity",
        "UMeany": "Mean spanwise velocity",
        "UMeanz": "Mean vertical velocity",
        "upMeanx": "Velocity-pressure correlation",
        "upMeany": "Velocity-pressure correlation",
        "upMeanz": "Velocity-pressure correlation",
        "omegax": "Streamwise vorticity",
        "omegay": "Spanwise vorticity",
        "omegaz": "Vertical vorticity",
        "uuMeanxx": "Mean velocity correlation",
        "uuMeanxy": "Mean velocity correlation",
        "uuMeanxz": "Mean velocity correlation",
        "uuMeanyy": "Mean velocity correlation",
        "uuMeanyz": "Mean velocity correlation",
        "uuMeanzz": "Mean velocity correlation",
        "RMeanxx": "Modeled turbulent stress",
        "RMeanxy": "Modeled turbulent stress",
        "RMeanxz": "Modeled turbulent stress",
        "RMeanyy": "Modeled turbulent stress",
        "RMeanyz": "Modeled turbulent stress",
        "RMeanzz": "Modeled turbulent stress",
        "uuRTotalxx": "Total turbulent stress",
        "uuRTotalxy": "Total turbulent stress",
        "uuRTotalxz": "Total turbulent stress",
        "uuRTotalyy": "Total turbulent stress",
        "uuRTotalyz": "Total turbulent stress",
        "uuRTotalzz": "Total turbulent stress",
        "UPrime2Meanxx": "Root mean squared velocity",
        "UPrime2Meanyy": "Root mean squared velocity",
        "UPrime2Meanzz": "Root mean squared velocity",
        # Global Aerodyn
        "B1Azimuth": "Azimuth angle - Blade 1",
        "B2Azimuth": "Azimuth angle - Blade 2",
        "B3Azimuth": "Azimuth angle - Blade 3",
        "B1Pitch": "Pitch angle - Blade 1",
        "B2Pitch": "Pitch angle - Blade 2",
        "B3Pitch": "Pitch angle - Blade 3",
        "RtSpeed": "Rotor speed",
        "RtTSR": "Rotor tip-speed ratio",
        "RtAeroFxh": "Total rotor aerodynamic load",
        "RtAeroMxh": "Total rotor aerodynamic load",
        "RtAeroPwr": "Rotor aerodynamic power",
        "RtAeroCp": "Rotor aerodynamic power coefficient",
        "RtAeroCq": "Rotor aerodynamic torque coefficient",
        "RtAeroCt": "Rotor aerodynamic thrust coefficient",
        # Local Aerodyn
        "bladeAlpha": "Angle of Attack",
        "bladeRe": "Reynolds Number",
        "bladeM": "Mach Number",
        "bladeCl": "Lift force coefficient",
        "bladeCd": "Drag force coefficient",
        "bladeCm": "Pitching moment coefficient",
        # Globa ElastoDyn
        "BlPitch1": "Pitch angle - Blade1",
        "BlPitch2": "Pitch angle - Blade2",
        "BlPitch3": "Pitch angle - Blade3",
        "TipDxc1": "Out-of-plane blade tip deflection",
        "TipDyc1": "In-plane blade tip deflection",
        "TipDzc1": "Axial blade tip deflection",
        "TipDxb1": "Flapwise blade tip deflection",
        "TipDyb1": "Edgewise blade tip deflection",
        "RootFxc1": "Out-of-plane blade root shear force",
        "RootFyc1": "In-plane blade root shear force",
        "RootFxb1": "Flapwise blade root shear force",
        "RootFyb1": "Edgewise blade root shear force",
        "RootMxc1": "In-plane blade root moment",
        "RootMyc1": "Out-of-plane blade root moment",
        "RootMxb1": "Edgewise blade root moment",
        "RootMyb1": "Flapwise blade root moment",
        "RootMzb1": "Pitching moment at blade root",
        "LSShftFxa": "Low-speed shaft Thrust force (Rotor Thrust force)",
        "LSShftMxa": "Low-speed shaft Torque (Rotor Torque)",
        "RotPwr": "Low-speed shaft Power (Rotor Power)",
        "HSShftPwr": "High-speed shaft Power",
        "GenSpeed": "Angular speed of the high-speed shaft and generator",
        "TwrBsFxt": "Tower base fore-aft shear force",
        "TwrBsFyt": "Tower base side-to-side shear force",
        "TwrBsMxt": "Tower base side-to-side moment",
        "TwrBsMyt": "Tower base fore-aft moment",
        "TwrBsMzt": "Tower base torsional moment",
        "TTDspFA": "Tower top fore-aft deflection",
        "TTDspSS": "Tower top side-to-side deflection",
        "YawBrFxp": "Tower top fore-aft shear force",
        "YawBrFyp": "Tower top side-to-side shear force",
        # SOWFA
        "bladePointAxialForce": "Blade axial force distribution",
        "bladePointHorizontalForce": "Blade horizontal force distribution",
        "bladePointTorqueForce": "Blade torque distribution",
        "bladePointVaxial": "Blade axial velocity distribution",
        "bladePointVerticalForce": "Blade vertical force distribution",
        "bladePointVmag": "Blade velocity magnitude distribution",
        "bladePointVradial": "Blade radial velocity distribution",
        "bladePointVtangential": "Blade tangential velocity distribution",
        "bladePointX": "Blade X distribution",
        "bladePointY": "Blade Y distribution",
        "bladePointZ": "Blade Z distribution",
        "nacelleAxialForce": "Nacelle axial force",
        "nacelleHorizontalForce": "Nacelle horizontal force",
        "nacellePointAxialForce": "Nacelle axial force distribution",
        "nacellePointDrag": "Nacelle drag distribution",
        "nacellePointHorizontalForce": "Nacelle horizontal force distribution",
        "nacellePointVaxial": "Nacelle axial velocity distribution",
        "nacellePointVerticalForce": "Nacelle vertical force distribution",
        "nacellePointVhorizontal": "Nacelle horizontal velocity distribution",
        "nacellePointVmag": "Nacelle velocity magnitude distribution",
        "nacellePointVvertical": "Nacelle vertical velocity distribution",
        "nacelleVerticalForce": "Nacelle vertical force",
        "rotorAxialForce": "Rotor axial force",
        "rotorHorizontalForce": "Rotor horizontal force",
        "rotorPower": "Rotor power",
        "rotorSpeed": "Rotor speed",
        "rotorTorque": "Rotor torque",
        "rotorVerticalForce": "Rotor vertical force",
        "towerAxialForce": "Tower axial force",
        "towerHorizontaklForce": "Tower horizontal force",
        "towerPointAlpha": "Tower alpha distribution",
        "towerPointAxialForce": "Tower axial force distribution",
        "towerPointDrag": "Tower drag distribution",
        "towerPointHorizontalForce": "Tower horizontal force distribution",
        "towerPointLift": "Tower lift distribution",
        "towerPointVaxial": "Tower axial velocity distribution",
        "towerPointVerticalForce": "Tower vertical force distribution",
        "towerPointVhorizontal": "Tower horizontal velocity distribution",
        "towerPointVmag": "Tower velocity magnitude distribution",
        "towerPointVvertical": "Tower vertical velocity distribution",
                }
    newName = titleDict[varName]
    return newName


def getAxes(varName):
    axesDict = {
        # OpenFOAM
        "Time": "Time [s]",
        "time": "Time [s]",
        "p": r"$p$",
        "pMean": r"$\langle p \rangle$",
        "ppMean": r"$\langle pp \rangle$",
        "k": r"$k$",
        "kMean": r"$\langle k \rangle$",
        "nuSgs": r"$\nu_{SGS}$",
        "Ux": r"$U$",
        "Uy": r"$V$",
        "Uz": r"$W$",
        "TI": r"$TI$",
        "TIx": r"$TI_x$",
        "TIy": r"$TI_y$",
        "TIz": r"$TI_z$",
        "UMeanx": r"$\langle U \rangle$",
        "UMeany": r"$\langle V \rangle$",
        "UMeanz": r"$\langle W \rangle$",
        "upMeanx": r"$\langle up \rangle$",
        "upMeany": r"$\langle vp \rangle$",
        "upMeanz": r"$\langle wp \rangle$",
        "omegax": r"$\omega_x$",
        "omegay": r"$\omega_y$",
        "omegaz": r"$\omega_z$",
        "uuMeanxx": r"$\langle uu \rangle$",
        "uuMeanxy": r"$\langle uv \rangle$",
        "uuMeanxz": r"$\langle uw \rangle$",
        "uuMeanyy": r"$\langle vv \rangle$",
        "uuMeanyz": r"$\langle vw \rangle$",
        "uuMeanzz": r"$\langle ww \rangle$",
        "RMeanxx": r"$\langle \tau_{xx} \rangle$",
        "RMeanxy": r"$\langle \tau_{xy} \rangle$",
        "RMeanxz": r"$\langle \tau_{xz} \rangle$",
        "RMeanyy": r"$\langle \tau_{yy} \rangle$",
        "RMeanyz": r"$\langle \tau_{yz} \rangle$",
        "RMeanzz": r"$\langle \tau_{zz} \rangle$",
        "uuRTotalxx": r"$\langle uu \rangle + \langle \tau_{xx} \rangle$",
        "uuRTotalxy": r"$\langle uv \rangle + \langle \tau_{xy} \rangle$",
        "uuRTotalxz": r"$\langle uw \rangle + \langle \tau_{xz} \rangle$",
        "uuRTotalyy": r"$\langle vv \rangle + \langle \tau_{yy} \rangle$",
        "uuRTotalyz": r"$\langle vw \rangle + \langle \tau_{yz} \rangle$",
        "uuRTotalzz": r"$\langle ww \rangle + \langle \tau_{zz} \rangle$",
        "UPrime2Meanxx": r"$U_{rms}$",
        "UPrime2Meanxz": r"$\langle u'w' \rangle$",
        "UPrime2Meanyy": r"$V_{rms}$",
        "UPrime2Meanzz": r"$W_{rms}$",
        # Global Aerodyn
        "B1Azimuth": "Azimuth angle [deg]",
        "B2Azimuth": "Azimuth angle [deg]",
        "B3Azimuth": "Azimuth angle [deg]",
        "B1Pitch": "Pitch angle [deg]",
        "B2Pitch": "Pitch angle [deg]",
        "B3Pitch": "Pitch angle [deg]",
        "RtSpeed": "Rotor speed [rpm]",
        "RtTSR": "Rotor tip-speed ratio",
        "RtAeroFxh": "Fx [N]",
        "RtAeroMxh": "Mx [N-m]",
        "RtAeroPwr": "Power [W]",
        "RtAeroCp": r"$C_p$ ",
        "RtAeroCq": r"$C_q$ ",
        "RtAeroCt": r"$C_t$ ",
        # Local Aerodyn
        "bladeAlpha": r"$\alpha$ [deg]",
        "bladeRe": "Re",
        "bladeM": "Mach",
        "bladeCl": r"$C_l$",
        "bladeCd": r"$C_d$",
        "bladeCm": r"$C_m$",
        # Globa ElastoDyn
        "BlPitch1": "Pitch angle [deg]",
        "BlPitch2": "Pitch angle [deg]",
        "BlPitch3": "Pitch angle [deg]",
        "TipDxc1": "Oop tip deflection [m]",
        "TipDyc1": "Ip tip deflection [m]",
        "TipDzc1": "Axial tip deflection [m]",
        "TipDxb1": "Flapwise tip deflection [m]",
        "TipDyb1": "Edgewise tip deflection [m]",
        "RootFxc1": "Out-of-plane shear force [kN]",
        "RootFyc1": "In-plane shear force [kN]",
        "RootFxb1": "Flapwise shear force [kN]",
        "RootFyb1": "Edgewise shear force [kN]",
        "RootMxc1": "In-plane moment [kN-m]",
        "RootMyc1": "Out-of-plane moment [kN-m]",
        "RootMxb1": "Edgewise moment [kN-m]",
        "RootMyb1": "Flapwise moment [kN-m]",
        "RootMzb1": "Pitching moment [kN-m]",
        "LSShftFxa": "Thrust force [kN]",
        "LSShftMxa": "Torque [kN-m]",
        "RotPwr": "Power [kW]",
        "HSShftPwr": "Power [kW]",
        "GenSpeed": "Generator speed [rpm]",
        "TwrBsFxt": "Fore-aft shear force [kN]",
        "TwrBsFyt": "Side-to-side shear force [kN]",
        "TwrBsMxt": "Side-to-side moment [kN-m]",
        "TwrBsMyt": "Fore-aft moment [kN-m]",
        "TwrBsMzt": "Torsional moment [kN-m]",
        "TTDspFA": "Fore-aft deflection [m]",
        "TTDspSS": "Side-to-side deflection [m]",
        "YawBrFxp": "Fore-aft shear force [kN]",
        "YawBrFyp": "Side-to-side shear force [kN]",
        # SOWFA
        "bladePointAxialForce": "Blade axial force [N]",
        "bladePointHorizontalForce": "Blade horizontal force [N]",
        "bladePointTorqueForce": "Blade torque [N-m]",
        "bladePointVaxial": "Blade axial velocity [m/s]",
        "bladePointVerticalForce": "Blade vertical force [N]",
        "bladePointVmag": "Blade velocity magnitude [m/s]]",
        "bladePointVradial": "Blade radial velocity [m/s]]",
        "bladePointVtangential": "Blade tangential velocity [m/s]]",
        "bladePointX": "Blade X [m]",
        "bladePointY": "Blade Y [m]]",
        "bladePointZ": "Blade Z [m]]",
        "nacelleAxialForce": "Nacelle axial force [N]",
        "nacelleHorizontalForce": "Nacelle horizontal force [N]",
        "nacellePointAxialForce": "Nacelle axial force [N]",
        "nacellePointDrag": "Nacelle drag",
        "nacellePointHorizontalForce": "Nacelle horizontal force [N]",
        "nacellePointVaxial": "Nacelle axial velocity [m/s]",
        "nacellePointVerticalForce": "Nacelle vertical force [N]",
        "nacellePointVhorizontal": "Nacelle horizontal velocity [m/s]",
        "nacellePointVmag": "Nacelle velocity magnitude [m/s]",
        "nacellePointVvertical": "Nacelle vertical velocity [m/s]",
        "nacelleVerticalForce": "Nacelle vertical force [N]",
        "rotorAxialForce": "Rotor axial force [N]",
        "rotorHorizontalForce": "Rotor horizontal force [N]",
        "rotorPower": "Rotor power [W]",
        "rotorSpeed": "Rotor speed [rpm]",
        "rotorTorque": "Rotor torque [N-m]",
        "rotorVerticalForce": "Rotor vertical [N]",
        "towerAxialForce": "Tower axial [N]",
        "towerHorizontaklForce": "Tower horizontal [N]",
        "towerPointAlpha": "Tower alpha [deg]",
        "towerPointAxialForce": "Tower axial force [N]",
        "towerPointDrag": "Tower drag",
        "towerPointHorizontalForce": "Tower horizontal force [N]",
        "towerPointLift": "Tower lift",
        "towerPointVaxial": "Tower axial velocity [m/s]",
        "towerPointVerticalForce": "Tower vertical force [N]",
        "towerPointVhorizontal": "Tower horizontal velocity [m/s]",
        "towerPointVmag": "Tower velocity magnitude [m/s]",
        "towerPointVvertical": "Tower vertical velocity [m/s]",
        }
    newName = axesDict[varName]
    return newName
