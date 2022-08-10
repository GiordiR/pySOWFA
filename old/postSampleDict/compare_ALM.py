from postProc import *

######################################
# 2
pathSets1 = './S_flex/sets/35'
case1 = '$\epsilon = 2\Delta_g$'

# U ref
pRef = post(pathSets1, 'probe_URef', case1)
pRef.readSets('UMeanx')

for i in range(1,3):
    probeName = 'probe_wake_cross'+str(i)
    pc = post(pathSets1, probeName, case1)
    pc.readSets('UMeanx')
    pc.readSets('UPrime2Meanxx')
    pc.readExp()
    pc.plotWakeH(i,'UMeanx',normVar=[pRef, 'UMeanx'])
    varExp = 'U'+str(i)+'Meanx'
    pc.plotWakeH(i, varExp, normVar='Exp', compareID=i)

for j in range(1,9):
    probeName = 'probe_wake_h'+str(j)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j+2, 'UMeanx',normVar=[pRef, 'UMeanx'])
    ph.plotWakeH(j+10, 'TIx')

for k in range(1,9):
    probeName = 'probe_wake_v'+str(k)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k+18, 'UMeanx',normVar=[pRef, 'UMeanx'])
    ph.plotWakeV(k+26, 'TIx')

######################################
# 125
pathSets2 = './S_125/sets/35'
case2 = '$\epsilon = 1.25\Delta_g$'

# U ref
pRef = post(pathSets2, 'probe_URef', case2)
pRef.readSets('UMeanx')

for i in range(1, 3):
    probeName = 'probe_wake_cross' + str(i)
    pc = post(pathSets2, probeName, case2)
    pc.readSets('UMeanx')
    pc.readSets('UPrime2Meanxx')
    pc.readExp()
    pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
    varExp = 'U' + str(i) + 'Meanx'

for j in range(1, 9):
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in range(1, 9):
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# chord
pathSets3 = './S_chord/sets/35'
case3 = '$\epsilon_{chord}$'

# U ref
pRef = post(pathSets3, 'probe_URef', case3)
pRef.readSets('UMeanx')

for i in range(1, 3):
    probeName = 'probe_wake_cross' + str(i)
    pc = post(pathSets3, probeName, case3)
    pc.readSets('UMeanx')
    pc.readSets('UPrime2Meanxx')
    pc.readExp()
    pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
    varExp = 'U' + str(i) + 'Meanx'

for j in range(1, 9):
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets3, probeName, case3)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in range(1, 9):
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets3, probeName, case3)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)