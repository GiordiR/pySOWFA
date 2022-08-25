from postProc import *

######################################
#COARSE
pathSets1 = './S_coarse/sets/35'
case1 = 'Coarse mesh'

# U ref
pRef = post(pathSets1, 'probe_URef', case1)
pRef.readSets('UMeanx')

for i in range(1,3):
    probeName = 'probe_wake_cross'+str(i)
    pc = post(pathSets1, probeName, case1)
    pc.readSets('UMeanx')
    pc.readExp()
    pc.plotWakeH(i,'UMeanx', normVar=[pRef, 'UMeanx'])
    varExp = 'U'+str(i)+'Meanx'
    pc.plotWakeH(i, varExp, normVar='Exp', compareID=i)

for j in range(1,9):
    probeName = 'probe_wake_h'+str(j)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.plotWakeH(j+2, 'UMeanx', normVar=[pRef, 'UMeanx'])

for k in range(1,9):
    probeName = 'probe_wake_v'+str(k)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.plotWakeV(k+18, 'UMeanx', normVar=[pRef, 'UMeanx'])

######################################
# FINE
pathSets2 = './S_flex/sets/35'
case2 = 'Fine mesh'

# U ref
pRef = post(pathSets2, 'probe_URef', case2)
pRef.readSets('UMeanx')

for i in range(1, 3):
    probeName = 'probe_wake_cross' + str(i)
    pc = post(pathSets2, probeName, case2)
    pc.readSets('UMeanx')
    pc.readExp()
    pc.plotWakeH(i, 'UMeanx', normVar=[pRef, 'UMeanx'], compareID=i)
    varExp = 'U' + str(i) + 'Meanx'

for j in range(1, 9):
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.plotWakeH(j + 2, 'UMeanx', normVar=[pRef, 'UMeanx'], compareID=j+2)

for k in range(1, 9):
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.plotWakeV(k + 18, 'UMeanx', normVar=[pRef, 'UMeanx'], compareID=k+18)