from postProc import *

######################################
# RIGID
pathSets1 = './WTM_TURBULENT_fine_Rated_rigid/postProcessing/sets/130'
case1 = 'Rigid'

# U ref
pRef = post(pathSets1, 'probe_URef', case1)
pRef.readSets('UMeanx')

# for i in range(1,3):
#     probeName = 'probe_wake_cross'+str(i)
#     pc = post(pathSets1, probeName, case1)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i,'UMeanx',normVar=[pRef, 'UMeanx'])
#     varExp = 'U'+str(i)+'Meanx'
#     pc.plotWakeH(i, varExp, normVar='Exp', compareID=i)

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h'+str(j)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j+2, 'UMeanx',normVar=[pRef, 'UMeanx'])
    ph.plotWakeH(j+10, 'TIx')

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v'+str(k)+'D'
    ph = post(pathSets1, probeName, case1)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k+18, 'UMeanx',normVar=[pRef, 'UMeanx'])
    ph.plotWakeV(k+26, 'TIx')

######################################
# FLEX
pathSets2 = './WTM_TURBULENT_fine_Rated_flex/postProcessing/sets/130'
case2 = 'Elastic Carbon'

# U ref
pRef = post(pathSets2, 'probe_URef', case2)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets2, probeName, case2)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets2, probeName, case2)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# FLEX 1
pathSets3 = './WTM_TURBULENT_fine_Rated_flex1/postProcessing/sets/130'
case3 = 'Elastic 1'

# U ref
pRef = post(pathSets3, 'probe_URef', case3)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets3, probeName, case3)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets3, probeName, case3)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets3, probeName, case3)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# FLEX 2
pathSets4 = './WTM_TURBULENT_fine_Rated_flex2/postProcessing/sets/130'
case4 = 'Elastic 2'

# U ref
pRef = post(pathSets4, 'probe_URef', case4)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets3, probeName, case3)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets4, probeName, case4)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets4, probeName, case4)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# FLEX 3
pathSets5 = './WTM_TURBULENT_fine_Rated_flex3/postProcessing/sets/130'
case5 = 'Elastic 3'

# U ref
pRef = post(pathSets5, 'probe_URef', case5)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets4, probeName, case3)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets5, probeName, case5)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets5, probeName, case5)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# FLEX 4
pathSets6 = './WTM_TURBULENT_fine_Rated_flex4/postProcessing/sets/130'
case6 = 'Elastic 4'

# U ref
pRef = post(pathSets6, 'probe_URef', case6)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets5, probeName, case3)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets6, probeName, case6)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets6, probeName, case6)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)

######################################
# FLEX 4
pathSets7 = './WTM_TURBULENT_fine_Rated_flex5/postProcessing/sets/130'
case7 = 'Elastic 5'

# U ref
pRef = post(pathSets7, 'probe_URef', case7)
pRef.readSets('UMeanx')

# for i in range(1, 3):
#     probeName = 'probe_wake_cross' + str(i)
#     pc = post(pathSets6, probeName, case3)
#     pc.readSets('UMeanx')
#     pc.readSets('UPrime2Meanxx')
#     pc.readExp()
#     pc.plotWakeH(i, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=i)
#     varExp = 'U' + str(i) + 'Meanx'

for j in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_h' + str(j) + 'D'
    ph = post(pathSets7, probeName, case7)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeH(j + 2, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=j+2)
    ph.plotWakeH(j + 10, 'TIx', compareID=j+10)

for k in [x for x in range(1,9) if x!=2]:
    probeName = 'probe_wake_v' + str(k) + 'D'
    ph = post(pathSets7, probeName, case7)
    ph.readSets('UMeanx')
    ph.readSets('UPrime2Meanxx')
    ph.getTurbulenceIntensity()
    ph.plotWakeV(k + 18, 'UMeanx',normVar=[pRef, 'UMeanx'], compareID=k+18)
    ph.plotWakeV(k + 26, 'TIx', compareID=k+26)