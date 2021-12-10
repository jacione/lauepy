
"""
A whole bunch of PV names used for cgi web display.

-- Ruqing Xu, Apr. 2012
"""

# server ip addresses
ip34ide = '164.54.125.207'
ip34ida = '164.54.125.208'
ip34idEPS = '164.54.125.40' #'164.54.125.11'
ipPE1 = '164.54.125.121'
ipPE2 = '164.54.125.122'
ipPE3 = '164.54.125.123'

# name prefixes
prefix34idEPS = '34idPLC'
prefix34ide = '34ide'
prefix34ida = '34ida'
prefixPE1 = '34idePE1'
prefixPE2 = '34idePE2'
prefixPE3 = '34idePE3'

# PVs for status tests
status34ide = '34ide:userName' #'34ide:Status'
status34idEPS = '34idPLC:eps:Y0' 

# Scan PVs
ide34scanpf = '34ide:scan'
ide34scanbusys = ['34ide:scan1.BUSY', '34ide:scan2.BUSY', \
                  '34ide:scan3.BUSY', '34ide:scan4.BUSY']
ide34scanPaused = '34ide:scanPause.VAL'
ide34scanDimension = '34ide:ScanDim.VAL'
ide34scanPtsNos = ['34ide:scan1.NPTS', '34ide:scan2.NPTS', \
                   '34ide:scan3.NPTS', '34ide:scan4.NPTS']

ide34scan1Mode = '34ide:scan1.P1SM'
ide34scan1TriggerPVs = ['34ide:scan1.T1PV', '34ide:scan1.T2PV', \
                        '34ide:scan1.T3PV', '34ide:scan1.T4PV']

ide34scanPercent = '34ide:scanProgress:percentDone'
ide34scanStartTimeStr = '34ide:scanProgress:startingTimeStr'
ide34scanEndTimeStr = '34ide:scanProgress:endingTimeStr'
ide34scanElapsedTimeStr = '34ide:scanProgress:totalElapsedTimeStr'
ide34scanPauseTimeStr = '34ide:scanProgress:pauseTimeStr'
ide34scanRemainingTimeStr = '34ide:scanProgress:remainingTimeStr'

# 34ide PVs
ide34monoMode = '34ide:monoT:transStatus'
ide34monoEnergy = '34ide:monoE:keV.RBV'
ide34monoIDSync = '34ide:UndSync.VAL' #'34ide:bit19.VAL'
ide34beambad = '34ide:beamBad.VAL' #'34ide:userCalc2.VAL'

# other 34id PVs
DShutterClosed = '34ide:beamBad.D' #'34ide:userCalc2.D'
FEShutterClosed = '34ide:beamBad.C' #'34ide:userCalc2.C'
ID34dsControlMode = '34ide:beamBad.F' #'34ide:userCalc2.F'
ID34dsGapmm = '34ide:beamBad.B' #'34ide:userCalc2.B'
ID34dsTapermm = '34ide:beamBad.H' #'34ide:userCalc2.H'
ID34dsHarmonic = '34ide:undSyncScan.D' #'34ide:userTran6.D'
ID34dsEnergy = '34ide:undSyncScan.M' #'34ide:userTran6.I'

# Storage ring info PVs
ringCurrent = '34ide:beamBad.A' #'34ide:userCalc2.A'

# detector PVs
scalerTrigger = '34ide:vsc:c0.CNT'
PEcommonStr = '34idePE'
PEfileNs = ['34idePE1:fileROI0:FileNumber_RBV', \
            '34idePE2:Nexus1:FileNumber_RBV',   \
            '34idePE3:fileROI0:FileNumber_RBV' ]

# EPS PVs
eps34iddOK = '34idPLC:eps:Y3'
eps34idVacOK = '34idPLC:eps:Y0'
eps34idaOK = '34idPLC:eps:Y1'
eps34ideHeOK = '34idPLC:eps:X46'

# user custom watch PVs
userPVstrings = ['34ide:userSeq4.DOL1','34ide:userSeq4.DOL2','34ide:userSeq4.DOL3']
#userPVenables = ['34ide:userCalc4.INAP','34ide:userCalc4.INBP','34ide:userCalc4.INCP']
userPVvalids = ['34ide:userSeq4.DOL1V','34ide:userSeq4.DOL2V','34ide:userSeq4.DOL3V']
userPVvalues = ['34ide:userSeq4.STR1','34ide:userSeq4.STR2','34ide:userSeq4.STR3']
userPVrefresh = '34ide:userSeq4.PROC'
