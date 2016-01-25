from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BsToPhiMuMu_MuOnia_data_2012D_v2'
config.General.workArea = 'crab3_BsToPhiMuMu_data_2012D_v2'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bstophimumu_Run2012.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = False
config.JobType.outputFiles = ['BsToPhiMuMu.root']

config.Data.inputDataset = '/MuOniaParked/Run2012D-22Jan2013-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
config.Data.runRange = '190456-208686' # '193093-194075'                                                                          

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = ''

config.Site.storageSite = 'T2_IN_TIFR'

