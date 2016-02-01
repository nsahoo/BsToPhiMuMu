from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Ntuple_BsToPhiMuMu_8TeV_MC_v3'
config.General.workArea = 'crab3_BsToPhiMuMu_MC_v3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bstophimumu_MC.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = False
config.JobType.outputFiles = ['BsToPhiMuMu.root']

config.Data.inputDataset = '/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/Summer12_DR53X-PU_RD2_START53_V19F-v3/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.lumiMask = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
#config.Data.runRange = '190456-208686' # '193093-194075'                                                                    
                                                                                                                                                 

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = ''

config.Site.storageSite = 'T2_IN_TIFR'

