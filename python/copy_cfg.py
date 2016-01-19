import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring (
#"/store/mc/Summer12_DR53X/BuToKstarMuMu_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/006C6D81-8B4D-E311-BEBC-E0CB4E1A1194.root")                   
'/store/mc/Summer12_DR53X/BsToJPsiPhi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/20000/0608FB5D-C477-E211-A453-002481E0DDBE.root')
                             
)

# tell the process to only run over 500 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32 (100)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string ("BsToPhiMuMu_MC_8TeV_aod_100.root")                       
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
