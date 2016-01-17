import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring('cerr', 'cout', 'message'),
    categories = cms.untracked.vstring('myHLT', 'myBeam', 'myBs'),
    cerr = cms.untracked.PSet(threshold = cms.untracked.string('WARNING')),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myHLT = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myBeam = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        myBs = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 

    ), 
     message = cms.untracked.PSet(
         threshold = cms.untracked.string('INFO'),
         INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
         myHLT = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
         myBeam = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
         myBs = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
     )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff') 
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# add track candidates
from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,
                    label        = 'TrackCands',                  
                    tracks       = cms.InputTag('generalTracks'), 
                    particleType = 'pi+',                         
                    preselection = 'pt > 0.1',                     
                    selection    = 'pt > 0.1',                     
                    isolation    = {},                            
                    isoDeposits  = [],                            
                    mcAs         = None          
)    

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'], outputModules=[])



process.ntuple = cms.EDAnalyzer(
    'BsToPhiMuMu',

    OutputFileName = cms.string("BsToPhiMuMu.root"),
    BuildBsToPhiMuMu = cms.untracked.bool(True), 

    MuonMass = cms.untracked.double(0.10565837), 
    MuonMassErr = cms.untracked.double(3.5e-9),   
    KaonMass = cms.untracked.double(0.493677), 
    KaonMassErr = cms.untracked.double(1.6e-5),
    BsMass = cms.untracked.double(5.36677),          ## put the Bs Mass (pdg value)

    # labels
    GenParticlesLabel = cms.InputTag("genParticles"),
    TriggerResultsLabel = cms.InputTag("TriggerResults","", 'HLT'),
    BeamSpotLabel = cms.InputTag('offlineBeamSpot'),
    VertexLabel = cms.InputTag('offlinePrimaryVertices'),
    MuonLabel = cms.InputTag('cleanPatMuonsTriggerMatch'),
    TrackLabel = cms.InputTag('cleanPatTrackCands'), 
    TriggerNames = cms.vstring([]),
    LastFilterNames = cms.vstring([]),

    # gen particle 
    IsMonteCarlo = cms.untracked.bool(False),
    KeepGENOnly  = cms.untracked.bool(False),
    TruthMatchMuonMaxR = cms.untracked.double(0.004), # [eta-phi]
    TruthMatchKaonMaxR = cms.untracked.double(0.3), # [eta-phi]


    # HLT-trigger cuts (for reference https://espace.cern.ch/cms-quarkonia/trigger-bph/SitePages/2012-LowMass.aspx)
    MuonMinPt = cms.untracked.double(3.5), # 3.0 [GeV]
    MuonMaxEta = cms.untracked.double(2.2),  
    MuonMaxDcaBs = cms.untracked.double(2.0), # [cm]

    MuMuMinPt = cms.untracked.double(6.9),      # [GeV/c]
    MuMuMinInvMass = cms.untracked.double(1.0), # [GeV/c2]
    MuMuMaxInvMass = cms.untracked.double(4.8), # [GeV/c2]

    MuMuMinVtxCl = cms.untracked.double(0.10), # 0.05
    MuMuMinLxySigmaBs = cms.untracked.double(3.0), 
    MuMuMaxDca = cms.untracked.double(0.5), # [cm]
    MuMuMinCosAlphaBs = cms.untracked.double(0.9),

    # pre-selection cuts 
    TrkMinPt = cms.untracked.double(0.4), # 0.4 [GeV/c]
    TrkMinDcaSigBs = cms.untracked.double(0.8), # 0.8 hadron DCA/sigma w/respect to BS (=>changed Max to Min)
    TrkMaxR = cms.untracked.double(110.0), # [cm] ==> size of tracker volume in radial direction
    TrkMaxZ = cms.untracked.double(280.0), # [cm] ==> size of tracker volume in Z direction

    ## phi(1020) mass = 1019.461 +/- 0.019 MeV, full width = 4.266 +/- 0.031 MeV
    PhiMinMass = cms.untracked.double(1.00), # [GeV/c2]  - 3 sigma of the width(~5MeV)
    PhiMaxMass = cms.untracked.double(1.04), # [GeV/c2]  + 3 sigma of the width

    BsMinVtxCl = cms.untracked.double(0.01), 
    BsMinMass = cms.untracked.double(4.5), # [GeV/c2] 
    BsMaxMass = cms.untracked.double(6.5), # [GeV/c2]  

)


# Remove not used from PAT 
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
#process.patDefaultSequence.remove(process.patJetPartonAssociation) #Uncomment this line for SW version <= 5_3_12
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.p = cms.Path(process.patDefaultSequence * process.ntuple)
