// -*- C++ -*-
//
// Package:    BsToPhiMuMu
// Class:      BsToPhiMuMu
// 
/**\class BsToPhiMuMu BsToPhiMuMu.cc BpHaNA/BsToPhiMuMu/src/BsToPhiMuMu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//=====================================================================
// original author:  Niladribihari Sahoo,42 3-024,+41227662373,        |
//        copyright  @ N.Sahoo, NISER, Bhubaneswar                     |
//         created:  Sat Nov 28 07:33:44 CET 2015                      |
//         added saveGenInfo  (sat 16 jan 2016)                        |
//         added soft muon id info, trigger matching info (sat 16 jan) | 
//=====================================================================
// $Id$
//
//


//-----------------------                                                                                                                      
// system include files                                                                                                                             
//-----------------------                                                                                                                     
#include <memory>

//----------------------                                                                                                                       
// user include files                                                                                                                                 
//----------------------                                                                                                                         
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>


using namespace std;


const int MUONMINUS_PDG_ID = 13;
const int KAONPLUS_PDG_ID = 321;
const int PHI_PDG_ID = 333;      // phi(1020)
const int BS_PDG_ID = 531;
const int JPSI_PDG_ID = 443;
const int PSI2S_PDG_ID = 100443;

const double PI = 3.141592653589793;


//----------------
//  structures
//----------------
struct HistArgs{
  char name[128];
  char title[128];
  int n_bins;
  double x_min;
  double x_max;
};

enum HistName{
  h_events,
  h_mupt,
  h_mueta,
  h_mumdcabs,
  h_mumutrkr,
  h_mumutrkz,

  h_mumudca,
  h_mumuvtxcl,
  h_mumupt,
  h_mumumass,
  h_mumulxybs,

  h_mumucosalphabs,
  h_trkpt,
  h_trkdcasigbs,
  h_bsvtxchisq,
  h_bsvtxcl,

  h_phimass,
  h_bsmass,                                                                                                                           

  kHistNameSize
};

//--------------------                                                                                                                           
// Global hist args                                                                                                                               
//--------------------                                                                                                                         
HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max                                                                                                              

  {"h_events", "Processed Events", 1, 0, 1},
  {"h_mupt", "Muon pT; [GeV]", 100, 0, 30},
  {"h_mueta", "Muon eta", 100, 0, 3},
  {"h_mumdcabs", "#mu^{-} DCA beam spot; DCA [cm]", 100, 0, 10},
  {"h_mumutrkr", "#mu^{+}#mu^{-} distance in phi-eta; [cm]", 100, 0, 50},
  {"h_mumutrkz", "#mu^{+}#mu^{-} distance in Z; [cm]", 100, 0, 100},

  {"h_mumudca",  "#mu^{+}#mu^{-} DCA; [cm]", 100, 0, 20},
  {"h_mumuvtxcl",  "#mu^{+}#mu^{-} vertex CL", 100, 0, 1},
  {"h_mumupt",    "#mu^{+}#mu^{-} pT ; pT [GeV]", 100, 0, 50},
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV/c^{2}]",
   100, 2, 20},
  {"h_mumulxybs", "#mu^{+}#mu^{-} Lxy #sigma beam spot", 100, 0, 100},

  {"h_mumucosalphabs", "#mu^{+}#mu^{-} cos #alpha beam spot", 100, 0, 1},
  {"h_trkpt", "Pion track pT; pT [GeV]", 100, 0, 20},
  {"h_trkdcasigbs", "Pion track DCA/#sigma beam spot; DCA/#sigma", 1000, 0, 100},
  {"h_bsvtxchisq", "B_{s} decay vertex chisq", 100, 0, 1000},
  {"h_bsvtxcl", "#B_{s} decay vertex CL", 100, 0, 1},

  {"h_phimass", "#phi(1020) mass; M(KK) [GeV/^{2}]", 100, 0, 20},   
  {"h_bsmass", "B_{s} mass; M(KK#mu#mu}) [GeV]", 100, 0, 20},                                                                                           

};

//--------------------                                                                                                                                    
// Define histograms                                                                                                                                       
//--------------------                                                                                                                                    
TH1F *histos[kHistNameSize];


//-----------------------
// class declaration
//-----------------------

class BsToPhiMuMu : public edm::EDAnalyzer {
   public:
      explicit BsToPhiMuMu(const edm::ParameterSet&);
      ~BsToPhiMuMu();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);



  bool buildBsToPhiMuMu(const edm::Event &);

  void calLS (double, double, double, double, double, double, double,
              double, double,  double, double, double, double, double,
              double, double, double, double, double*, double*);

  void calCosAlpha (double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double,
                    double, double, double, double,
                    double*, double*);

  void calCosAlpha2d (double, double, double, double, double,
                      double, double, double, double, double,
                      double, double, double, double,
                      double, double, double, double,
                      double*, double*);

  void calCtau(RefCountedKinematicTree, double &, double &);
  double calEta(double, double, double);
  double calPhi(double, double, double);
  double calEtaPhiDistance (double, double, double, double, double, double);
  void clearVariables();

  bool hasBeamSpot(const edm::Event&);

  bool calClosestApproachTracks(const reco::TransientTrack,
                                const reco::TransientTrack,
                                double&, double &, double &);

  bool hasGoodPhiVertex(const reco::TransientTrack, const reco::TransientTrack, 
			reco::TransientTrack &, reco::TransientTrack &, 
			double &, double &);

  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint,
                             double, double &, double &);

  bool hasGoodBsMass(RefCountedKinematicTree, double &);  

  bool hasGoodBsVertex(const reco::TransientTrack, const reco::TransientTrack,
		       const reco::TransientTrack, const reco::TransientTrack,
		       double &, double &, double &, RefCountedKinematicTree &);

  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
                          reco::TransientTrack &, reco::TransientTrack &,
                          double &, double &, double &, double &, double &,
                          double &, double &, double &);

  bool hasGoodTrack(const edm::Event&, const pat::GenericParticle, double &);

  bool hasPrimaryVertex(const edm::Event &);

  void hltReport(const edm::Event&);

  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  bool matchMuonTracks (const edm::Event&, const vector<reco::TrackRef>);
  bool matchPrimaryVertexTracks ();

  void saveBsToPhiMuMu(const RefCountedKinematicTree);
  void saveBsVertex(RefCountedKinematicTree);
  void saveBsCosAlpha(RefCountedKinematicTree);
  void saveBsCosAlpha2d(RefCountedKinematicTree);
  void saveBsLsig(RefCountedKinematicTree);
  void saveBsCtau(RefCountedKinematicTree);

  void saveGenInfo(const edm::Event&);
  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef);
  void saveDimuVariables(double, double, double, double, double, double,
                         double, double, double, double, double, double,
                         double, double);
  void saveMuonTriggerMatches(const pat::Muon, const pat::Muon);
  void saveTruthMatch(const edm::Event& iEvent);


      // ----------member data ---------------------------


  // --- begin input from python file ---                                                                                                           
  string OutputFileName_;
  bool BuildBsToPhiMuMu_;

  //----------------------                                                                                                                                 
  // particle properties                                                                                                                                   
  //----------------------                                                                                                                                    
  ParticleMass MuonMass_;
  float MuonMassErr_;
  ParticleMass KaonMass_;
  float KaonMassErr_;
  ParticleMass PhiMass_;
  float PhiMassErr_;
  double BsMass_;

  //----------                                                                                                                                                
  // labels                                                                                                                                               
  //----------                                                                                                                                               
  edm::InputTag GenParticlesLabel_;
  edm::InputTag TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  edm::InputTag MuonLabel_;
  edm::InputTag LambdaLabel_;
  edm::InputTag TrackLabel_;
  vector<string> TriggerNames_;
  vector<string> LastFilterNames_;

  //---------------                                                                                                                                        
  // gen particle                                                                                                                                                  
  //---------------                                                                                                                                          
  bool   IsMonteCarlo_;
  bool   KeepGENOnly_;
  double TruthMatchMuonMaxR_;
  double TruthMatchKaonMaxR_;

  //---------------------                                                                                                                                           
  // pre-selection cuts                                                                                                                                   
  //---------------------                                                                                                                              
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;
  double TrkMinPt_;
  double TrkMinDcaSigBs_;
  double TrkMaxR_;
  double TrkMaxZ_;
  double MuMuMaxDca_;
  double MuMuMinVtxCl_;
  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;
  double MuMuMinLxySigmaBs_;
  double MuMuMinCosAlphaBs_;
  double PhiMinMass_;
  double PhiMaxMass_;
  double BsMinVtxCl_;
  double BsMinMass_;
  double BsMaxMass_;


  //--------------------                                                                                                                       
  // Across the event                                                                                                                                   
  //--------------------                                                                                                                                     
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  //-----------------                                                                                                                                          
  // Root Variables                                                                                                                                          
  //-----------------                                                                                                                                           
  TFile* fout_;
  TTree* tree_;

  unsigned int run, event, lumiblock, nprivtx;
  vector<string> *triggernames;
  vector<int> *triggerprescales;

  //----------                                                                                                                                              
  // dimuon                                                                                                                                                 
  //----------                                                                                                                                             
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz;
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz;
  vector<double> *mumutrkr, *mumutrkz , *mumudca;
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;
  vector<double> *mumucosalphabs, *mumucosalphabserr;
  vector<double> *mumumass, *mumumasserr;

  //-----------------------                                                                                                                                 
  // soft muon variables                                                                                                                                     
  //-----------------------                                                                                                                                  
  vector<bool>   *mumisgoodmuon, *mupisgoodmuon ;
  vector<int>    *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers;
  vector<int>    *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers;
  vector<double> *mumnormchi2, *mupnormchi2;
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx;
  vector<string> *mumtriglastfilter, *muptriglastfilter;
  vector<double> *mumpt, *muppt, *mumeta, *mupeta;

  //--------------
  // kaon track   
  //--------------                                                                                                                                                          
  vector<int> *trkchg; // +1 for K+, -1 for K-                                                                                                                            
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *trkdcabs, *trkdcabserr;

  //--------------
  // phi(1020)
  //--------------
  vector<int>    *kpchg;
  vector<double> *kppx, *kppy, *kppz ;
  vector<int>    *kmchg;
  vector<double> *kmpx, *kmpy, *kmpz ;
  vector<double> *phipx, *phipy, *phipz;
  vector<double> *phivtxx, *phivtxy, *phivtxz, *phivtxcl, *philsbs, *philsbserr;

  vector<double> *phimass, *phimasserr, *phibarmass, *phibarmasserr;

  //-----------------
  // Bs and Bsbar
  //-----------------
  int nb;
  vector<double> *bpx, *bpxerr, *bpy, *bpyerr, *bpz, *bpzerr ;
  vector<double> *bmass, *bmasserr;
  vector<double> *bvtxcl, *bvtxx, *bvtxxerr, *bvtxy, *bvtxyerr, *bvtxz, *bvtxzerr;
  vector<double> *bcosalphabs, *bcosalphabserr, *bcosalphabs2d, *bcosalphabs2derr, *blsbs, *blsbserr, *bctau, *bctauerr; 

  vector<double> *bbarmass, *bbarmasserr;

  //----------
  // For MC   
  //----------                                                                                                                                                              
  double genbpx, genbpy, genbpz;
  double genphipx, genphipy, genphipz;
  double genphivtxx, genphivtxy, genphivtxz;

  //                                                                                                                                                                        
  int genkpchg;
  double genkppx, genkppy, genkppz;
  int genkmchg;
  double genkmpx, genkmpy, genkmpz;

  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;
  
  string decname;

  vector<bool> *istruemum, *istruemup, *istruekp, *istruekm, *istruebs;

  //-----------------------
  // variables to monitor  
  //-----------------------                                                                                                                                                 
  TDatime t_begin_ , t_now_ ;
  int n_processed_, n_selected_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BsToPhiMuMu::BsToPhiMuMu(const edm::ParameterSet& iConfig):

  OutputFileName_(iConfig.getParameter<string>("OutputFileName")),
  BuildBsToPhiMuMu_(iConfig.getUntrackedParameter<bool>("BuildBsToPhiMuMu")),

  // particle properties                                                                                                                                                    
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  KaonMass_(iConfig.getUntrackedParameter<double>("KaonMass")),
  KaonMassErr_(iConfig.getUntrackedParameter<double>("KaonMassErr")),
  BsMass_(iConfig.getUntrackedParameter<double>("BsMass")),

  // labels                                                                                                                                                                 
  GenParticlesLabel_(iConfig.getParameter<edm::InputTag>("GenParticlesLabel")),
  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
  TrackLabel_(iConfig.getParameter<edm::InputTag>("TrackLabel")),
  TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
  LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),

  // gen particle                                                                                                                                                           
  IsMonteCarlo_(iConfig.getUntrackedParameter<bool>("IsMonteCarlo")),
  KeepGENOnly_(iConfig.getUntrackedParameter<bool>("KeepGENOnly")),
  TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),
  TruthMatchKaonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchKaonMaxR")),

  // pre-selection cuts                                                                                                                                                     
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMinDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMinDcaSigBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),

  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),

  PhiMinMass_(iConfig.getUntrackedParameter<double>("PhiMinMass")),
  PhiMaxMass_(iConfig.getUntrackedParameter<double>("PhiMaxMass")),
  BsMinVtxCl_(iConfig.getUntrackedParameter<double>("BsMinVtxCl")),
  BsMinMass_(iConfig.getUntrackedParameter<double>("BsMinMass")),
  BsMaxMass_(iConfig.getUntrackedParameter<double>("BsMaxMass")),

  tree_(0),
  triggernames(0), triggerprescales(0),
  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0),
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0), muppz(0),
  mumutrkr(0), mumutrkz(0), mumudca(0),  mumuvtxcl(0),  mumulsbs(0),
  mumulsbserr(0), mumucosalphabs(0),  mumucosalphabserr(0),
  mumumass(0), mumumasserr(0),
  mumisgoodmuon(0), mupisgoodmuon(0),
  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0),
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0),
  mumnormchi2(0), mupnormchi2(0), mumdxyvtx(0), mupdxyvtx(0),
  mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0),
  mumpt(0), muppt(0), mumeta(0), mupeta(0),

  trkchg(0), trkpx(0), trkpy(0), trkpz(0), trkpt(0),
  trkdcabs(0), trkdcabserr(0),

  kpchg(0),
  kppx(0), kppy(0), kppz(0),
  
  kmchg(0),
  kmpx(0), kmpy(0), kmpz(0),


  phipx(0), phipy(0), phipz(0),
  phivtxx(0), phivtxy(0), phivtxz(0),

  phimass(0), phimasserr(0), phibarmass(0), phibarmasserr(0),

  nb(0), bpx(0), bpxerr(0), bpy(0), bpyerr(0), bpz(0), bpzerr(0), bmass(0), bmasserr(0),
  bvtxcl(0), bvtxx(0), bvtxxerr(0), bvtxy(0), bvtxyerr(0), bvtxz(0), bvtxzerr(0),
  bcosalphabs(0), bcosalphabserr(0), bcosalphabs2d(0), bcosalphabs2derr(0), blsbs(0), blsbserr(0), bctau(0), bctauerr(0),

  bbarmass(0), bbarmasserr(0),

  genbpx(0), genbpy(0), genbpz(0),
  genphipx(0), genphipy(0), genphipz(0), genphivtxx(0), genphivtxy(0), genphivtxz(0),

  genkpchg(0),
  genkppx(0), genkppy(0), genkppz(0),

  genkmchg(0),
  genkmpx(0), genkmpy(0), genkmpz(0),

  genmumpx(0), genmumpy(0), genmumpz(0),
  genmuppx(0), genmuppy(0), genmuppz(0),

  decname(""),

  istruemum(0), istruemup(0), istruekp(0), istruekm(0), istruebs(0)


{
   //now do what ever initialization is needed
  assert(TriggerNames_.size() == LastFilterNames_.size());
  for (size_t i = 0; i < TriggerNames_.size(); ++i)
    mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];

}


BsToPhiMuMu::~BsToPhiMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BsToPhiMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  clearVariables();

  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock();

  n_processed_ += 1;
  histos[h_events]->Fill(0);

  if (IsMonteCarlo_) saveGenInfo(iEvent);

  hltReport(iEvent);

  if ( KeepGENOnly_){
    tree_->Fill();
    n_selected_ += 1;
  }else{
    if ( hasBeamSpot(iEvent) ) {
      iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
      if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) {
        buildBsToPhiMuMu(iEvent) ;
        if (IsMonteCarlo_) saveTruthMatch(iEvent);
        n_selected_ += 1;
      }
    }

    if (IsMonteCarlo_ || nb > 0){ // Keep failed events for MC to calculate reconstruction efficiency.                                                               
      tree_->Fill();
    }
  }


  clearVariables();


}


// ------------ method called once each job just before starting event loop  ------------
void 
BsToPhiMuMu::beginJob()
{

  t_begin_.Set();
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;


  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();

  for(int i=0; i<kHistNameSize; i++) {
    histos[i] = new TH1F(hist_args[i].name, hist_args[i].title,
                         hist_args[i].n_bins,
                         hist_args[i].x_min, hist_args[i].x_max);

  }

  tree_ = new TTree ("tree", "BsToPhiMuMu");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);
  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);
  tree_->Branch("trkchg", &trkchg);
  tree_->Branch("trkpx", &trkpx);
  tree_->Branch("trkpy", &trkpy);
  tree_->Branch("trkpz", &trkpz);
  tree_->Branch("trkpt", &trkpt);
  tree_->Branch("trkdcabs", &trkdcabs);
  tree_->Branch("trkdcabserr", &trkdcabserr);
  tree_->Branch("kpchg", &kpchg);
  tree_->Branch("kppx", &kppx);
  tree_->Branch("kppy", &kppy);
  tree_->Branch("kppz", &kppz);
  tree_->Branch("kmchg", &kmchg);
  tree_->Branch("kmpx", &kmpx);
  tree_->Branch("kmpy", &kmpy);
  tree_->Branch("kmpz", &kmpz);

  tree_->Branch("phipx", &phipx);
  tree_->Branch("phipy", &phipy);
  tree_->Branch("phipz", &phipz);
  tree_->Branch("phivtxx", &phivtxx);
  tree_->Branch("phivtxy", &phivtxy);
  tree_->Branch("phivtxz", &phivtxz);
  tree_->Branch("phimass", &phimass);
  tree_->Branch("phimasserr", &phimasserr);
  tree_->Branch("phibarmass", &phibarmass);
  tree_->Branch("phibarmasserr", &phibarmasserr);
  tree_->Branch("nb", &nb, "nb/I");
  tree_->Branch("bpx", &bpx);
  tree_->Branch("bpxerr", &bpxerr);
  tree_->Branch("bpy", &bpy);
  tree_->Branch("bpyerr", &bpyerr);
  tree_->Branch("bpz", &bpz);
  tree_->Branch("bpzerr", &bpzerr);
  tree_->Branch("bmass", &bmass);
  tree_->Branch("bmasserr", &bmasserr);
  tree_->Branch("bvtxcl", &bvtxcl);
  tree_->Branch("bvtxx", &bvtxx);
  tree_->Branch("bvtxxerr", &bvtxxerr);
  tree_->Branch("bvtxy", &bvtxy);
  tree_->Branch("bvtxyerr", &bvtxyerr);
  tree_->Branch("bvtxz", &bvtxz);
  tree_->Branch("bvtxzerr", &bvtxzerr);
  tree_->Branch("bcosalphabs", &bcosalphabs);
  tree_->Branch("bcosalphabserr", &bcosalphabserr);
  tree_->Branch("bcosalphabs2d", &bcosalphabs2d);
  tree_->Branch("bcosalphabs2derr", &bcosalphabs2derr);
  tree_->Branch("blsbs", &blsbs);
  tree_->Branch("blsbserr", &blsbserr);
  tree_->Branch("bctau", &bctau);
  tree_->Branch("bctauerr", &bctauerr);
  tree_->Branch("bbarmass", &bbarmass);
  tree_->Branch("bbarmasserr", &bbarmasserr);

  if (IsMonteCarlo_) {
    tree_->Branch("genbpx",      &genbpx     , "genbpx/D"    );
    tree_->Branch("genbpy",      &genbpy     , "genbpy/D"    );
    tree_->Branch("genbpz",      &genbpz     , "genbpz/D"    );
    tree_->Branch("genphipx",    &genphipx   , "genphipx/D"  );
    tree_->Branch("genphipy",    &genphipy   , "genphipy/D"  );
    tree_->Branch("genphipz",    &genphipz   , "genphipz/D"  );
    tree_->Branch("genphivtxx",  &genphivtxx    , "genphivtxx/D"   );
    tree_->Branch("genphivtxy",  &genphivtxy    , "genphivtxy/D"   );
    tree_->Branch("genphivtxz",  &genphivtxz    , "genphivtxz/D"   );
    tree_->Branch("genkpchg",   &genkpchg   , "genkpchg/I"   );
    tree_->Branch("genkppx",    &genkppx    , "genkppx/D"   );
    tree_->Branch("genkppy",    &genkppy    , "genkppy/D"   );
    tree_->Branch("genkppz",    &genkppz    , "genkppz/D"   );
    tree_->Branch("genkmchg",   &genkmchg   , "genkmchg/I"   );
    tree_->Branch("genkmpx",    &genkmpx    , "genkmpx/D"   );
    tree_->Branch("genkmpy",    &genkmpy    , "genkmpy/D"   );
    tree_->Branch("genkmpz",    &genkmpz    , "genkmpz/D"   );
    tree_->Branch("genmumpx",   &genmumpx   , "genmumpx/D"  );
    tree_->Branch("genmumpy",   &genmumpy   , "genmumpy/D"  );
    tree_->Branch("genmumpz",   &genmumpz   , "genmumpz/D"  );
    tree_->Branch("genmuppx",   &genmuppx   , "genmuppx/D"  );
    tree_->Branch("genmuppy",   &genmuppy   , "genmuppy/D"  );
    tree_->Branch("genmuppz",   &genmuppz   , "genmuppz/D"  );

    tree_->Branch("decname",    &decname);
    tree_->Branch("istruemum",  &istruemum );
    tree_->Branch("istruemup",  &istruemup );
    tree_->Branch("istruekp",   &istruekp  );
    tree_->Branch("istruekm",   &istruekm  );
    tree_->Branch("istruebs",   &istruebs  );

  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BsToPhiMuMu::endJob() 
{


  fout_->cd();
  tree_->Write();

  for(int i = 0; i < kHistNameSize; i++) {
    histos[i]->Write();
    histos[i]->Delete();
  }
  fout_->Close();

  t_now_.Set();
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();
  printf(" processed: %i \n selected: %i \n \
 duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_,
	 t_now_.Convert() - t_begin_.Convert(),
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

}

// ------------ method called when starting to processes a run  ------------
void 
BsToPhiMuMu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BsToPhiMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BsToPhiMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BsToPhiMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BsToPhiMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
BsToPhiMuMu::clearVariables(){

  run = 0;
  event = 0;
  lumiblock = 0;
  nprivtx = 0;
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();  mumpz->clear();
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear();  muppz->clear();
  mumutrkr->clear(); mumutrkz->clear();
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear();
  mumucosalphabs->clear();  mumucosalphabserr->clear();
  mumumass->clear(); mumumasserr->clear();
  mumisgoodmuon->clear();  mupisgoodmuon->clear();
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear();  mupnpixlayers->clear();
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();  mupntrklayers->clear();

  mumnormchi2->clear(); mupnormchi2->clear();
  mumdxyvtx->clear(); mupdxyvtx->clear();
  mumdzvtx->clear(); mupdzvtx->clear();
  mumtriglastfilter->clear(); muptriglastfilter->clear();
  mumpt->clear(); muppt->clear();
  mumeta->clear(); mupeta->clear();

  trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear(); trkpt->clear();
  trkdcabs->clear(); trkdcabserr->clear();

  kpchg->clear();
  kppx->clear(); kppy->clear(); kppz->clear();
 
  kmchg->clear();
  kmpx->clear(); kmpy->clear(); kmpz->clear();

  phipx->clear(); phipy->clear(); phipz->clear();
  phivtxx->clear(); phivtxy->clear(); phivtxz->clear();

  phimass->clear(); phimasserr->clear();
  phibarmass->clear(); phibarmasserr->clear();

  nb = 0;

  bpx->clear(); bpxerr->clear(); bpy->clear();  bpyerr->clear();
  bpz->clear(); bpzerr->clear();

  bmass->clear(); bmasserr->clear();
  bvtxcl->clear(); bvtxx->clear(); bvtxxerr->clear(); bvtxy->clear(); bvtxyerr->clear();
  bvtxz->clear(); bvtxzerr->clear(); bcosalphabs->clear(); bcosalphabserr->clear();
  bcosalphabs2d->clear(); bcosalphabs2derr->clear();
  blsbs->clear(); blsbserr->clear(); bctau->clear(); bctauerr->clear();

  bbarmass->clear(); bbarmasserr->clear();

  if (IsMonteCarlo_) {

    genbpx = 0;  genbpy = 0;  genbpz = 0;
    genphipx = 0;  genphipy = 0;  genphipz = 0;
    genphivtxx = 0; genphivtxy = 0; genphivtxz = 0;

    genkpchg = 0;
    genkppx = 0;  genkppy = 0;  genkppz = 0;
    genkmchg = 0;
    genkmpx = 0;  genkmpy = 0;  genkmpz = 0;

    genmumpx = 0;  genmumpy = 0;  genmumpz = 0;
    genmuppx = 0;  genmuppy = 0;  genmuppz = 0;

    decname = "";
    istruemum->clear(); istruemup->clear(); istruekp->clear();
    istruekm->clear(); istruebs->clear();


  }

}


void
BsToPhiMuMu::hltReport(const edm::Event& iEvent)
{

  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT")
      << __LINE__ << " : couldn't get handle on HLT Trigger" ; }

  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case.                                                                                                                          
      if ((*hltTriggerResults)[itrig].accept() == 1){

        string triggername = triggerNames_.triggerName(itrig);
        int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

        // Loop over our interested HLT trigger names to find if this event contains.                                                                               
        for (unsigned int it=0; it<TriggerNames_.size(); it++){
          if (triggername.find(TriggerNames_[it]) != string::npos) {
            // save the no versioned case                                                                                                                            
            triggernames->push_back(TriggerNames_[it]);
            triggerprescales->push_back(triggerprescale);

          }}}}}

}

bool
BsToPhiMuMu::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);

  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }

  beamSpot_ = *beamSpotHandle;
  return true;
}

bool
BsToPhiMuMu::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size();

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;

  return true;
}

//-------------------------------------------------------
//  main function to build BsToPhiMuMu candidate
//-------------------------------------------------------

bool
BsToPhiMuMu::buildBsToPhiMuMu(const edm::Event& iEvent)
{

  // init variables
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return false;

  edm::Handle< vector<pat::GenericParticle> >thePATTrackHandle;
  iEvent.getByLabel(TrackLabel_, thePATTrackHandle);

  bool passed;
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu;
  double trk_R, trk_Z, trk_DCA;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  double MuMuLSBS, MuMuLSBSErr;
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  double trk_pt; 
  reco::TransientTrack refitKmTT, refitKpTT;
  double phi_mass, phi_vtx_cl;
  double b_vtx_chisq, b_vtx_cl, b_mass;
  //double bbar_mass;
  double DCAPhiTrkBS, DCAPhiTrkBSErr;
  RefCountedKinematicTree vertexFitTree, barVertexFitTree;


  // ---------------------------------
  // loop 1: mu-
  // ---------------------------------
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin();
       iMuonM != patMuonHandle->end(); iMuonM++){

    reco::TrackRef muTrackm = iMuonM->innerTrack();
    if ( muTrackm.isNull() ) continue;

    histos[h_mupt]->Fill(muTrackm->pt());
    histos[h_mueta]->Fill(muTrackm->eta());

    if ( (muTrackm->charge() != -1) ||
	 (muTrackm->pt() < MuonMinPt_) ||
	 (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;

    // check mu- DCA to beam spot
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));
    passed = hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ;
    histos[h_mumdcabs]->Fill(DCAmumBS);
    if ( ! passed ) continue;


    // ---------------------------------
    // loop 2: mu+
    // ---------------------------------
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin();
	 iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack();
      if ( muTrackp.isNull() ||
	   (muTrackp->charge() != 1) ||
	   (muTrackp->pt() < MuonMinPt_) ||
	   (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;

      // check mu+ DCA to beam spot
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));
      passed = hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr);
      if ( ! passed ) continue;

      // check goodness of muons closest approach and the 3D-DCA
      // passed = hasGoodClosestApproachTracks(muTrackpTT, muTrackmTT,
      //     mumutrk_R, mumutrk_Z, DCAmumu);
      if ( !calClosestApproachTracks(muTrackpTT, muTrackmTT,
				     mumutrk_R, mumutrk_Z, DCAmumu)) continue;
      histos[h_mumutrkr]->Fill(mumutrk_R);
      histos[h_mumutrkz]->Fill(mumutrk_Z);
      histos[h_mumudca]->Fill(DCAmumu);
      // if ( !passed ) continue;
      if ( mumutrk_R > TrkMaxR_ ||
	      mumutrk_Z > TrkMaxZ_ ||
	   DCAmumu > MuMuMaxDca_ ) continue;

      // check dimuon vertex
      passed = hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT,
				 mu_mu_vtx_cl, mu_mu_pt,
				 mu_mu_mass, mu_mu_mass_err,
				 MuMuLSBS, MuMuLSBSErr,
				 MuMuCosAlphaBS, MuMuCosAlphaBSErr);

      ///cout << "hasGoodMuMuVertex: " << boolalpha << passed << endl;

      histos[h_mumuvtxcl]->Fill(mu_mu_vtx_cl);
      histos[h_mumupt]->Fill(mu_mu_pt);
      histos[h_mumumass]->Fill(mu_mu_mass);
      histos[h_mumulxybs]->Fill(MuMuLSBS/MuMuLSBSErr);
      histos[h_mumucosalphabs]->Fill(MuMuCosAlphaBS);
      if ( !passed) continue;


      // ---------------------------------
      // loop 3: track-
      // ---------------------------------
      for ( vector<pat::GenericParticle>::const_iterator iTrackM
	      = thePATTrackHandle->begin();
	    iTrackM != thePATTrackHandle->end(); ++iTrackM ) {
	
	reco::TrackRef Trackm = iTrackM->track();
	if ( Trackm.isNull() || (Trackm->charge() != -1) ) continue;

	passed = hasGoodTrack(iEvent, *iTrackM, trk_pt);
	histos[h_trkpt]->Fill(trk_pt);
	if (!passed) continue;
	
	// compute track DCA to beam spot
	const reco::TransientTrack theTrackmTT(Trackm, &(*bFieldHandle_));
	passed = hasGoodTrackDcaBs(theTrackmTT, DCAPhiTrkBS, DCAPhiTrkBSErr);
	histos[h_trkdcasigbs]->Fill(DCAPhiTrkBS/DCAPhiTrkBSErr);
	if (!passed) continue;

	// ---------------------------------
	// loop 4: track+
	// ---------------------------------
	for ( vector<pat::GenericParticle>::const_iterator iTrackP
		= thePATTrackHandle->begin();
	      iTrackP != thePATTrackHandle->end(); ++iTrackP ) {
	  
	  reco::TrackRef Trackp = iTrackP->track();
	  if ( Trackp.isNull() || (Trackp->charge() != 1) ) continue;

	  passed = hasGoodTrack(iEvent, *iTrackP, trk_pt);
	  // histos[h_trkpt]->Fill(trk_pt);
	  if (!passed) continue;
	  
	  // compute track DCA to beam spot
	  const reco::TransientTrack theTrackpTT(Trackp, &(*bFieldHandle_));
	  passed = hasGoodTrackDcaBs(theTrackpTT, DCAPhiTrkBS, DCAPhiTrkBSErr);
	  if (!passed) continue;


	  // check goodness of two tracks closest approach and the 3D-DCA
	  if (! calClosestApproachTracks(theTrackpTT, theTrackmTT,
					 trk_R, trk_Z, trk_DCA)) continue ;
	  if ( trk_R > TrkMaxR_ || trk_Z > TrkMaxZ_ ) continue;
	  


	  // check two tracks vertex for Bs
	  if ( ! hasGoodPhiVertex(theTrackmTT, theTrackpTT, refitKmTT, refitKpTT, phi_vtx_cl, phi_mass) ) continue;
	  if ( phi_mass < PhiMinMass_ || phi_mass > PhiMaxMass_ ) continue;

	  // check two tracks vertex for Bsbar
	  if ( ! hasGoodPhiVertex(theTrackpTT, theTrackmTT, refitKpTT, refitKmTT, phi_vtx_cl, phi_mass) ) continue;
	  if ( phi_mass < PhiMinMass_ || phi_mass > PhiMaxMass_ ) continue;




	  // fit Bs vertex  mu- mu+ K- K+
	  if ( ! hasGoodBsVertex(muTrackmTT, muTrackpTT, theTrackmTT, theTrackpTT,
				 b_vtx_chisq, b_vtx_cl, b_mass,
				 vertexFitTree) ) continue;

	  if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) continue;

	  // fit Bsbar vertex mu- mu+ K+ K-
	  if ( ! hasGoodBsVertex(muTrackmTT, muTrackpTT, theTrackpTT, theTrackmTT,
				 b_vtx_chisq, b_vtx_cl, b_mass,
				 barVertexFitTree) ) continue;

	  if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) continue;



	  // need to check with primaryVertex tracks?

	  nb++;

	  
	  // save the tree variables
	  saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
			    mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl,
			    MuMuLSBS, MuMuLSBSErr,
			    MuMuCosAlphaBS, MuMuCosAlphaBSErr,
			    mu_mu_mass, mu_mu_mass_err);
	  

	  saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp);
	  // saveBsToPhiMuMu(vertexFitTree);

	  phimass->push_back(phi_mass);
	  bmass->push_back(b_mass);

	} // close track+ loop
      } // close track- loop
    } // close mu+ loop
  } // close mu- loop

  if ( nb > 0) {
    edm::LogInfo("myBs") << "Found " << nb << " Bs -> phi(KK) mu+ mu- ";
    printf("------------------------------------------------------------\n");
    return true;
  }
  return false;


}

void
BsToPhiMuMu::calLS (double Vx, double Vy, double Vz,
		     double Wx, double Wy, double Wz,
		     double VxErr2, double VyErr2, double VzErr2,
		     double VxyCov, double VxzCov, double VyzCov,
		     double WxErr2, double WyErr2, double WzErr2,
		     double WxyCov, double WxzCov, double WyzCov,
		     double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


void
BsToPhiMuMu::calCosAlpha (double Vx, double Vy, double Vz,
			   double Wx, double Wy, double Wz,
			   double VxErr2, double VyErr2, double VzErr2,
			   double VxyCov, double VxzCov, double VyzCov,
			   double WxErr2, double WyErr2, double WzErr2,
			   double WxyCov, double WxzCov, double WyzCov,
			   double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			 (Wnorm*Wnorm*Wnorm*Wnorm) +
			 
			 ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			  
			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			 (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }

}


void
BsToPhiMuMu::calCosAlpha2d (double Vx, double Vy, double Vz,
			     double Wx, double Wy, double Wz,
			     double VxErr2, double VyErr2, double VzErr2,
			     double VxyCov, double VxzCov, double VyzCov,
			     double WxErr2, double WyErr2, double WzErr2,
			     double WxyCov, double WxzCov, double WyzCov,
			     double* cosAlpha2d, double* cosAlpha2dErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha2d = VdotW / (Vnorm * Wnorm);
    *cosAlpha2dErr = sqrt( (
			    (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +

			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +

			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha2d = 0.;
    *cosAlpha2dErr = 0.;
  }

}


bool
BsToPhiMuMu::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT,
				double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    muTrackTT.trajectoryStateClosestToPoint(
					    GlobalPoint(beamSpot_.position().x(),
							beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(muDcaBs) > MuonMaxDcaBs_ )   return false;
  return true;
}

bool
BsToPhiMuMu::hasGoodTrackDcaBs (const reco::TransientTrack TrackTT,
				 double &DcaBs, double &DcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    TrackTT.trajectoryStateClosestToPoint(
					  GlobalPoint(beamSpot_.position().x(),
						      beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  DcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(DcaBs/DcaBsErr) < TrkMinDcaSigBs_ )   return false;
  return true;
}


bool
BsToPhiMuMu::hasGoodTrackDcaPoint (const reco::TransientTrack track,
				    const GlobalPoint p,
				    double maxdca, double &dca, double &dcaerr)
{
  TrajectoryStateClosestToPoint theDCAX = track.trajectoryStateClosestToPoint(p);
  if ( !theDCAX.isValid() ) return false;

  dca = theDCAX.perigeeParameters().transverseImpactParameter();
  dcaerr = theDCAX.perigeeError().transverseImpactParameterError();
  if ( dca > maxdca ) return false;

  return true;
}


bool
BsToPhiMuMu::calClosestApproachTracks (const reco::TransientTrack trackpTT,
					const reco::TransientTrack trackmTT,
					double & trk_R,
					double & trk_Z,
					double & trk_DCA)
{
  ClosestApproachInRPhi ClosestApp;
  ClosestApp.calculate(trackpTT.initialFreeState(),
		       trackmTT.initialFreeState());
  if (! ClosestApp.status() )  return false ;

  GlobalPoint XingPoint = ClosestApp.crossingPoint();

  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
  trk_Z = fabs(XingPoint.z());

  // if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) >
  //      TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_))  return false;

  trk_DCA = ClosestApp.distance();
  // if (DCAmumu > MuMuMaxDca_) return false;

  return true;
}


bool
BsToPhiMuMu::hasGoodMuMuVertex ( const reco::TransientTrack muTrackpTT,
				 const reco::TransientTrack muTrackmTT,
				 reco::TransientTrack &refitMupTT,
				 reco::TransientTrack &refitMumTT,
				 double & mu_mu_vtx_cl, double & mu_mu_pt,
				 double & mu_mu_mass, double & mu_mu_mass_err,
				 double & MuMuLSBS, double & MuMuLSBSErr,
				 double & MuMuCosAlphaBS,
				 double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,
					       MuonMass_,chi,ndf,MuonMassErr_));

  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);

  if ( !mumuVertexFitTree->isValid())  return false;

  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();

  if ( !mumu_KV->vertexIsValid()) return false;

  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
			     int(rint(mumu_KV->degreesOfFreedom())));

  if (mu_mu_vtx_cl < MuMuMinVtxCl_)  return false;

  // extract the re-fitted tracks
  mumuVertexFitTree->movePointerToTheTop();

  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();

  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();

  TLorentzVector mymum, mymup, mydimu;

  mymum.SetXYZM(refitMumTT.track().momentum().x(),
		refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_);

  mymup.SetXYZM(refitMupTT.track().momentum().x(),
		refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_);

  mydimu = mymum + mymup;
  mu_mu_pt = mydimu.Perp();

  mu_mu_mass = mumu_KP->currentState().mass();
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			matrix()(6,6));

  if ((mu_mu_pt < MuMuMinPt_) || (mu_mu_mass < MuMuMinInvMass_) ||
      (mu_mu_mass > MuMuMaxInvMass_))  return false;

  // compute the distance between mumu vtx and beam spot
  calLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	 mumu_KV->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &MuMuLSBS,&MuMuLSBSErr);

  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  calCosAlpha(mumu_KP->currentState().globalMomentum().x(),
	      mumu_KP->currentState().globalMomentum().y(),
	      0.0,
	      mumu_KV->position().x() - beamSpot_.position().x(),
	      mumu_KV->position().y() - beamSpot_.position().y(),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
	      mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
	      0.0,
	      0.0,
	      mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
	      mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
	      0.0,
	      mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
	      0.0,
	      0.0,
	      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);

  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;

  return true;
}


bool
BsToPhiMuMu::matchMuonTrack (const edm::Event& iEvent,
			      const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel(MuonLabel_, thePATMuonHandle);

  reco::TrackRef muTrackRef;
  for (vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
       iMuon != thePATMuonHandle->end(); iMuon++){

    muTrackRef = iMuon->innerTrack();
    if ( muTrackRef.isNull() ) continue;

    if (muTrackRef == theTrackRef) return true;
  }

  return false;
}


bool
BsToPhiMuMu::matchMuonTracks (const edm::Event& iEvent,
			       const vector<reco::TrackRef> theTracks)
{
  reco::TrackRef theTrackRef;
  for(unsigned int j = 0; j < theTracks.size(); ++j) {
    theTrackRef = theTracks[j];
    if ( matchMuonTrack(iEvent, theTrackRef) ) return true;
  }
  return false;
}


bool
BsToPhiMuMu::hasGoodTrack(const edm::Event& iEvent,
			  const pat::GenericParticle iTrack,
			  double & kaon_trk_pt)
{
  reco::TrackRef theTrackRef = iTrack.track();
  if ( theTrackRef.isNull() ) return false;

  // veto muon tracks
  if ( matchMuonTrack(iEvent, theTrackRef) ) return false;

  // veto pion tracks from Kshort
  //if ( matchKshortTrack(iEvent, theTrackRef) ) return false;

  // check the track kinematics
  kaon_trk_pt = theTrackRef->pt();

  if ( theTrackRef->pt() < TrkMinPt_ ) return false;

  return true;
}


bool
BsToPhiMuMu::hasGoodPhiVertex( const reco::TransientTrack kaonmTT,
			       const reco::TransientTrack kaonpTT,
			       reco::TransientTrack &refitKmTT,
			       reco::TransientTrack &refitKpTT,
			       double & phi_vtx_cl, double & phi_mass)
{
  KinematicParticleFactoryFromTransientTrack pFactory;

  float chi = 0.;
  float ndf = 0.;

  vector<RefCountedKinematicParticle> phiParticles;
  phiParticles.push_back(pFactory.particle(kaonmTT,KaonMass_,chi,ndf,KaonMassErr_));
  phiParticles.push_back(pFactory.particle(kaonpTT,KaonMass_,chi,ndf,KaonMassErr_));

  KinematicParticleVertexFitter fitter;
  RefCountedKinematicTree phiVertexFitTree = fitter.fit(phiParticles);
  if ( ! phiVertexFitTree->isValid() ) return false ;

  phiVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle phi_KP = phiVertexFitTree->currentParticle();
  RefCountedKinematicVertex phi_KV   = phiVertexFitTree->currentDecayVertex();
  if ( !phi_KV->vertexIsValid() ) return false;

  phi_vtx_cl = TMath::Prob((double)phi_KV->chiSquared(),
			   int(rint(phi_KV->degreesOfFreedom())));


  // extract the re-fitted tracks
  phiVertexFitTree->movePointerToTheTop();

  phiVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitKm = phiVertexFitTree->currentParticle();
  refitKmTT = refitKm->refittedTransientTrack();

  phiVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitKp = phiVertexFitTree->currentParticle();
  refitKpTT = refitKp->refittedTransientTrack();

  TLorentzVector mykm, mykp, myphi;

  mykm.SetXYZM(refitKmTT.track().momentum().x(),
	       refitKmTT.track().momentum().y(),
	       refitKmTT.track().momentum().z(), KaonMass_);

  mykp.SetXYZM(refitKpTT.track().momentum().x(),
	       refitKpTT.track().momentum().y(),
	       refitKpTT.track().momentum().z(), KaonMass_);

  myphi = mykm + mykp;
  //phi_pt = myphi.Perp();

  phi_mass = phi_KP->currentState().mass();
  //phi_mass_err = sqrt(phi_KP->currentState().kinematicParametersError().matrix()(6,6));

  //if ( (phi_pt < PhiMinPt_) || (phi_mass < PhiMinMass_) || (phi_mass > PhiMaxMass_) )  return false;

  return true;

}

bool
BsToPhiMuMu::hasGoodBsVertex(const reco::TransientTrack mu1TT,
			     const reco::TransientTrack mu2TT,
			     const reco::TransientTrack kaonmTT,
			     const reco::TransientTrack kaonpTT,
			     double & b_vtx_chisq, double & b_vtx_cl,
			     double & b_mass,
			     RefCountedKinematicTree & vertexFitTree)
{

  KinematicParticleFactoryFromTransientTrack pFactory;
  float chi = 0.;
  float ndf = 0.;

  // Bs -> mu+ mu- phi (K+ K-)
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(kaonmTT, KaonMass_, chi,
					      ndf, KaonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(kaonpTT, KaonMass_, chi,
					      ndf, KaonMassErr_));


  KinematicParticleVertexFitter fitter;
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false;
  ////cout << "particles fitted to a single vertex found: " << boolalpha << vertexFitTree->isValid() << endl; 

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();

  if ( !b_KV->vertexIsValid()) return false;
  ////cout << "Bs decay vertex found: " << boolalpha << b_KV->vertexIsValid() << endl;

  b_vtx_cl = TMath::Prob((double)b_KV->chiSquared(),
			 int(rint(b_KV->degreesOfFreedom())));

  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  b_mass = b_KP->currentState().mass();

  //if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) return false;

  //printf("reco Bs cand vtxcl: %6.4f , mass: %6.4f \n", b_vtx_cl, b_mass); 

  return true;

}


void
BsToPhiMuMu::saveGenInfo(const edm::Event& iEvent){
  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel(GenParticlesLabel_, genparticles );

  // loop over all gen particles                                                                                                                      
  for( size_t i = 0; i < genparticles->size(); ++i ) {
    const reco::GenParticle & b = (*genparticles)[i];

    // only select Bs candidate                                                                                                                 
    if ( abs(b.pdgId()) != BS_PDG_ID ) continue;

    int imum(-1), imup(-1), iphi(-1), ikp(-1), ikm(-1);
    int ijpsi(-1), ipsi2s(-1);

    // loop over all Bs daughters                                                                                                       
    for ( size_t j = 0; j < b.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(b.daughter(j));

      if (dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
      if (dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      if (abs(dau.pdgId()) == PHI_PDG_ID) iphi = j;
      if (dau.pdgId() == JPSI_PDG_ID ) ijpsi = j;
      if (dau.pdgId() == PSI2S_PDG_ID ) ipsi2s = j;

    }

    if ( iphi == -1 ) continue;

    const reco::Candidate & phi = *(b.daughter(iphi));

    for ( size_t j = 0; j < phi.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(phi.daughter(j));
      if (dau.pdgId() == KAONPLUS_PDG_ID) ikp = j;
      if (dau.pdgId() == -KAONPLUS_PDG_ID) ikm = j;

    }

    if (ikp == -1 || ikm == -1) continue;


    // store the Bs and phi vars                                                                                                                        
    const reco::Candidate & kp = *(phi.daughter(ikp));
    const reco::Candidate & km = *(phi.daughter(ikm));

    const reco::Candidate *mum = NULL;
    const reco::Candidate *mup = NULL;

    //--------------------
    // Bs -> phi mu mu    
    //--------------------                                                                                                                       
    if (imum != -1 && imup != -1) {
      // cout << "Found GEN Bs -> phi mu mu " << endl;                                                                                                    
      mum = b.daughter(imum);
      mup = b.daughter(imup);
      decname = "BsToPhiMuMu";
    }

    //----------------------------
    // Bs -> phi J/psi(->mu mu)   
    //----------------------------                                                                                                            
    else if ( ijpsi != -1 ) {
      // cout << "Found GEN Bs --> phi J/psi " << endl;                                                                                                  
      const reco::Candidate & jpsi = *(b.daughter(ijpsi));
      for ( size_t j = 0; j < jpsi.numberOfDaughters(); ++j){
        const reco::Candidate  &dau = *(jpsi.daughter(j));
        if ( dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
        if ( dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      }
      if (imum != -1 && imup != -1) {
        mum = jpsi.daughter(imum);
        mup = jpsi.daughter(imup);
        decname = "BsToPhiJPsi";
      }
    }

    //------------------------------
    // Bs -> phi psi(2S)(->mu mu)   
    //------------------------------                                                                                                                    
    else if ( ipsi2s != -1) {
      // cout << "Found GEN Bs --> phi psi(2S) " << endl;                                                                                                 
      const reco::Candidate & psi2s = *(b.daughter(ipsi2s));
      for ( size_t j = 0; j < psi2s.numberOfDaughters(); ++j){
        const reco::Candidate  &dau = *(psi2s.daughter(j));
        if ( dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
        if ( dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      }
      if (imum != -1 && imup != -1) {
        mum = psi2s.daughter(imum);
        mup = psi2s.daughter(imup);
        decname = "BsToPhiPsi2S";
      }
    }

    if ( mum == NULL || mup == NULL) continue;

    // save gen info                                                                                                                                     
    genbpx = b.px();
    genbpy = b.py();
    genbpz = b.pz();

    genphipx = phi.px();
    genphipy = phi.py();
    genphipz = phi.pz();

    genphivtxx = phi.vx();
    genphivtxy = phi.vy();
    genphivtxz = phi.vz();

    genkpchg = kp.charge();
    genkppx = kp.px();
    genkppy = kp.py();
    genkppz = kp.pz();

    genkmchg = km.charge();
    genkmpx = km.px();
    genkmpy = km.py();
    genkmpz = km.pz();

    genmumpx = mum->px();
    genmumpy = mum->py();
    genmumpz = mum->pz();

    genmuppx = mup->px();
    genmuppy = mup->py();
    genmuppz = mup->pz();

  }

}



void
BsToPhiMuMu::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP,
				   reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{

  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight));
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight));
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());

  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());

  mumnormchi2->push_back(muTrackm->normalizedChi2());
  mupnormchi2->push_back(muTrackp->normalizedChi2());

  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position()));
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position()));

  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position()));
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position()));

  mumpt->push_back(muTrackm->pt());
  muppt->push_back(muTrackp->pt());

  mumeta->push_back(muTrackm->eta());
  mupeta->push_back(muTrackp->eta());

}


void
BsToPhiMuMu::saveDimuVariables(double DCAmumBS, double DCAmumBSErr,
			       double DCAmupBS, double DCAmupBSErr,
			       double mumutrk_R, double mumutrk_Z,
			       double DCAmumu,  double mu_mu_vtx_cl,
			       double MuMuLSBS, double MuMuLSBSErr,
			       double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
			       double mu_mu_mass, double mu_mu_mass_err)

{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z);
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl);
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr);

  mumumass->push_back(mu_mu_mass);
  mumumasserr->push_back(mu_mu_mass_err);
}


void
BsToPhiMuMu::saveMuonTriggerMatches(const pat::Muon iMuonM, const pat::Muon iMuonP)
{
  string mum_matched_lastfilter_name = "";
  string mup_matched_lastfilter_name = "";

  for(vector<string>::iterator it = triggernames->begin();
      it != triggernames->end(); ++it) {

    string hltLastFilterName = mapTriggerToLastFilter_[*it] ;

    const pat::TriggerObjectStandAloneCollection mumHLTMatches
      = iMuonM.triggerObjectMatchesByFilter( hltLastFilterName );
    const pat::TriggerObjectStandAloneCollection mupHLTMatches
      = iMuonP.triggerObjectMatchesByFilter( hltLastFilterName );

    if ( mumHLTMatches.size() > 0 )
      mum_matched_lastfilter_name.append(hltLastFilterName+" ") ;

    if ( mupHLTMatches.size() > 0 )
      mup_matched_lastfilter_name.append(hltLastFilterName+" ") ;
  }

  mumtriglastfilter->push_back(mum_matched_lastfilter_name);
  muptriglastfilter->push_back(mup_matched_lastfilter_name);
}



void
BsToPhiMuMu::saveTruthMatch(const edm::Event& iEvent){

}

//define this as a plug-in
DEFINE_FWK_MODULE(BsToPhiMuMu);
