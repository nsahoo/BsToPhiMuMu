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
//====================================================================
// original author:  Niladribihari Sahoo,42 3-024,+41227662373,
//        copyright  @ N.Sahoo, NISER, Bhubaneswar
//         created:  Sat Nov 28 07:33:44 CET 2015
//====================================================================
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

  bool hasGoodPhiVertex(const vector<reco::TrackRef>,
		       RefCountedKinematicTree &);

  bool hasGoodPhiVertexMKC(const vector<reco::TrackRef>,
			  RefCountedKinematicTree &);


  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint,
                             double, double &, double &);

  bool hasGoodBsMass(RefCountedKinematicTree, double &);  

  bool hasGoodBsVertex(const reco::TransientTrack, const reco::TransientTrack,   
                       const vector<reco::TrackRef>, double &, double &, double &,      
                       RefCountedKinematicTree &, RefCountedKinematicTree &);

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
  //ParticleMass ProtonMass_;
  //float ProtonMassErr_;
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
  //double TruthMatchProtonMaxR_;

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
    tree_->Branch("genbpx",        &genbpx     , "genbpx/D"    );
    tree_->Branch("genbpy",        &genbpy     , "genbpy/D"    );
    tree_->Branch("genbpz",        &genbpz     , "genbpz/D"    );
    tree_->Branch("genphipx",      &genphipx   , "genphipx/D"  );
    tree_->Branch("genphipy",      &genphipy   , "genphipy/D"  );
    tree_->Branch("genphipz",      &genphipz   , "genphipz/D"  );
    tree_->Branch("genphivtxx",    &genphivtxx    , "genphivtxx/D"   );
    tree_->Branch("genphivtxy",    &genphivtxy    , "genphivtxy/D"   );
    tree_->Branch("genphivtxz",    &genphivtxz    , "genphivtxz/D"   );
    tree_->Branch("genkpchg",    &genkpchg   , "genkpchg/I"   );
    tree_->Branch("genkppx",     &genkppx    , "genkppx/D"   );
    tree_->Branch("genkppy",     &genkppy    , "genkppy/D"   );
    tree_->Branch("genkppz",     &genkppz    , "genkppz/D"   );
    tree_->Branch("genkmchg",    &genkmchg   , "genkmchg/I"   );
    tree_->Branch("genkmpx",     &genkmpx    , "genkmpx/D"   );
    tree_->Branch("genkmpy",     &genkmpy    , "genkmpy/D"   );
    tree_->Branch("genkmpz",     &genkmpz    , "genkmpz/D"   );
    tree_->Branch("genmumpx",    &genmumpx   , "genmumpx/D"  );
    tree_->Branch("genmumpy",    &genmumpy   , "genmumpy/D"  );
    tree_->Branch("genmumpz",    &genmumpz   , "genmumpz/D"  );
    tree_->Branch("genmuppx",    &genmuppx   , "genmuppx/D"  );
    tree_->Branch("genmuppy",    &genmuppy   , "genmuppy/D"  );
    tree_->Branch("genmuppz",    &genmuppz   , "genmuppz/D"  );

    tree_->Branch("decname",  &decname);
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



//define this as a plug-in
DEFINE_FWK_MODULE(BsToPhiMuMu);
