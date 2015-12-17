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
// Original Author:  Niladribihari Sahoo,42 3-024,+41227662373,
//       copyright @ N.Sahoo, NISER, Bhubaneswar
//         Created:  Sat Nov 28 07:33:44 CET 2015
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
  bool BuildLbToLzMuMu_;

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

  //--------------------                                                                                                                                   
  // kaon track                                                                                                                                    
  //--------------------                                                                                                                                    
  //vector<int> *prchg;
  //vector<double> *prpx, *prpy, *prpz;
  //vector<int> *pichg;
  //vector<double> *pipx, *pipy, *pipz;

  //--------------
  // kaon track   
  //--------------                                                                                                                                                          
  vector<int> *trkchg; // +1 for K+, -1 for K-                                                                                                                            
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *trkdcabs, *trkdcabserr;

  // phi(1020)
  vector<double> *phimass;

  int nb;

  vector<double> *bpx, *bpxerr, *bpy, *bpyerr, *bpz, *bpzerr ;
  vector<double>  *bmass, *bmasserr;
  vector<double> *bvtxcl, *bvtxx, *bvtxxerr, *bvtxy, *bvtxyerr, *bvtxz, *bvtxzerr;
  vector<double> *bcosalphabs, *bcosalphabserr, *bcosalphabs2d, *bcosalphabs2derr, *blsbs, *blsbserr, *bctau, *bctauerr; 

  // Bs and Bsbar                                                                                                                                                           
  vector<double> *bbarmass, *bbarmasserr;

  // For MC                                                                                                                                                                 
  double genbpx, genbpy, genbpz;
  double genphipx, genphipy, genphipz;
  double genphivtxx, genphivtxy, genphivtxz;

  //                                                                                                                                                                        
  double genkppx, genkppy, genkppz;
  double genkmpx, genkmpy, genkmpz;

  //int gentrkchg;
  //double gentrkpx, gentrkpy, gentrkpz;
  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;
  
  string decname;

  vector<bool> *istruemum, *istruemup, *istruekp, *istruekm, *istruebs;


  // variables to monitor                                                                                                                                                   
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
BsToPhiMuMu::BsToPhiMuMu(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

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
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
BsToPhiMuMu::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BsToPhiMuMu::endJob() 
{
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

//define this as a plug-in
DEFINE_FWK_MODULE(BsToPhiMuMu);
