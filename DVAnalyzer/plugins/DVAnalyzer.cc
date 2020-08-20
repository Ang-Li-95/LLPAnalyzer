// -*- C++ -*-
//
// Package:    LLPAnalyzer/DVAnalyzer
// Class:      DVAnalyzer
//
/**\class DVAnalyzer DVAnalyzer.cc LLPAnalyzer/DVAnalyzer/plugins/DVAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ang Li
//         Created:  Thu, 11 Jun 2020 14:08:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "LLPAnalyzer/Formats/interface/TrackRefMap.h"
#include "LLPAnalyzer/DVAnalyzer/interface/TrackRescaler.h"
#include "LLPAnalyzer/Formats/interface/TrackAssociation.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;
using reco::TransientTrack;

typedef std::set<reco::TrackRef> track_set;
typedef std::vector<reco::TrackRef> track_vec;


struct eventInfo
{
  int evt;
  std::vector <int> vtx_track_size;
  std::vector <double> vtx_dBV;
  std::vector <double> vtx_sigma_dBV;
  std::vector <double> vtx_x;
  std::vector <double> vtx_y;
  std::vector <double> vtx_z;
};

class DVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
   public:
      explicit DVAnalyzer(const edm::ParameterSet&);
      ~DVAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void beginRun(edm::Run const& , edm::EventSetup const&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const& , edm::EventSetup const&) override;
      virtual void endJob() override;
      void initEventStructure();
      std::vector<TransientVertex> refitVtx(std::vector<reco::TransientTrack>& ttks);
      reco::TrackRefVector SharedTracks(track_set t1, track_set t2);
      std::vector<TransientTrack> combineTracks(track_set t1, track_set t2, std::map<reco::TrackRef, size_t> seed_track_ref_map, std::vector<TransientTrack> tTracks);
      track_set vertexTracks(const reco::Vertex v) const;
      track_vec vertex_track_vec(const reco::Vertex& v) const;
      bool is_track_subset(const track_set& a, const track_set& b) const;
      void print_vtx_track(const reco::Vertex v) const;
      track_set getJetTrack(const pat::Jet& jet, DVAna::UnpackedCandidateTracksMap& tRefMap, DVAna::RescaledTrackMap& tRescaledMap);
      std::vector<int> jet_vtx_nTracks(const pat::Jet& jet, double track_pt_min, DVAna::TrackAssociation vtx_tk_map, DVAna::UnpackedCandidateTracksMap& tRefMap, DVAna::RescaledTrackMap& tRescaledMap);
      double mag(double x, double y);


      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::BeamSpot> bsToken_;
      edm::Handle<reco::BeamSpot> bsHandle_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
      edm::Handle<std::vector<pat::Jet>> jetHandle_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
      edm::Handle<std::vector<pat::Muon>> muonHandle_;
      //edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tracksToken_;
      //edm::Handle<std::vector<pat::PackedCandidate>> tracksHandle_;
      edm::EDGetTokenT<reco::TrackRefVector> trackRefToken_;
      edm::Handle<reco::TrackRefVector> trackRefHandle_;
      edm::EDGetTokenT<DVAna::UnpackedCandidateTracksMap> trackRefMapToken_;
      edm::Handle<DVAna::UnpackedCandidateTracksMap> trackRefMapHandle_;
      edm::EDGetTokenT<DVAna::RescaledTrackMap> trackRescaledToken_;
      edm::Handle<DVAna::RescaledTrackMap> trackRescaledHandle_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
      edm::Handle<std::vector<reco::Vertex>> vtxHandle_;
      edm::EDGetTokenT<DVAna::TrackAssociation> trackAssoToken_;
      edm::Handle<DVAna::TrackAssociation> trackAssoHandle_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::Handle<edm::TriggerResults> triggerResultsHandle_;

      HLTConfigProvider hltConfig_;

      std::string processName_;
      std::string triggerName_;

      std::unique_ptr<KalmanVertexFitter> fitter;

      bool debug;

      TTree *eventTree;
      eventInfo *evInfo;
      TH1F* h_PFJet_PT_;
      TH1F* h_PFJet_eta_;
      TH1F* h_PFJet_HT_;
      TH1F* h_PFJet_N_;
      TH1F* h_muon_PT_;
      TH1F* h_muon_eta_;
      TH1F* h_track_PT_;
      TH1F* h_n_1_track_PT_;
      TH1F* h_n_1_track_npxl_;
      TH1F* h_n_1_track_nmin_;
      TH1F* h_n_1_track_nstl_;
      TH1F* h_n_1_Nsig_;
      TH1* h_Event_cutflow_;
      TH1* h_PFJet_cutflow_;
      TH1* h_Track_cutflow_;
      TH1* h_NTrackPerVtx_;
      TH1* h_dBV_;
      TH1* h_sigma_dBV_;
      TH2* h_vtx_XY_;
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
DVAnalyzer::DVAnalyzer(const edm::ParameterSet& iConfig)
 :
  bsToken_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamspotTag"))),
  jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))),
  muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))),
  //tracksToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  trackRefToken_(consumes<reco::TrackRefVector>(iConfig.getUntrackedParameter<edm::InputTag>("trackRefTag"))),
  trackRefMapToken_(consumes<DVAna::UnpackedCandidateTracksMap>(iConfig.getUntrackedParameter<edm::InputTag>("trackRefMapTag"))),
  trackRescaledToken_(consumes<DVAna::RescaledTrackMap>(iConfig.getUntrackedParameter<edm::InputTag>("trackRescaleMapTag"))),
  vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))),
  trackAssoToken_(consumes<DVAna::TrackAssociation>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssoTag"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag"))),
  processName_(iConfig.getUntrackedParameter<std::string>("processName")),
  triggerName_(iConfig.getUntrackedParameter<std::string>("triggerNameTag")),
  fitter(new KalmanVertexFitter()),
  debug(iConfig.getUntrackedParameter<bool>("debug"))

{
   //now do what ever initialization is needed
  evInfo = new eventInfo;

}


DVAnalyzer::~DVAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  int n_evt = iEvent.id().event();

  // assign approprate Handles and get the object collections
  //iEvent.getByToken(genJetToken_, genJetHandle_);

  iEvent.getByToken(bsToken_, bsHandle_);
  reco::BeamSpot bs = *bsHandle_.product();
  const reco::Vertex fake_bs_vtx(bs.position(), bs.covariance3D());
  VertexState theBeam = VertexState(bs);

  iEvent.getByToken(jetToken_, jetHandle_);
  auto jets = *jetHandle_.product();

  iEvent.getByToken(muonToken_, muonHandle_);
  auto muons = *muonHandle_.product();
  
  //iEvent.getByToken(tracksToken_, tracksHandle_);
  //auto tracks = *tracksHandle_.product();

  iEvent.getByToken(trackRefToken_, trackRefHandle_);
  auto tRef = *trackRefHandle_.product();

  //std::vector<reco::TrackRef> tRef;
  iEvent.getByToken(trackRefMapToken_, trackRefMapHandle_);
  auto tRefMap = *trackRefMapHandle_.product();

  iEvent.getByToken(trackRescaledToken_, trackRescaledHandle_);
  auto tRescaledMap = *trackRescaledHandle_.product();

  iEvent.getByToken(vtxToken_, vtxHandle_);
  const auto& vertex = (*vtxHandle_)[0];

  iEvent.getByToken(trackAssoToken_, trackAssoHandle_);
  const DVAna::TrackAssociation& vtx_tk_map = *trackAssoHandle_.product();

  // load Trigger paths
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle_);

  if (!triggerResultsHandle_.isValid()) {
    std::cout << "****Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }
  //auto trigNames = hltConfig_.triggerNames();
  //for (auto& tName:trigNames){
  //  std::cout << tName << std::endl;
  //}

  assert(triggerResultsHandle_->size()==hltConfig_.size());
  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIdx(hltConfig_.triggerIndex(triggerName_));
  assert(triggerIdx==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));
  if(triggerIdx>=ntrigs){
    std::cout << "Trigger path: " << triggerName_ << " not available! " << std::endl;
  }
  bool accept = triggerResultsHandle_->accept(triggerIdx);
  if(!accept){
    return;
  }

  // get transient track builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);


  initEventStructure();

  evInfo->evt = n_evt;

  // loop over objects in iEvent
  
  h_Event_cutflow_->Fill("allEvents",1);

  int pv = -1; // the vertex that passes the cuts and with the largest summed pT^2
  for (size_t iv=0; iv<vtxHandle_->size(); ++iv){
    //std::cout << "track size: " << (*vtxHandle_)[iv].tracksSize() << std::endl;
    reco::VertexRef v(vtxHandle_, iv);
    reco::TrackRefVector tks_vtx = vtx_tk_map[v];
    //std::cout << "tracks from vtx: " << v.key() << " : ";
    //for (auto tk_vtx=tks_vtx.begin(); tk_vtx!=tks_vtx.end(); ++tk_vtx){
    //  std::cout << tk_vtx->key() << " ";
    //}
    //std::cout <<  std::endl;

    if (!v->isFake() && v->ndof() > 4 && fabs(v->z()) <= 24 && v->position().rho() < 2) {
      pv = iv;
      break;
    }
  }

  if(pv==-1)  return;
  reco::VertexRef primary_vertex(vtxHandle_, pv);
  //std::cout << "PV: ";
  //print_vtx_track(*primary_vertex.get());

  h_Event_cutflow_->Fill("have good PV",1);

  //PFJets
  double pfJet_HT = 0;
  int pfJet_N = 0;

  for (const auto & jet:jets){
    h_PFJet_cutflow_->Fill("all jets",1);
    // PF Jet
    if(!jet.isPFJet()) continue;
    h_PFJet_cutflow_->Fill("PFJet",1);
    // pt>20GeV
    if(jet.pt()<=20) continue;
    h_PFJet_cutflow_->Fill("pT<20GeV",1);
    // |eta|<2.5
    if(std::abs(jet.eta())>2.5) continue;
    h_PFJet_cutflow_->Fill("|eta|<2.5",1);
    // # constituents > 1
    if( jet.numberOfDaughters() <= 1 ) continue;
    h_PFJet_cutflow_->Fill("#const>1",1);
    // Neutral hadron energy fraction < 0.9
    if( (jet.neutralHadronEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_NH<0.9",1);
    // Neutral EM energy fraction <0.9
    if( (jet.neutralEmEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_EM<0.9",1);
    // muon energy fraction < 0.8
    if( (jet.muonEnergyFraction()) > 0.8 ) continue;
    h_PFJet_cutflow_->Fill("f_E_MU<0.8",1);
    // cuts for |eta|<2.4
    if(std::abs(jet.eta())<2.4){
      // charged hadron energy fraction >0
      if( (jet.chargedHadronEnergyFraction()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("f_HE_c>0",1);
      // charged multiplicity>0
      if( (jet.chargedMultiplicity()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("M_c>0",1);
      // charged EM energy fraction < 0.8
      if( (jet.chargedEmEnergyFraction()) > 0.8 ) continue;
      h_PFJet_cutflow_->Fill("f_EM_c>0.8",1);
    }

    //get tracks from jets daughters
    //
    /*
    for (const reco::CandidatePtr& p : jet.daughterPtrVector()){
      auto findMap = tRefMap.find(p);
      if(findMap!=tRefMap.end()){
        reco::TrackRef tk = findMap->second;
        if(tk.isNonnull()){
          auto tFind = tRescaledMap.find(tk);
          if(tFind!=tRescaledMap.end()){
            tk = tFind->second;
            if(tk.isNonnull())
              tRef.push_back(tk);
          }
        }
      }
    }
    */


    int jet_vtx_idx = -1;
    std::vector<int> jetAssociation = jet_vtx_nTracks(jet, 0.0, vtx_tk_map, tRefMap, tRescaledMap);
    if(jetAssociation[0]!=-1)
      jet_vtx_idx = jetAssociation[0];
    else
      jet_vtx_idx = jetAssociation[1];
    //std::cout << "pv: " << pv << " jetass: " << jet_vtx_idx << std::endl;
    if(jet_vtx_idx!=pv) continue;

    h_PFJet_PT_->Fill(jet.pt());
    h_PFJet_eta_->Fill(jet.eta());
    if(jet.pt()>=40){
      pfJet_HT += jet.pt();
    }
    ++pfJet_N;
  }
  // Event preselection:
  // HT > 1200 GeV
  if(pfJet_HT<=1200) return;
  h_Event_cutflow_->Fill("HT>1200", 1);
  h_PFJet_HT_->Fill(pfJet_HT);
  // At least four jet satisfying above criteria
  if(pfJet_N<4) return;
  h_Event_cutflow_->Fill("NJet>=4", 1);
  h_PFJet_N_->Fill(pfJet_N);

  //muons
  for (const auto & muon:muons){
    if(muon.pt()<27) continue;
    if(std::abs(muon.eta())>2.4) continue;
    if(muon::isTightMuon(muon, vertex)){
      h_muon_PT_->Fill(muon.pt());
      h_muon_eta_->Fill(muon.eta());
    }
  }

  //tracks

  std::vector<reco::TransientTrack> tTracks;
  std::map<reco::TrackRef, size_t> seed_track_ref_map;

  for (auto tk = tRef.begin(); tk!=tRef.end();++tk){
    auto track = *tk;
    //calculate basic variables for tracks
    const double pt = track->pt();
    const double dxybs = track->dxy(bs);
    const double dxypv = track->dxy(primary_vertex->position());
    const double dxyerr = track->dxyError();
    //const double rescaled_dxyerr = rs.rescaled_tk.dxyError();
    const double sigmadxybs = dxybs / dxyerr;
    //const double rescaled_sigmadxybs = dxybs / rescaled_dxyerr;
    const double sigmadxypv = dxypv / dxyerr;
    const int nhits = track->hitPattern().numberOfValidHits();
    const int npxhits = track->hitPattern().numberOfValidPixelHits();
    //const int nsthits = track->hitPattern().numberOfValidStripHits();
    const int npxlayers = track->hitPattern().pixelLayersWithMeasurement();
    const int nstlayers = track->hitPattern().stripLayersWithMeasurement();
    int min_r = 2000000000;
    for (int i = 1; i <= 4; ++i)
      if (track->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,i)) {
        min_r = i;
        break;
      }

    if(fabs(dxybs)<=0 || dxyerr>=1e9 || fabs(sigmadxypv)<=0 || nhits<0 || npxhits<0) continue;

    if( min_r<=1 && npxlayers>=2 && nstlayers >= 6 && fabs(sigmadxybs) > 4)
      h_n_1_track_PT_->Fill(pt);
    if( min_r<=1 && pt>1 && nstlayers >= 6 && fabs(sigmadxybs) > 4)
      h_n_1_track_npxl_->Fill(npxlayers);
    if( pt>1 && npxlayers>=2 && nstlayers >= 6 && fabs(sigmadxybs) > 4)
      h_n_1_track_nmin_->Fill(min_r);
    if( min_r<=1 && pt>1 && fabs(sigmadxybs) > 4 && npxlayers>=2 )
      h_n_1_track_nstl_->Fill(nstlayers);
    if( min_r<=1 && npxlayers>=2 && pt>1 && nstlayers >= 6)
      h_n_1_Nsig_->Fill(fabs(sigmadxybs));

    if( min_r>1 || npxlayers<2 || nstlayers < 6 || fabs(sigmadxybs) <= 4 || pt<=1)
      continue;


    //// track pt > 1 GeV
    //if(track->pt()<=1) continue;
    //h_Track_cutflow_->Fill("pt>1GeV", 1);
    //
    ////rmin = 1
    //if(!track->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1)) continue;
    //h_Track_cutflow_->Fill("rmin==1", 1);

    //// npxl >= 2
    //if(track->hitPattern().pixelLayersWithMeasurement()<2) continue;
    //h_Track_cutflow_->Fill("npxl>=2", 1);

    //// nstl >= 6
    //if(track->hitPattern().stripLayersWithMeasurement() < 6) continue;
    //h_Track_cutflow_->Fill("nstl>=6", 1);

    //// |dxy|/sigma_dxy > 4 
    //if( (std::abs(track->dxy())/track->dxyError()) <= 4 ) continue;
    //h_Track_cutflow_->Fill("dxy/sigma_dxy>4", 1);

    //h_track_PT_->Fill(track->pt());

    //reco::TransientTrack t(track, &(*bFieldHandle));
    reco::TransientTrack t = (*theB).build(track);
    tTracks.push_back(t);
    seed_track_ref_map[track] = tTracks.size()-1;
  }

  /*
  for (const auto & track:tracks){
    h_Track_cutflow_->Fill("all tracks", 1);

    // track pt > 1 GeV
    if(track.pt()<=1) continue;
    h_Track_cutflow_->Fill("pt>1GeV", 1);
    
    //rmin = 1
    if(track.lostInnerHits()!=-1) continue;
    h_Track_cutflow_->Fill("rmin==1", 1);

    // npxl >= 2
    if(track.pixelLayersWithMeasurement()<2) continue;
    h_Track_cutflow_->Fill("npxl>=2", 1);

    // nstl >= 6
    if(track.stripLayersWithMeasurement() < 6) continue;
    h_Track_cutflow_->Fill("nstl>=6", 1);

    // |dxy|/sigma_dxy > 4 
    if( (std::abs(track.dxy())/track.dxyError()) <= 4 ) continue;
    h_Track_cutflow_->Fill("dxy/sigma_dxy>4", 1);

    h_track_PT_->Fill(track.pt());

    if(track.hasTrackDetails()){
      reco::TransientTrack t(track.pseudoTrack(), &(*bFieldHandle));
      tTracks.push_back(t);
      seed_track_ref_map[track] = tTracks.size()-1;
      //reco::TransientTrack t = (*theB).build(track.bestTrack());
    }
  }*/

  // refit vertices

  if(debug){
    std::cout << "Transient Track Size: " << tTracks.size() << std::endl;
  }
  if(tTracks.size()<2) return;
  std::vector<reco::Vertex> v_refit;
  for(auto tTrack = tTracks.cbegin(); tTrack!=tTracks.cend()-1; ++tTrack){
    const auto& tT1 = *tTrack;
    for(auto tTrack2 = tTrack+1; tTrack2!=tTracks.cend(); ++tTrack2){
      const auto& tT2 = *tTrack2;
      std::vector<reco::TransientTrack> tTrackVector({tT1, tT2});
      TransientVertex tv = fitter->vertex(tTrackVector); 
      if(tv.isValid() && (tv.normalisedChiSquared())<5){
        v_refit.push_back(reco::Vertex(tv));
      }
    }
  }
  if(debug){
    std::cout << "seed vtx refit done, number of vtx: " << v_refit.size() << std::endl;
    for (auto iv=v_refit.begin(); iv!=v_refit.end(); ++iv){
      std::cout << "vtx: ";
      print_vtx_track(*iv);
    }
  }
  VertexDistance3D vtxDistance;
  std::vector<reco::Vertex>::iterator v[2];
  for(v[0]=v_refit.begin(); v[0]!=v_refit.end(); ++v[0]) {
    if(debug){
      std::cout << "v0: ";
      print_vtx_track(*v[0]);
    }
    const auto& v1 = *v[0];
    track_set track[2]; //tracks in vertices
    track[0] = vertexTracks(v1);
    if(track[0].size()<2) {
      v[0] = v_refit.erase(v[0])-1;
      continue;
    }

    bool duplicate = false;
    bool merge = false;
    bool refit = false;
    track_set tk_rm[2]; //tracks that are going to be removed from vertices

    for(v[1]=v[0]+1; v[1]!=v_refit.end(); ++v[1]){
      const auto& v2 = *v[1];

      track[1] = vertexTracks(v2);
      if(track[1].size()<2){
        v[1] = v_refit.erase(v[1])-1;
        continue;
      }

      // if one vertex have all tracks included in the other vertex
      // remove the second vertex
      if(is_track_subset(track[0],track[1])){
        duplicate = true;
        if(debug){
          std::cout << "duplicate, delete v1: ";
          print_vtx_track(*v[1]);
        }
        break;
      }

      auto sharedTracks = SharedTracks(track[0], track[1]);
      if(sharedTracks.size()==0) continue;
      if(debug){
        std::cout << "v1: ";
        print_vtx_track(*v[1]);
      }
      //std::cout << "# Stracks: " << sharedTracks.size() << std::endl;
      //std::cout << "t1: " << track[0].size() << std::endl;
      //std::cout << "t2: " << track[1].size() << std::endl;
      if(vtxDistance.distance(v1, v2).significance()<4){
        // combine the two vertex tracks and refit
        merge = true;
        if(debug){
          std::cout << "dist_sig: " << vtxDistance.distance(v1, v2).significance() << " will merge" << std::endl;
        }
      }
      else{
        refit = true;
        if(debug){
          std::cout << "will refit" << std::endl;
        }
      }
      for (auto tkIter=sharedTracks.begin(); tkIter!=sharedTracks.end(); ++tkIter) {
        auto tk = *tkIter;
        TransientTrack ttk = tTracks[seed_track_ref_map[tk]];
        auto ip3D_pair1 = IPTools::absoluteImpactParameter3D(ttk, v1);
        auto ip3D_pair2 = IPTools::absoluteImpactParameter3D(ttk, v2);
        double ip3D_1_sig = ip3D_pair1.second.significance();
        double ip3D_2_sig = ip3D_pair2.second.significance();
        ip3D_pair1.first = ip3D_pair1.first && (ip3D_pair1.second.significance()<5);
        ip3D_pair2.first = ip3D_pair2.first && (ip3D_pair2.second.significance()<5);
        bool rm_from_1 = !ip3D_pair1.first;
        bool rm_from_2 = !ip3D_pair2.first;

        if(debug){
          std::cout << "tk " << tk.key() << " sig1: " << ip3D_1_sig << " sig2: " << ip3D_2_sig << std::endl;
        }
        if(ip3D_1_sig<1.5 && ip3D_2_sig<1.5){
          if(v1.tracksSize()>v2.tracksSize()){
            rm_from_2 = true;
          }
          else{
            rm_from_1 = true;
          }
        }
        else if(ip3D_1_sig<=ip3D_2_sig){
          rm_from_2 = true;
        }
        else{
          rm_from_1 = true;
        }
        if(rm_from_1) tk_rm[0].insert(tk);
        if(rm_from_2) tk_rm[1].insert(tk);
        break; // remove one track at a time
      }
      //std::cout << "tks to rm from 0: " << tk_rm[0].size() << std::endl;
      //std::cout << "tks to rm from 1: " << tk_rm[1].size() << std::endl;
      break;
    }
    if(duplicate){
      //std::cout << "duplicate" << std::endl;
      v_refit.erase(v[1]);
    }
    else if(merge){
      //std::cout << "merge"<<track[0].size() <<"  " << track[1].size() << std::endl;
      track_set combinedSet;
      for (auto tk:track[0]){
        combinedSet.insert(tk);
      }
      for (auto tk:track[1]){
        combinedSet.insert(tk);
      }
      //std::cout << "combined: " << combinedSet.size() << std::endl;
      std::vector<TransientTrack> cT;
      for (auto tk:combinedSet){
        cT.push_back(tTracks[seed_track_ref_map[tk]]);
      }

      reco::VertexCollection new_vertices;
      for (const TransientVertex& tv:refitVtx(cT))
        new_vertices.push_back(reco::Vertex(tv));

      if(new_vertices.size() > 1){
        if(debug){
          std::cout << "get 2 refitted vtx..." << std::endl; 
        }
        assert(new_vertices.size() == 2);
        *v[1] = reco::Vertex(new_vertices[1]);
        *v[0] = reco::Vertex(new_vertices[0]);
      }
      else if(new_vertices.size()==1 && vertexTracks(new_vertices[0])==combinedSet){
        v_refit.erase(v[1]);
        *v[0] = reco::Vertex(new_vertices[0]);
      }
      else{
        if(debug)
          std::cout << "merge didn't work, trying refits..." << std::endl;
        refit = true;
      }
      //auto v_new = refitVtx(cT);
      //if(v_new.size()==1 && vertexTracks(v_new[0])==combinedSet){
      //  v_refit.erase(v[1]);
      //  *v[0] = reco::Vertex(v_new[0]);
      //}
      //else refit = true;
      //std::cout << "merge done" << std::endl;
    }
    if(refit){
      bool erase[2] = {false};
      reco::Vertex vsave[2] = {*v[0], *v[1]};
      //std::cout << "refit" << std::endl;
      for (int i = 0; i<2; ++i){
        if(debug){
          std::cout << tk_rm[i].size() << " track to remove from " << i << " ";
          for(auto t:tk_rm[i]){
            std::cout << t.key() << " ";
          }
          std::cout << std::endl;
        }
        if(tk_rm[i].empty()) continue;
        std::vector<reco::TransientTrack> ttks;
        for (auto tk:track[i]){
          if(tk_rm[i].count(tk)==0)
            ttks.push_back(tTracks[seed_track_ref_map[tk]]);
        }

        reco::VertexCollection new_vertices;
        for ( const TransientVertex& tv : refitVtx(ttks))
          new_vertices.push_back(reco::Vertex(tv));
        if(debug){
          std::cout << "refit " << new_vertices.size() << " vertices " << std::endl;
          for (const auto& nv: new_vertices){
            std::cout << "vtx: ";
            print_vtx_track(nv);
          }
        }
        if(new_vertices.size()==1)
          *v[i] = new_vertices[0];
        else
          erase[i] = true;

        //auto v_new = refitVtx(ttks);
        //if(v_new.size()==1){
        //  *v[i] = reco::Vertex(v_new[0]);
        //  if(debug){
        //    std::cout << "refit: " << i << "  tks: ";
        //    print_vtx_track(*v[i]);
        //  }
        //}
        //else{
        //  if(debug){
        //    std::cout << "refit rm: " << i << "  tks: ";
        //    print_vtx_track(*v[i]);
        //  }
        //  v_refit.erase(v[i]);
        //}
      }
      if(erase[1]) v_refit.erase(v[1]);
      if(erase[0]) v_refit.erase(v[0]);
      //std::cout << "refit done" << std::endl;
    }

    // start over if vertices are changed
    if(duplicate || merge || refit){
      //std::cout << "start over" << std::endl;
      v[0] = v_refit.begin()-1;
    }
  }

  // merge vertices that are still too close
  //
  /*
  if(debug){
    std::cout << "merge vtx that are still close" << std::endl;
  }
  for (v[0] = v_refit.begin(); v[0]!=v_refit.end(); ++v[0]){
    if(debug){
      std::cout << "v1: ";
      print_vtx_track(*v[0]);
    }
    bool merge = false;
    for (v[1] = v[0]+1; v[1]!=v_refit.end(); ++v[1]){
      if(debug){
        std::cout << "v2: ";
        print_vtx_track(*v[1]);
        std::cout << "dist_sig: " << vtxDistance.distance(*v[0], *v[1]).significance() << std::endl;
      }
      if(vtxDistance.distance(*v[0], *v[1]).significance()<4){
        if(debug){
          std::cout << "dist_sig: " << vtxDistance.distance(*v[0], *v[1]).significance() << "merging..." << std::endl;
        }
        merge = true;
        break;
      }
    }
    if(merge){
      track_set combinedSet;
      for (auto tk:vertexTracks(*v[0])){
        combinedSet.insert(tk);
      }
      for (auto tk:vertexTracks(*v[1])){
        combinedSet.insert(tk);
      }
      std::vector<TransientTrack> cT;
      for (auto tk:combinedSet){
        cT.push_back(tTracks[seed_track_ref_map[tk]]);
      }
      auto v_new = refitVtx(cT);
      if(v_new.size()==1 && vertexTracks(v_new[0])==combinedSet){
        v_refit.erase(v[1]);
        *v[0] = reco::Vertex(v_new[0]);
        if(debug){
          std::cout << "merged vertex: ";
          print_vtx_track(*v[0]);
        }
        v[0] = v_refit.begin()-1;
      }
    }
  }
  */

  if(debug){
    std::cout << "revertex done " << std::endl;
    for (auto iv=v_refit.begin(); iv!=v_refit.end(); ++iv){
      std::cout << "vtx: ";
      print_vtx_track(*iv);
      for (auto iv2=iv+1; iv2!=v_refit.end(); ++iv2){
        auto sharedTracks_after = SharedTracks(vertexTracks(*iv), vertexTracks(*iv2));
        if(sharedTracks_after.size()!=0){
          std::cout << "warning!!! still shared tracks." << std::endl;
        }
      }
    }
  }

  for (v[0]=v_refit.begin(); v[0]!=v_refit.end(); ++ v[0]){
    if(debug){
      std::cout << "vtx_ref: ";
      print_vtx_track(*v[0]);
    }
    const track_vec tks = vertex_track_vec(*v[0]);
    const size_t ntks = tks.size();
    if(ntks<3)
      continue;

    std::vector<TransientTrack> ttks(ntks-1);
    for(size_t i = 0; i<ntks; ++i){
      for(size_t j = 0; j<ntks; ++j){
        if(j!=i)
          //ttks[j-(j>=i)] = TransientTrack(tks[j], &(*bFieldHandle));
          ttks[j-(j>=i)] = (*theB).build(tks[j]);
      }
      //std::cout << "rm" << std::endl;
      reco::Vertex v_rm_track(TransientVertex(fitter->vertex(ttks)));
      const double dz = std::abs(v_rm_track.z()-v[0]->z());
      if(debug){
        std::cout << "rm: " << tks[i].key() << " dz: " << dz << " chi2: " << v_rm_track.chi2() << " v: ";
        print_vtx_track(v_rm_track);

      }
      if(v_rm_track.chi2()<0 || 
          dz > 0.005){
        if(debug){
          std::cout << "replacing..." << std::endl;
        }
        *v[0] = v_rm_track;
        --v[0];
        break;
      }
    }
  }

  if(debug){
    std::cout << "remove done " << std::endl;
    for (auto iv=v_refit.begin(); iv!=v_refit.end(); ++iv){
      std::cout << "vtx: ";
      print_vtx_track(*iv);
      for (auto iv2=iv+1; iv2!=v_refit.end(); ++iv2){
        auto sharedTracks_after = SharedTracks(vertexTracks(*iv), vertexTracks(*iv2));
        if(sharedTracks_after.size()!=0){
          std::cout << "warning!!! still shared tracks." << std::endl;
        }
      }
    }
  }

  VertexDistanceXY dBVCalc = VertexDistanceXY();
  if(debug){
    std::cout << "final number of vtx: " << v_refit.size() << std::endl;
  }
  for (auto vIter=v_refit.begin(); vIter!=v_refit.end();++vIter){
    //std::cout << "Beam x: " << bs.x(vIter->z()) << " y: " << bs.y(vIter->z()) << std::endl;
    if(debug){
      std::cout << "final vtx with " << vIter->tracksSize() << " tracks: ";
      print_vtx_track(*vIter);
    }
    //if( vIter->normalizedChi2()<=0 || vIter->normalizedChi2()>5) continue;
    //if(vIter->tracksSize()<5) continue;

    const double vtxZ = vIter->z();
    const double bsxAtZ = bs.x(vtxZ);
    const double bsyAtZ = bs.y(vtxZ);
    double dBV_calc = mag(vIter->x()-bsxAtZ, vIter->y()-bsyAtZ);
    //Measurement1D dBV = dBVCalc.distance(*vIter, theBeam);
    Measurement1D dBV = dBVCalc.distance(*vIter, fake_bs_vtx);
    //if(dBV.value()<=0.01) continue;
    //if(dBV.error()>=0.0025) continue;
    //if(dBV.value()>2.09) continue;
    //evInfo->evt.push_back(n_evt);
    evInfo->vtx_track_size.push_back(vIter->tracksSize());
    evInfo->vtx_dBV.push_back(dBV_calc);
    evInfo->vtx_sigma_dBV.push_back(dBV.error());
    evInfo->vtx_x.push_back(vIter->x());
    evInfo->vtx_y.push_back(vIter->y());
    evInfo->vtx_z.push_back(vIter->z());

    h_dBV_->Fill(dBV.value());
    h_sigma_dBV_->Fill(dBV.error());
    h_NTrackPerVtx_->Fill(vIter->tracksSize());
    h_vtx_XY_->Fill(vIter->x(), vIter->y());
  }

  eventTree->Fill();


  //for(auto vIter=v_refit.begin(); vIter!=v_refit.end(); ++vIter){
  //  auto v_after = *vIter;
  //  if(! (v_after.isValid() && (v_after.totalChiSquared()/v_after.degreesOfFreedom())<5)){
  //    v_refit.erase(vIter);
  //    --vIter;
  //  }
  //  std::cout << "# tks: " << v_after.tracks().size() << std::endl;
  //}
  //std::cout << "# vtx: " << v_refit.size() << std::endl;

    //Handle<TrackCollection> tracks;
    //iEvent.getByToken(tracksToken_, tracks);
    //for(TrackCollection::const_iterator itTrack = tracks->begin();
    //    itTrack != tracks->end();
    //    ++itTrack) {
    //  // do something with track parameters, e.g, plot the charge.
    //  // int charge = itTrack->charge();
    //}

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
DVAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");

  eventTree = fileService->make<TTree>( "tree_DV", "tree_DV" );
  eventTree->Branch( "evt",                  &evInfo->evt);
  eventTree->Branch( "vtx_track_size",       &evInfo->vtx_track_size);
  eventTree->Branch( "vtx_dBV",              &evInfo->vtx_dBV);
  eventTree->Branch( "vtx_sigma_dBV",        &evInfo->vtx_sigma_dBV);  
  eventTree->Branch( "vtx_x",                &evInfo->vtx_x);
  eventTree->Branch( "vtx_y",                &evInfo->vtx_y);
  eventTree->Branch( "vtx_z",                &evInfo->vtx_z);
  //eventTree->Branch( "bsAtZ_x",              &evInfo->bsAtZ_x);
  //eventTree->Branch( "bsAtZ_y",              &evInfo->bsAtZ_y);

  h_PFJet_PT_ = fileService->make<TH1F>("h_PFJet_PT","h_PFJet_PT",100,0,100);
  h_PFJet_eta_ = fileService->make<TH1F>("h_PFJet_eta","h_PFJet_eta",30,-3,3);
  h_PFJet_HT_ = fileService->make<TH1F>("h_PFJet_HT","h_PFJet_HT",35,0,3500);
  h_PFJet_N_ = fileService->make<TH1F>("h_PFJet_N","h_PFJet_N",30,0,30);
  h_muon_PT_ = fileService->make<TH1F>("h_muon_PT","h_muon_PT",100,0,100);
  h_muon_eta_ = fileService->make<TH1F>("h_muon_eta","h_muon_eta",100,-5,5);
  h_track_PT_ = fileService->make<TH1F>("h_track_PT","h_track_PT",100,0,25);
  h_n_1_track_PT_ = fileService->make<TH1F>("h_n_1_track_PT","h_n_1_track_PT",200,0,20);
  h_n_1_track_npxl_ = fileService->make<TH1F>("h_n_1_track_npxl","h_n_1_track_npxl",20,0,20);
  h_n_1_track_nmin_ = fileService->make<TH1F>("h_n_1_track_nmin","h_n_1_track_nmin",20,0,20);
  h_n_1_track_nstl_ = fileService->make<TH1F>("h_n_1_track_nstl","h_n_1_track_nstl",20,0,20);
  h_n_1_Nsig_ = fileService->make<TH1F>("h_n_1_Nsig","h_n_1_Nsig",100,0,10);
  h_Event_cutflow_ = fileService->make<TH1D>("h_Event_cutflow","h_Event_cutflow",10,0,10);
  h_PFJet_cutflow_ = fileService->make<TH1D>("h_PFJet_cutflow","h_PFJet_cutflow",10,0,10);
  h_Track_cutflow_ = fileService->make<TH1D>("h_Track_cutflow","h_Track_cutflow",10,0,10);
  h_NTrackPerVtx_ = fileService->make<TH1D>("h_NTrackPerVtx","h_NTrackPerVtx",50,0,50);
  h_dBV_ = fileService->make<TH1D>("h_dBV","h_dBV",100,0,0.4);
  h_sigma_dBV_ = fileService->make<TH1D>("h_sigma_dBV","h_sigma_dBV",500,0,0.05);
  h_vtx_XY_ = fileService->make<TH2D>("h_vtx_XY","h_vtx_XY",200,-4,4,200,-4,4);
}

void DVAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if(hltConfig_.init(iRun, iSetup, processName_, changed)){
    if(changed) {
      const unsigned int n(hltConfig_.size());
      unsigned int triggerIdx(hltConfig_.triggerIndex(triggerName_));
      if (triggerIdx>=n) {
        std::cout << "TriggerName: " << triggerName_ << " not available in config!" << std::endl;
      }
    }
  }
  else {
    std::cout << "Warning, didn't find trigger process HLT,\t" << processName_ << std::endl;
  }
}

std::vector<TransientVertex> DVAnalyzer::refitVtx(std::vector<reco::TransientTrack>& ttks) {
  if (ttks.size() < 2)
    return std::vector<TransientVertex>();
  std::vector<TransientVertex> v(1, fitter->vertex(ttks));
  if (v[0].normalisedChiSquared() > 5)
    return std::vector<TransientVertex>();
  return v;
}

reco::TrackRefVector DVAnalyzer::SharedTracks(track_set t1, track_set t2){
  reco::TrackRefVector sT;
  for(auto tk:t1){
    if(t2.count(tk)>0){
      sT.push_back(tk);
    }
  }
  return sT;
}

std::vector<TransientTrack> DVAnalyzer::combineTracks(track_set t1, track_set t2, std::map<reco::TrackRef, size_t> seed_track_ref_map, std::vector<TransientTrack> tTracks) {
  track_set combinedSet;
  for (auto tk:t1){
    combinedSet.insert(tk);
  }
  for (auto tk:t2){
    combinedSet.insert(tk);
  }
  std::vector<TransientTrack> tks;
  for (auto tk:combinedSet){
    tks.push_back(tTracks[seed_track_ref_map[tk]]);
  }
  return tks;
}


track_set DVAnalyzer::vertexTracks(const reco::Vertex v) const
{
  track_set result;
  for(auto t=v.tracks_begin(); t!=v.tracks_end(); ++t){
    result.insert(t->castTo<reco::TrackRef>());
  }
  return result;
}

track_vec DVAnalyzer::vertex_track_vec(const reco::Vertex& v) const {
    track_set s = vertexTracks(v);
    return track_vec(s.begin(), s.end());
}

bool DVAnalyzer::is_track_subset(const track_set& a, const track_set& b) const {
  bool is_subset = true;
  const track_set& s = a.size() <= b.size() ? a : b;
  const track_set& l = a.size() <= b.size() ? b : a;

  for (auto t : s)
    if (l.count(t)<1){
      is_subset = false;
      break;
    }
  return is_subset;
}

void DVAnalyzer::print_vtx_track(const reco::Vertex v) const {
  for (auto i = v.tracks_begin(); i!=v.tracks_end(); ++i){
    std::cout << i->key() << "  ";
  }
  std::cout << std::endl;
}

track_set DVAnalyzer::getJetTrack(const pat::Jet& jet, DVAna::UnpackedCandidateTracksMap& tRefMap, DVAna::RescaledTrackMap& tRescaledMap){
  track_set tRef;
  for (const reco::CandidatePtr& p : jet.daughterPtrVector()){
    auto findMap = tRefMap.find(p);
    if(findMap!=tRefMap.end()){
      reco::TrackRef tk = findMap->second;
      if(tk.isNonnull()){
        auto tFind = tRescaledMap.find(tk);
        if(tFind!=tRescaledMap.end()){
          tk = tFind->second;
          if(tk.isNonnull())
            tRef.insert(tk);
        }
      }
    }
  }
  return tRef;
}

std::vector<int> DVAnalyzer::jet_vtx_nTracks(const pat::Jet& jet, double track_pt_min, DVAna::TrackAssociation vtx_tk_map, DVAna::UnpackedCandidateTracksMap& tRefMap, DVAna::RescaledTrackMap& tRescaledMap){
  std::vector<int> jetAssociation(4,-1); 
  //jet association vector:
  // 0: association by the number of tracks from vertices
  // 1: association by the cos angle between SV and PV
  // 2: association by the miss distance
  track_set jet_tracks = getJetTrack(jet, tRefMap, tRescaledMap);
  std::set<int> jet_tracks_key;

  //std::cout << "jet tracks: ";
  for(auto tjet=jet_tracks.begin(); tjet!=jet_tracks.end(); ++tjet){
    jet_tracks_key.insert(tjet->key());
    //std::cout << tjet->key() << " ";
  }
  //std::cout << std::endl;

  const reco::SecondaryVertexTagInfo* jet_tag = jet.tagInfoSecondaryVertex("secondaryVertexMaxDR6p0");
  //const reco::SecondaryVertexTagInfo* jet_tag = jet.tagInfoSecondaryVertex();
  const reco::Vertex* jet_tag_vtx = 0;
  TVector3 jet_tag_vtx_pos;
  if(jet_tag && jet_tag->nVertices() > 0){
    jet_tag_vtx = &jet_tag->secondaryVertex(0);
    jet_tag_vtx_pos.SetXYZ(jet_tag_vtx->x(), jet_tag_vtx->y(), jet_tag_vtx->z());
  }
  const TVector3 jet_mom_dir = TVector3(jet.px(), jet.py(), jet.pz()).Unit();
  int ntracks_max = 0;
  double cos_angle_min = 1e9;
  Measurement1D miss_dist_min(1e9, 1);

  const double max_cos_angle_diff = 2;
  const double max_miss_dist = 1e6;
  const double max_miss_sig = 2;
  for (size_t ivtx=0; ivtx<vtxHandle_->size(); ++ivtx){

    const double min_vertex_track_weight = 0.5;

    const reco::VertexRef vtx = reco::VertexRef(vtxHandle_, ivtx);
    //const reco::Vertex& vtx = vtxHandle_->at(ivtx);
    if( vtx_tk_map.find(vtx) == vtx_tk_map.end() ) continue;
    reco::TrackRefVector v_tks = vtx_tk_map[vtx];

    int ntracks = 0;
    double cos_angle = 1e9;
    Measurement1D miss_dist(1e9, 1);

    //std::cout << "vtx " << ivtx << " has tracks: ";
    //for(auto itkp = v_tks.begin(); itkp!=v_tks.end(); ++itkp){
    //  std::cout << itkp->key() << " ";
    //}
    //std::cout << std::endl;
    //std::cout << "vtx " << ivtx << " has tracks: ";
    for(auto itk = v_tks.begin(); itk!=v_tks.end(); ++itk){
      //if(vtx->trackWeight(*itk)>=min_vertex_track_weight){
      if(1){
        //reco::TrackRef tk = (*itk)->castTo<reco::TrackRef>();
        reco::TrackRef tk = *itk;
        //std::cout << "ref: " << tk.key() << " pt: " << tk->pt() << " count: " << jet_tracks.count(tk) << std::endl;
        if(jet_tracks_key.count(tk.key())>0 && tk->pt()>track_pt_min){
          //std::cout << "Found matched track " << tk.key() << std::endl;
          ++ntracks;
        }
      }
    }
    //std::cout << std::endl;
    if(jet_tag_vtx){
      const TVector3 sv_to_tv = jet_tag_vtx_pos - TVector3(vtx->x(), vtx->y(), vtx->z());
      cos_angle = sv_to_tv.Dot(jet_mom_dir) / sv_to_tv.Mag();
      // miss distance is magnitude of (jet direction cross (tv - sv))
      // to calculate uncertainty, use |n X d|^2 = (|n||d|)^2 - (n . d)^2
      const TVector3& n = jet_mom_dir;
      const TVector3& d = sv_to_tv;
      const double n_dot_d = n.Dot(d);
      const TVector3 n_cross_d = n.Cross(d);
      typedef math::VectorD<3>::type vec_t;
      typedef math::ErrorD<3>::type mat_t;
      vec_t jacobian(2*d.x() - 2*n_dot_d*n.x(),
                     2*d.y() - 2*n_dot_d*n.y(),
                     2*d.z() - 2*n_dot_d*n.z());
      mat_t sv_to_tv_cov_matrix = vtx->covariance() + jet_tag_vtx->covariance();
      double sigma_f2 = sqrt(ROOT::Math::Similarity(jacobian, sv_to_tv_cov_matrix));
      double miss_dist_value = n_cross_d.Mag();
      double miss_dist_err = sigma_f2 / 2 / miss_dist_value;
      miss_dist = Measurement1D(miss_dist_value, miss_dist_err);
      //std::cout << "miss_dist: " << miss_dist_value << " cos_angle: " << cos_angle << std::endl;
    }
    if(ntracks>=1 && ntracks>ntracks_max){
      ntracks_max = ntracks;
      jetAssociation[0] = ivtx;
    }
    if(fabs(cos_angle - 1) <= max_cos_angle_diff && fabs(cos_angle - 1) < fabs(cos_angle_min - 1)){
      cos_angle_min = cos_angle;
      jetAssociation[1] = ivtx;
    }
    if(miss_dist.value() <= max_miss_dist && miss_dist.significance() <= max_miss_sig && miss_dist.value() < miss_dist_min.value()){
      miss_dist_min = miss_dist;
      jetAssociation[2] = ivtx;
    }
  }
  return jetAssociation;
}

double DVAnalyzer::mag(double x, double y){
  return sqrt(x*x + y*y);
}

void
DVAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{}

void DVAnalyzer::initEventStructure()
{
  evInfo->evt=-1;
  evInfo->vtx_track_size.clear();
  evInfo->vtx_dBV.clear();
  evInfo->vtx_sigma_dBV.clear();
  evInfo->vtx_x.clear();
  evInfo->vtx_y.clear();
  evInfo->vtx_z.clear();
}

// ------------ method called once each job just after ending the event loop  ------------
void
DVAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DVAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DVAnalyzer);
