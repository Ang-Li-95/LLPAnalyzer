// -*- C++ -*-
//
// Package:    LLPAnalyzer/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc LLPAnalyzer/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ang Li
//         Created:  Fri, 29 May 2020 02:07:02 GMT
//
//


// system include files
#include <memory>
#include <cassert>
#include <map>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TH1.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      //virtual void beginRun(edm::Run const &, edm::EventSetup const&) override;
      void beginRun(edm::Run const& , edm::EventSetup const&) override ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endRun(edm::Run const& , edm::EventSetup const&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      //edm::EDGetTokenT<std::vector<reco::GenJet>> genJetToken_;
      //edm::Handle<std::vector<reco::GenJet>> genJetHandle_;
      edm::EDGetTokenT<reco::PFJetCollection> pfJetToken_;
      edm::Handle<reco::PFJetCollection> pfJetHandle_;
      edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
      edm::Handle<std::vector<reco::Muon>> muonHandle_;
      edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
      edm::Handle<std::vector<reco::Track>> tracksHandle_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
      edm::Handle<std::vector<reco::Vertex>> vtxHandle_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::Handle<edm::TriggerResults> triggerResultsHandle_;

      HLTConfigProvider hltConfig_;

      std::string processName_;
      std::string triggerName_;

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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  //genJetToken_(consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetTag")))
  pfJetToken_(consumes<reco::PFJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfJetTag"))),
  muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))),
  tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag"))),
  processName_(iConfig.getUntrackedParameter<std::string>("processName")),
  triggerName_(iConfig.getUntrackedParameter<std::string>("triggerNameTag"))
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed

}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // assign approprate Handles and get the object collections
  //iEvent.getByToken(genJetToken_, genJetHandle_);
  iEvent.getByToken(pfJetToken_, pfJetHandle_);
  auto pfJets = *pfJetHandle_.product();

  iEvent.getByToken(muonToken_, muonHandle_);
  auto muons = *muonHandle_.product();
  
  iEvent.getByToken(tracksToken_, tracksHandle_);
  auto tracks = *tracksHandle_.product();

  iEvent.getByToken(vtxToken_, vtxHandle_);
  const auto& vertex = (*vtxHandle_)[0];

  // load Trigger paths
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle_);

  if (!triggerResultsHandle_.isValid()) {
    std::cout << "****Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }

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

  //std::vector<std::string> triggerNames = hltConfig_.triggerNames();
  //std::cout << "# paths: " << hltConfig_.size() << "  " << triggerResultsHandle_->size() << std::endl;
  //for (unsigned int iPath=0; iPath<triggerNames.size();++iPath){
  //  std::string pathName = triggerNames[iPath];
  //  std::cout << "HLT path name: " << pathName << std::endl;
  //}

  // loop over objects in iEvent
  
  //PFJets
  double pfJet_HT = 0;
  int pfJet_N = 0;

  h_Event_cutflow_->Fill("allEvents",1);

  for (const auto & pfJet:pfJets){
    h_PFJet_cutflow_->Fill("all jets",1);
    // pt>20GeV
    if(pfJet.pt()<=20) continue;
    h_PFJet_cutflow_->Fill("pT<20GeV",1);
    // |eta|<2.5
    if(std::abs(pfJet.eta())>2.5) continue;
    h_PFJet_cutflow_->Fill("|eta|<2.5",1);
    // # constituents > 1
    if( ((pfJet.getPFConstituents()).size()) <= 1 ) continue;
    h_PFJet_cutflow_->Fill("#const>1",1);
    // Neutral hadron energy fraction < 0.9
    if( (pfJet.neutralHadronEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_NH<0.9",1);
    // Neutral EM energy fraction <0.9
    if( (pfJet.neutralEmEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_EM<0.9",1);
    // muon energy fraction < 0.8
    if( (pfJet.muonEnergyFraction()) > 0.8 ) continue;
    h_PFJet_cutflow_->Fill("f_E_MU<0.8",1);
    // cuts for |eta|<2.4
    if(std::abs(pfJet.eta())<2.4){
      // charged hadron energy fraction >0
      if( (pfJet.chargedHadronEnergyFraction()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("f_HE_c",1);
      // charged multiplicity>0
      if( (pfJet.chargedMultiplicity()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("M_c>0.8",1);
      // charged EM energy fraction < 0.8
      if( (pfJet.chargedEmEnergyFraction()) > 0.8 ) continue;
      h_PFJet_cutflow_->Fill("f_EM_c>0.8",1);
    }
    h_PFJet_PT_->Fill(pfJet.pt());
    h_PFJet_eta_->Fill(pfJet.eta());
    if(pfJet.pt()>=40){
      pfJet_HT += pfJet.pt();
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
  //

  for(auto & track:tracks){
    if( track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1) && track.hitPattern().pixelLayersWithMeasurement()>=2 && track.hitPattern().stripLayersWithMeasurement() >= 6 && (std::abs(track.dxy()/track.dxyError())) > 4)
      h_n_1_track_PT_->Fill(track.pt());
    if( track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1) && track.pt()>1 && track.hitPattern().stripLayersWithMeasurement() >= 6 && (std::abs(track.dxy()/track.dxyError())) > 4)
      h_n_1_track_npxl_->Fill(track.hitPattern().pixelLayersWithMeasurement());
    if( track.pt()>1 && track.hitPattern().pixelLayersWithMeasurement()>=2 && track.hitPattern().stripLayersWithMeasurement() >= 6 && (std::abs(track.dxy()/track.dxyError())) > 4)
      h_n_1_track_nmin_->Fill(0);
    if( track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1) && track.hitPattern().pixelLayersWithMeasurement()>=2 && track.pt()>1 && (std::abs(track.dxy()/track.dxyError())) > 4)
      h_n_1_track_nstl_->Fill(track.hitPattern().stripLayersWithMeasurement());
    if( track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1) && track.hitPattern().pixelLayersWithMeasurement()>=2 && track.pt()>1 && track.hitPattern().stripLayersWithMeasurement() >= 6)
      h_n_1_Nsig_->Fill((std::abs(track.dxy()/track.dxyError())));
    
  }

  /*
  for (const auto & track:tracks){
    h_Track_cutflow_->Fill("allTracks", 1);

    // track pt > 1 GeV
    if(track.pt()<=1) continue;
    h_Track_cutflow_->Fill("pt>1GeV", 1);

    uint32_t npxl = 0;
    uint32_t nstl = 0;
    uint32_t rmin = 999;
    std::map<uint32_t, std::set<uint32_t>> pxlLayers;
    std::map<uint32_t, std::set<uint32_t>> stlLayers;
    // hit pattern of the track
    const reco::HitPattern &p = track.hitPattern();

    //loop over the hits of the track
    for (int i=0; i<p.numberOfAllHits(reco::HitPattern::TRACK_HITS); ++i){
      uint32_t hit = p.getHitPattern(reco::HitPattern::TRACK_HITS, i);
      // check if hit is valid
      if(p.validHitFilter(hit)){
        uint32_t subStructure = p.getSubStructure(hit);
        uint32_t layer = p.getLayer(hit);
        // hits from pixel layers
        if(p.pixelHitFilter(hit)){
          // find rmin by only looking at hit in PXB (subStructure=1)
          if(subStructure==1){
            if(rmin>layer){
              rmin = layer;
            }
          }
          if(pxlLayers.find(subStructure)==pxlLayers.end()){
            pxlLayers[subStructure] = std::set<uint32_t>{layer};
          }
          else{
            pxlLayers[subStructure].insert(layer);
          }
        }
        // hits from strip layers
        else if(p.stripHitFilter(hit)){
          if(stlLayers.find(subStructure)==stlLayers.end()){
            stlLayers[subStructure] = std::set<uint32_t>{layer};
          }
          else{
            stlLayers[subStructure].insert(layer);
          }
        }
      }
    }
    for(std::map<uint32_t, std::set<uint32_t>>::iterator iter=pxlLayers.begin(); iter != pxlLayers.end(); ++iter){
      npxl += iter->second.size();
    }
    for(std::map<uint32_t, std::set<uint32_t>>::iterator iter=stlLayers.begin(); iter != stlLayers.end(); ++iter){
      nstl += iter->second.size();
    }

    // rmin = 1
    if(rmin!=1) {
      continue;
    }
    h_Track_cutflow_->Fill("rmin=1", 1);
  
    // npxl >= 2
    if(npxl < 2){
      continue;
    }
    h_Track_cutflow_->Fill("npxl>=2", 1);

    // nstl >= 6
    if(nstl < 6){
      continue;
    }
    h_Track_cutflow_->Fill("nstl>=6", 1);

    // |dxy|/sigma_dxy > 4 
    if( (std::abs(track.dxy())/track.dxyError()) <= 4 ){
      continue;
    }
    h_Track_cutflow_->Fill("dxy/sigma_dxy>4", 1);


    h_track_PT_->Fill(track.pt());
  }
  */


  //for(size_t j=0; j < genJetHandle_->size(); ++j){
  //  const auto& genJet = (*genJetHandle_)[j];
  //  if(genJet.pt()>=40){
  //    HT += genJet.pt();
  //  }
  //}
  //h_jet_HT_->Fill(HT);
  //h_jet_N_->Fill(genJetHandle_->size());
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
DemoAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
  h_PFJet_PT_ = fileService->make<TH1F>("h_PFJet_PT","h_PFJet_PT",100,0,100);
  h_PFJet_eta_ = fileService->make<TH1F>("h_PFJet_eta","h_PFJet_eta",100,-10,10);
  h_PFJet_HT_ = fileService->make<TH1F>("h_PFJet_HT","h_PFJet_HT",50,0,5000);
  h_PFJet_N_ = fileService->make<TH1F>("h_PFJet_N","h_PFJet_N",60,0,60);
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
}

void DemoAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
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

void
DemoAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("genJetTag",edm::InputTag("ak4GenJets"));
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
