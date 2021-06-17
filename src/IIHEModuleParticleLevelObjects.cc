#include "UserCode/IIHETree/interface/IIHEModuleParticleLevelObjects.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleParticleLevelObjects::IIHEModuleParticleLevelObjects(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  particleLevelJetsToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelJetsCollection"));
  particleLevelPhotonToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelPhotonCollection"));
  particleLevelak1DressedLeptonToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelak1DressedLeptonCollection"));
  particleLevelMETToken_ = iC.consumes<std::vector<reco::MET>>(iConfig.getParameter<edm::InputTag>("particleLevelMETCollection"));
}
IIHEModuleParticleLevelObjects::~IIHEModuleParticleLevelObjects(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleParticleLevelObjects::beginJob(){
  setBranchType(kVectorFloat) ;
  addBranch("pl_jet_pt"   ) ;
  addBranch("pl_jet_eta"   ) ;
  addBranch("pl_jet_phi"   ) ;
  addBranch("pl_lep_pt"   ) ;
  addBranch("pl_lep_eta"   ) ;
  addBranch("pl_lep_phi"   ) ;
  addBranch("pl_ph_pt"   ) ;
  addBranch("pl_ph_eta"   ) ;
  addBranch("pl_ph_phi"   ) ;
  addBranch("pl_MET_pt"   ) ;
  addBranch("pl_MET_phi"   ) ;

  setBranchType(kVectorInt) ;
  addBranch("pl_lep_pdgid"   ) ;
  addBranch("pl_lep_charge"   ) ;
  addBranch("pl_jet_pdgid"   ) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleParticleLevelObjects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<std::vector<reco::GenJet>> particleLevelJetsHandle_;
  iEvent.getByToken(particleLevelJetsToken_, particleLevelJetsHandle_);

  edm::Handle<std::vector<reco::GenJet>> particleLevelak1DressedLeptonHandle_;
  iEvent.getByToken(particleLevelak1DressedLeptonToken_, particleLevelak1DressedLeptonHandle_);

  edm::Handle<std::vector<reco::MET>> particleLevelMETHandle_;
  iEvent.getByToken(particleLevelMETToken_, particleLevelMETHandle_);

  edm::Handle<std::vector<reco::GenJet>> particleLevelPhotonHandle_;
  iEvent.getByToken(particleLevelJetsToken_, particleLevelPhotonHandle_);

  for (size_t j = 0; j < particleLevelJetsHandle_->size();++j){
    if(particleLevelJetsHandle_->at(j).pt()<10.0) continue;
    store("pl_jet_pt"     , particleLevelJetsHandle_->at(j).pt()) ;
    store("pl_jet_eta"    , particleLevelJetsHandle_->at(j).eta()) ;
    store("pl_jet_phi"    , particleLevelJetsHandle_->at(j).phi()) ;
    store("pl_jet_pdgid"    , particleLevelJetsHandle_->at(j).pdgId()) ;
  }


  for (size_t j = 0; j < particleLevelak1DressedLeptonHandle_->size();++j){
    if(particleLevelak1DressedLeptonHandle_->at(j).pt()<10.0) continue;
    store("pl_lep_pt"     , particleLevelak1DressedLeptonHandle_->at(j).pt()) ;
    store("pl_lep_eta"    , particleLevelak1DressedLeptonHandle_->at(j).eta()) ;
    store("pl_lep_phi"    , particleLevelak1DressedLeptonHandle_->at(j).phi()) ;
    store("pl_lep_charge"    , particleLevelak1DressedLeptonHandle_->at(j).charge()) ;
    store("pl_lep_pdgid"    , particleLevelak1DressedLeptonHandle_->at(j).pdgId()) ;
  }



  for (size_t j = 0; j < particleLevelPhotonHandle_->size();++j){
    if(particleLevelPhotonHandle_->at(j).pt()<10.0) continue;
    store("pl_ph_pt"        , particleLevelPhotonHandle_->at(j).pt()) ;
    store("pl_ph_eta"       , particleLevelPhotonHandle_->at(j).eta()) ;
    store("pl_ph_phi"       , particleLevelPhotonHandle_->at(j).phi()) ;
  }

    store("pl_MET_pt"     , particleLevelMETHandle_->at(0).pt()) ;
    store("pl_MET_phi"    , particleLevelMETHandle_->at(0).phi()) ;

}
void IIHEModuleParticleLevelObjects::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleParticleLevelObjects::beginEvent(){}
void IIHEModuleParticleLevelObjects::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleParticleLevelObjects::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleParticleLevelObjects);
