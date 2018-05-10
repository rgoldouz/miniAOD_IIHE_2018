#include "UserCode/IIHETree/interface/IIHEModuleL1.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleL1::IIHEModuleL1(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  l1noniso_ = iC.consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("l1NonIsoCollection"));
  glbalgblk_token_ =  iC.consumes<BXVector<GlobalAlgBlk>>(edm::InputTag("gtStage2Digis"));
}
IIHEModuleL1::~IIHEModuleL1(){}

void IIHEModuleL1::beginJob(){

  setBranchType(kVectorFloat);
  addBranch("L1_EG_pt");
  addBranch("L1_EG_eta");
  addBranch("L1_EG_phi");
  setBranchType(kVectorInt);
  addBranch("L1_EG_Iso");
  addBranch("L1_pass_final");
}

void IIHEModuleL1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<BXVector<l1t::EGamma>> egammas;
  iEvent.getByToken( l1noniso_, egammas );
  for(int i = egammas->getFirstBX(); i <= egammas->getLastBX(); ++i) {
    for(std::vector<l1t::EGamma>::const_iterator eg = egammas->begin(i); eg != egammas->end(i); ++eg) {
      if (i != 0) continue;
        store("L1_EG_pt" , eg->pt());
        store("L1_EG_eta", eg->eta());
        store("L1_EG_phi", eg->phi());
        store("L1_EG_Iso", eg->hwIso());
    }
  }

  edm::Handle<BXVector<GlobalAlgBlk>> l1algos;
  iEvent.getByToken(glbalgblk_token_, l1algos);

  for(int i =0; i <512; i++){
    store("L1_pass_final", l1algos->at(0,0).getAlgoDecisionFinal(i));
  }
}
void IIHEModuleL1::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleL1::beginEvent(){}
void IIHEModuleL1::endEvent(){}


void IIHEModuleL1::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleL1);



