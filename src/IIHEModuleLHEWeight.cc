#include "UserCode/IIHETree/interface/IIHEModuleLHEWeight.h"
#include <iostream>
#include <TMath.h>
#include <vector>
#include <typeinfo>
#include <sstream>
#include <string>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"


using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleLHEWeight::IIHEModuleLHEWeight(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
	lheEventLabel_ = iC.consumes<LHEEventProduct> (iConfig.getParameter<InputTag>("LHELabel"));
        addSysLHE_       = iConfig.getUntrackedParameter<bool>("addSysLHE_"     ) ;
}
IIHEModuleLHEWeight::~IIHEModuleLHEWeight(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleLHEWeight::beginJob(){
	setBranchType(kVectorFloat) ;
	addBranch("LHE_Pt");
	addBranch("LHE_Eta");
	addBranch("LHE_Phi");
	addBranch("LHE_E");
	setBranchType(kVectorInt) ;
	addBranch("LHE_pdgid");
	addBranch("LHE_status");
        addBranch("LHE_motherIndex1");
        addBranch("LHE_motherIndex2");

	setBranchType(kFloat) ;
	addBranch("LHE_weight_nominal");
	setBranchType(kVectorFloat) ;
	addBranch("LHE_weight_sys");

}

// ------------ method called to for each event  ------------
void IIHEModuleLHEWeight::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	edm::Handle<LHEEventProduct> lhe_handle;
	iEvent.getByToken(lheEventLabel_, lhe_handle);


	if (lhe_handle.isValid()){

		std::vector<lhef::HEPEUP::FiveVector> lheParticles = lhe_handle->hepeup().PUP;
		ROOT::Math::PxPyPzEVector cand_;
		for (unsigned i = 0; i < lheParticles.size(); ++i) {
			cand_ = ROOT::Math::PxPyPzEVector(lheParticles[i][0],lheParticles[i][1],lheParticles[i][2],lheParticles[i][3]);
			store("LHE_Pt",(cand_).Pt());
			store("LHE_Eta",(cand_).Eta());
			store("LHE_Phi",(cand_).Phi());
			store("LHE_E",(cand_).E());
			store("LHE_pdgid",lhe_handle->hepeup().IDUP[i]);
			store("LHE_status",lhe_handle->hepeup().ISTUP[i]);
                        store("LHE_motherIndex1",lhe_handle->hepeup().MOTHUP[i].first);
                        store("LHE_motherIndex2",lhe_handle->hepeup().MOTHUP[i].second);
		}

		if(!(lhe_handle->weights().empty())) {
			store("LHE_weight_nominal",(float) lhe_handle->weights().at(0).wgt);
                        if(addSysLHE_){
	  	        	for (unsigned i = 0; i < lhe_handle->weights().size(); ++i) {
	  	          		string target(lhe_handle->weights().at(i).id.data());
	  		 		if(sumofLHEWeights_.size() != lhe_handle->weights().size()) {
	  					sumofLHEWeights_.push_back(0);
	  					LHEweightsId_.push_back(target);
	  				}
	  				store("LHE_weight_sys",(float) lhe_handle->weights().at(i).wgt);
	  				sumofLHEWeights_[i] += (float) lhe_handle->weights().at(i).wgt;     
	  			}
			}
		} else {
			store("LHE_weight_nominal",1);
			store("LHE_weight_sys", 1);

		}
	}
}

void IIHEModuleLHEWeight::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleLHEWeight::beginEvent(){}
void IIHEModuleLHEWeight::endEvent(){}
void IIHEModuleLHEWeight::endJob(){
	addFVValueToMetaTree("mc_sumofLHEWeights", sumofLHEWeights_) ;
	addCVValueToMetaTree("mc_LHEweightsId",LHEweightsId_) ;
}

DEFINE_FWK_MODULE(IIHEModuleLHEWeight);
