#include "UserCode/IIHETree/interface/IIHEModuleMET.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

//////////////////////////////////////////////////////////////////////////////////////////
//                             IIHEMETVariable classes                            
//////////////////////////////////////////////////////////////////////////////////////////
IIHEMETVariableBase::IIHEMETVariableBase(std::string prefix, std::string name, int type){
  name_       = name ;
  branchName_ = prefix + "_" + name_ ;
  branchType_ = type ;
}
bool IIHEMETVariableBase::addBranch(IIHEAnalysis* analysis){
  return analysis->addBranch(branchName_, branchType_) ;
}

IIHEMETVariableInt::IIHEMETVariableInt(std::string prefix, std::string name):
IIHEMETVariableBase(prefix, name, kVectorInt){
  reset() ;
}
void IIHEMETVariableInt::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_) ;
}

IIHEMETVariableFloat::IIHEMETVariableFloat(std::string prefix, std::string name):
IIHEMETVariableBase(prefix, name, kFloat){
  reset() ;
}
void IIHEMETVariableFloat::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_ ) ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                  IIHEMET class                                 
//////////////////////////////////////////////////////////////////////////////////////////
IIHEMETWrapper::IIHEMETWrapper(std::string prefix){
  prefix_ = prefix ;

  Pt_               = new IIHEMETVariableFloat  (prefix_, "Pt"                 ) ;
  Px_               = new IIHEMETVariableFloat  (prefix_, "Px"                 ) ;
  Py_               = new IIHEMETVariableFloat  (prefix_, "Py"                 ) ;
  phi_              = new IIHEMETVariableFloat  (prefix_, "phi"                ) ;
  significance_     = new IIHEMETVariableFloat  (prefix_, "significance"       ) ;

  variables_.push_back((IIHEMETVariableBase*) Pt_                   ) ;
  variables_.push_back((IIHEMETVariableBase*) Px_                   ) ;
  variables_.push_back((IIHEMETVariableBase*) Py_                   ) ;
  variables_.push_back((IIHEMETVariableBase*) phi_                  ) ;
  variables_.push_back((IIHEMETVariableBase*) significance_         ) ;
}

void IIHEMETWrapper::addBranches(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->addBranch(analysis) ;
  }
}
void IIHEMETWrapper::reset(){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->reset() ;
  }
}

void IIHEMETWrapper::fill(pat::MET MET){
  Pt_        ->fill(MET.pt()                ) ;
  Px_        ->fill(MET.px()                ) ;
  Py_        ->fill(MET.py()                ) ;
  phi_       ->fill(MET.phi()               ) ;
  significance_   ->fill(MET.metSignificance()   ) ;
}
void IIHEMETWrapper::store(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->store(analysis) ;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                  Main IIHEMETModule                                 
//////////////////////////////////////////////////////////////////////////////////////////


IIHEModuleMET::IIHEModuleMET(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig),
  metnominalWrapper_(new IIHEMETWrapper("MET_nominal"))
//  metWrapper_(new IIHEMETWrapper("MET")),
//  metT1Wrapper_(new IIHEMETWrapper("MET_T1")),
//  metT1JetEnDownWrapper_(new IIHEMETWrapper("MET_T1JetEnDown")),
//  metT1JetEnUpWrapper_(new IIHEMETWrapper("MET_T1JetEnUp")),
//  metT1SmearWrapper_(new IIHEMETWrapper("MET_T1Smear")),
//  metT1SmearJetEnDownWrapper_(new IIHEMETWrapper("MET_T1SmearJetEnDown")),
//  metT1SmearJetEnUpWrapper_(new IIHEMETWrapper("MET_T1SmearJetEnUp")),
//  metT1SmearJetResDownWrapper_(new IIHEMETWrapper("MET_T1SmearJetResDown")),
//  metT1SmearJetResUpWrapper_(new IIHEMETWrapper("MET_T1SmearJetResUp")),
//  metT1TxyWrapper_(new IIHEMETWrapper("MET_T1Txy")),
//  metFinalWrapper_(new IIHEMETWrapper("MET_FinalCollection"))
{
  pfMETToken_                               =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("METCollection"));
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC") ;
}
IIHEModuleMET::~IIHEModuleMET(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMET::beginJob(){
  IIHEAnalysis* analysis = parent_ ;
  metnominalWrapper_->addBranches(analysis) ;
  if(isMC_){
    setBranchType(kFloat) ;
    addBranch("MET_gen_pt"   ) ;
    addBranch("MET_gen_phi"   ) ;
    setBranchType(kVectorFloat) ;
    addBranch("MET_Type1Unc_Px") ;
    addBranch("MET_Type1Unc_Py") ;
    addBranch("MET_Type1Unc_Pt") ;

  }
}

// ------------ method called to for each event  ------------
void IIHEModuleMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::MET> > pfMETHandle_;
  iEvent.getByToken(pfMETToken_, pfMETHandle_);



  metnominalWrapper_->reset() ;

  IIHEAnalysis* analysis = parent_ ;  
  Ptr<pat::MET> pfMET = pfMETHandle_->ptrAt( 0 );

  metnominalWrapper_->fill(pfMETHandle_->front()) ;
  metnominalWrapper_->store(analysis) ;


  if (isMC_){

    store("MET_gen_pt"    , pfMET->genMET()->pt()     ) ;
    store("MET_gen_phi"   , pfMET->genMET()->phi()     ) ;

  }
}
void IIHEModuleMET::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMET::beginEvent(){}
void IIHEModuleMET::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMET::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMET);
