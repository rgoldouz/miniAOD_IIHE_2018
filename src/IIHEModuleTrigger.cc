#include "UserCode/IIHETree/interface/IIHEModuleTrigger.h"
#include "UserCode/IIHETree/interface/TriggerObject.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleTrigger::IIHEModuleTrigger(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
//  hlTriggerResultsTag_ = iConfig.getParameter<edm::InputTag>("TriggerResults") ;
  nEvents_ = 0 ;
  nWasRun_ = 0 ;
  nAccept_ = 0 ;
  nErrors_ = 0 ;

  triggerBitsLabel_       = iConfig.getParameter<edm::InputTag>("triggerResultsCollectionHLT") ;
  triggerObjectsLabel_    = iConfig.getParameter<edm::InputTag>("triggerObjectStandAloneCollection") ;
  triggerPrescalesLabel_  = iConfig.getParameter<edm::InputTag>("patTriggerCollection") ;

  triggerBits_ = iC.consumes<edm::TriggerResults>(InputTag(triggerBitsLabel_));
  triggerObjects_ = iC.consumes<pat::TriggerObjectStandAloneCollection>(InputTag(triggerObjectsLabel_));
  triggerPrescales_ = iC.consumes<pat::PackedTriggerPrescales>(InputTag(triggerPrescalesLabel_));

  triggerResultsLabel_                      = iConfig.getParameter<edm::InputTag>("triggerResultsCollectionPAT") ;
  triggerResultsToken_                      = iC.consumes<edm::TriggerResults>(triggerResultsLabel_);

  triggerResultsLabelRECO_                      = iConfig.getParameter<edm::InputTag>("triggerResultsCollectionRECO") ;
  triggerResultsTokenRECO_                      = iC.consumes<edm::TriggerResults>(triggerResultsLabelRECO_);
  
  std::string triggersIn = iConfig.getUntrackedParameter<std::string>("triggers") ;
  triggerNamesFromPSet_ = splitString(triggersIn, ",") ;
 
  std::cout << triggersIn << std::endl ;
  includeSingleElectronTriggers_ = (triggersIn.find("singleElectron")!=std::string::npos) ;
  includeDoubleElectronTriggers_ = (triggersIn.find("doubleElectron")!=std::string::npos) ;
  includeTripleElectronTriggers_ = (triggersIn.find("tripleElectron")!=std::string::npos) ;
  includeSingleMuonTriggers_     = (triggersIn.find("singleMuon"    )!=std::string::npos) ;
  includeDoubleMuonTriggers_     = (triggersIn.find("doubleMuon"    )!=std::string::npos) ;
  includeTripleMuonTriggers_     = (triggersIn.find("tripleMuon"    )!=std::string::npos) ;
  includeSingleElectronSingleMuonTriggers_ = (triggersIn.find("singleElectronSingleMuon")!=std::string::npos) ;
  includeSingleElectronDoubleMuonTriggers_ = (triggersIn.find("singleElectronDoubleMuon")!=std::string::npos) ;
  includeDoubleElectronSingleMuonTriggers_ = (triggersIn.find("doubleElectronSingleMuon")!=std::string::npos) ;
  includeSinglePhotonTriggers_ = (triggersIn.find("singlePhoton" )!=std::string::npos) ; 
  includeMETTriggers_ = (triggersIn.find("MET" )!=std::string::npos) ; 
  includeSingleTauTriggers_ = (triggersIn.find("singleTau" )!=std::string::npos) ;

  std::cout << "Including single electron triggers:            " << includeSingleElectronTriggers_ << std::endl ;
  std::cout << "Including double electron triggers:            " << includeDoubleElectronTriggers_ << std::endl ;
  std::cout << "Including triple electron triggers:            " << includeTripleElectronTriggers_ << std::endl ;
  std::cout << "Including single muon triggers:                " << includeSingleMuonTriggers_     << std::endl ;
  std::cout << "Including double muon triggers:                " << includeDoubleMuonTriggers_     << std::endl ;
  std::cout << "Including triple muon triggers:                " << includeTripleMuonTriggers_     << std::endl ;
  std::cout << "Including single electron single muon triggers:" << includeSingleElectronSingleMuonTriggers_ << std::endl ;
  std::cout << "Including single electron double muon triggers:" << includeSingleElectronDoubleMuonTriggers_ << std::endl ;
  std::cout << "Including double electron single muon triggers:" << includeDoubleElectronSingleMuonTriggers_ << std::endl ;
  std::cout << "Including single photon triggers:              " << includeSinglePhotonTriggers_ << std::endl ; 
  std::cout << "Including MET triggers:                        " << includeMETTriggers_ << std::endl ;
  std::cout << "Including single tau triggers:                 " << includeSingleTauTriggers_ << std::endl ;
}
IIHEModuleTrigger::~IIHEModuleTrigger(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTrigger::beginJob(){
}
//remove some HLT paths 
bool IIHEModuleTrigger::addHLTrigger(HLTrigger* hlt){
  if(   hlt->nSubstringInString(hlt->name(), "Upsilon")
     || hlt->nSubstringInString(hlt->name(), "LowMass")
     || hlt->nSubstringInString(hlt->name(), "Tau")
     || hlt->nSubstringInString(hlt->name(), "MET")
     || hlt->nSubstringInString(hlt->name(), "HT")
     || hlt->nSubstringInString(hlt->name(), "Jpsi")
     || hlt->nSubstringInString(hlt->name(), "PsiPrime")
     || hlt->nSubstringInString(hlt->name(), "Scouting")
     || hlt->nSubstringInString(hlt->name(), "Jet")
     || hlt->nSubstringInString(hlt->name(), "IP")
     || hlt->nSubstringInString(hlt->name(), "PF")
  ) {
      delete hlt;
      return false ;
  }

  for(unsigned int i=0 ; i<HLTriggers_.size() ; ++i){
    if(HLTriggers_.at(i)->name()==hlt->name()){
      return false ;
    }
  }
  if(   hlt->nSubstringInString(hlt->name(), "Ele35_WPTight") 
     || hlt->nSubstringInString(hlt->name(), "DoubleEle33_CaloIdL_MW_v") 
     || hlt->nSubstringInString(hlt->name(), "DoubleEle25_CaloIdL_MW_v")
     || hlt->nSubstringInString(hlt->name(), "Ele23_Ele12" )
     || hlt->nSubstringInString(hlt->name(), "Ele32_WPTight" )
     || hlt->nSubstringInString(hlt->name(), "DiEle27_WPTightCaloOnly" )
     || hlt->nSubstringInString(hlt->name(), "Ele27_eta2p1_WPTight" ) 
  ) {
    hlt->saveFilters();
  }
//  hlt->savePrescale();
  HLTriggers_.push_back(hlt) ;
  return true ;
}

int IIHEModuleTrigger::addBranches(){
// Just loop over the HLTrigger objects and make branches for each one, including its
// filters.  Return the number of branches we've added so that we can get an idea of
//the trigger overhead.  (This is not used in the main analysis, but can be trivially
// printed to screen.)
  int result = 0 ;
  IIHEAnalysis* analysis = parent_ ;
  for(unsigned int i=0 ; i<HLTriggers_.size() ; i++){
    if (std::find(savedHLTriggers_.begin(), savedHLTriggers_.end(), HLTriggers_.at(i)->name().substr(0, HLTriggers_.at(i)->name().find("_v"))) != savedHLTriggers_.end()) continue;
    savedHLTriggers_.push_back(HLTriggers_.at(i)->name().substr(0, HLTriggers_.at(i)->name().find("_v"))); 
    result += HLTriggers_.at(i)->createBranches(analysis) ;
  }
  return result ;
}


// ------------ method called to for each event  ------------
void IIHEModuleTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Trigger information
  //MET filter 
  // get hold of TriggerResults
  edm::Handle<TriggerResults> HLTR ;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, HLTR) ;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<TriggerResults> triggerResultsCollection_ ;
  iEvent.getByToken(triggerResultsToken_, triggerResultsCollection_);

  edm::Handle<TriggerResults> triggerResultsCollectionRECO_ ;
  iEvent.getByToken(triggerResultsTokenRECO_, triggerResultsCollectionRECO_);

  // Now fill the values
  IIHEAnalysis* analysis = parent_ ;
  for(unsigned int i=0 ; i<HLTriggersPAT_.size() ; i++){
    HLTrigger* hlt = HLTriggersPAT_.at(i) ;
    hlt->status(triggerResultsCollection_) ;
    hlt->store(analysis) ;
  } 

  for(unsigned int i=0 ; i<HLTriggersRECO_.size() ; i++){
    HLTrigger* hlt = HLTriggersRECO_.at(i) ;
    hlt->status(triggerResultsCollectionRECO_) ;
    hlt->store(analysis) ;
  }

  for(unsigned int i=0 ; i<HLTriggers_.size() ; i++){
    HLTrigger* hlt = HLTriggers_.at(i) ;
    hlt->fullStatus(iEvent, iSetup, hltConfig_, HLTR, triggerObjects, triggerPrescales ,analysis) ;
    hlt->store(analysis) ;
  }
  nEvents_++ ;
}

void IIHEModuleTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
  IIHEAnalysis* analysis = parent_ ;
  if(changed_){
    hltConfigPAT_.init(iRun, iSetup, triggerResultsLabel_.process(), changed_);
    HLTNamesFromConfigPAT_ = hltConfigPAT_.triggerNames() ;
    for(unsigned int i=0 ; i<HLTNamesFromConfigPAT_.size() ; ++i){
      std::string namePAT = HLTNamesFromConfigPAT_.at(i) ;
      HLTrigger* hltPAT = new HLTrigger(namePAT, hltConfigPAT_) ;
      HLTriggersPAT_.push_back(hltPAT) ;
      hltPAT->createBranches(analysis) ;
    }

    hltConfigRECO_.init(iRun, iSetup, triggerResultsLabelRECO_.process(), changed_);
    HLTNamesFromConfigRECO_ = hltConfigRECO_.triggerNames() ;
    for(unsigned int i=0 ; i<HLTNamesFromConfigRECO_.size() ; ++i){
      std::string nameRECO = HLTNamesFromConfigRECO_.at(i) ;
      HLTrigger* hltRECO = new HLTrigger(nameRECO, hltConfigRECO_) ;
      HLTriggersRECO_.push_back(hltRECO) ;
      hltRECO->createBranches(analysis) ;
    }

  parent_->configureBranches() ;
//  changed_ = false ;
  }

  bool changed = true ;
  if(hltConfig_.init(iRun, iSetup, triggerBitsLabel_.process(), changed)){
    if(changed_){
      clearHLTrigger();
      if(false) hltConfig_.dump("Modules") ;
      
      // Get the updated list of trigger names
      HLTNamesFromConfig_ = hltConfig_.triggerNames() ;
      for(unsigned int i=0 ; i<HLTNamesFromConfig_.size() ; ++i){
        std::string name = HLTNamesFromConfig_.at(i) ;
        // Attempt to add the trigger
        bool addThisTrigger = false ;
        
        HLTrigger* hlt = new HLTrigger(name, hltConfig_) ;
        
        // First check to see if it's in the list of requested triggers
        if(hlt->isOnlySingleElectron()           && includeSingleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlyDoubleElectron()           && includeDoubleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlyTripleElectron()           && includeTripleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlySingleMuon()               && includeSingleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlyDoubleMuon()               && includeDoubleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlyTripleMuon()               && includeTripleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlySingleElectronSingleMuon() && includeSingleElectronSingleMuonTriggers_) addThisTrigger = true ;
        if(hlt->isOnlySingleElectronDoubleMuon() && includeSingleElectronDoubleMuonTriggers_) addThisTrigger = true ;
        if(hlt->isOnlyDoubleElectronSingleMuon() && includeDoubleElectronSingleMuonTriggers_) addThisTrigger = true ;
        if(hlt->isSinglePhoton() && includeSinglePhotonTriggers_ ) addThisTrigger = true ;       
        if(hlt->isSingleTau() && includeSingleTauTriggers_ ) addThisTrigger = true ; 
        if(hlt->isMET() && includeMETTriggers_) addThisTrigger = true ;
        // Only loop over trigger names if we have to
        if(addThisTrigger==false){
          for(unsigned int j=0 ; j<triggerNamesFromPSet_.size() ; ++j){
            if(triggerNamesFromPSet_.at(j)==name){
              addThisTrigger = true ;
              break ;
            }
          }
        }
        
        if(addThisTrigger==false) {
          delete hlt;
          continue ;}
//        hlt->savePrescale();
        addHLTrigger(hlt) ;
      }
      
      // Now we need to re-map the indices to the names, given that some new triggers may have been inserted to the menu
      for(unsigned int i=0 ; i<HLTriggers_.size() ; ++i){
        HLTriggers_.at(i)->beginRun(hltConfig_) ;
      }
      
      // Attempt to add branches
      addBranches() ;
      parent_->configureBranches() ;
      // Now reset things to 0
      nEvents_ = 0 ;
      nWasRun_ = 0 ;
      nAccept_ = 0 ;
      nErrors_ = 0 ;
      changed_ = false ;
    }
  }
  else{
    std::cout << "Failed to init hltConfig" << std::endl ;
  }

}

void IIHEModuleTrigger::beginEvent(){}
void IIHEModuleTrigger::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTrigger::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTrigger);
