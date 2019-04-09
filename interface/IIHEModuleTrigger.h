#ifndef UserCode_IIHETree_IIHEModuleTrigger_h
#define UserCode_IIHETree_IIHEModuleTrigger_h

#include "UserCode/IIHETree/interface/IIHEModule.h"

#include "UserCode/IIHETree/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// class decleration
class IIHEModuleTrigger : public IIHEModule {
public:
  explicit IIHEModuleTrigger(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModuleTrigger(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModuleTrigger();
  
  void   pubBeginJob(){   beginJob() ; } ;
  void pubBeginEvent(){ beginEvent() ; } ;
  void   pubEndEvent(){   endEvent() ; } ;
  virtual void pubAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ analyze(iEvent, iSetup) ; } ;
  
  virtual void beginEvent() ;
  virtual void endEvent() ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
private:
  int addBranches() ;
  
  bool addHLTrigger(HLTrigger*) ;
  std::vector<L1Trigger*> L1Triggers_ ;
  std::vector<HLTrigger*> HLTriggers_ ;
  std::vector<HLTrigger*> HLTriggersPAT_ ;
  std::vector<HLTrigger*> HLTriggersRECO_ ;
  bool changed_ = true ;

  HLTConfigProvider hltConfig_ ;
  HLTConfigProvider hltConfigPAT_ ;
  HLTConfigProvider hltConfigRECO_ ;
  edm::InputTag hlTriggerResultsTag_ ;
  std::vector<std::string> HLTNamesFromConfig_ ;
  std::vector<std::string> HLTNamesFromConfigPAT_ ;
  std::vector<std::string> HLTNamesFromConfigRECO_ ;
  std::vector<std::string> triggerNamesFromPSet_ ;
  std::vector<std::string> savedHLTriggers_ ; 
  void clearHLTrigger(){
  for (unsigned int k=0;k<HLTriggers_.size();++k){
    delete HLTriggers_[k];
  }
  HLTriggers_.clear();} ;


  edm::InputTag  triggerBitsLabel_;
  edm::InputTag  triggerObjectsLabel_;
  edm::InputTag  triggerPrescalesLabel_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;


  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsTokenRECO_;
  edm::InputTag triggerResultsLabel_;
  edm::InputTag triggerResultsLabelRECO_;

 
  bool isSingleElectonTriggerName(std::string) ;
  bool isDoubleElectonTriggerName(std::string) ;
  bool isTripleElectonTriggerName(std::string) ;
  bool    isSingleMuonTriggerName(std::string) ;
  bool    isDoubleMuonTriggerName(std::string) ;
  bool    isTripleMuonTriggerName(std::string) ;

  bool isSingleElectronSingleMuonTriggerName(std::string) ;
  bool isDoubleElectronSingleMuonTriggerName(std::string) ;
  bool isSingleElectronDoubleMuonTriggerName(std::string) ;
  
  int nEvents_ ;
  int nWasRun_ ;
  int nAccept_ ;
  int nErrors_ ;
  
  bool includeSingleElectronTriggers_ ;
  bool includeDoubleElectronTriggers_ ;
  bool includeTripleElectronTriggers_ ;
  bool includeSingleMuonTriggers_   ;
  bool includeDoubleMuonTriggers_   ;
  bool includeTripleMuonTriggers_   ;
  bool includeSingleElectronSingleMuonTriggers_ ;
  bool includeSingleElectronDoubleMuonTriggers_ ;
  bool includeDoubleElectronSingleMuonTriggers_ ;
  bool includeSinglePhotonTriggers_ ;
  bool includeSingleTauTriggers_ ;
  bool includeMETTriggers_ ;
};
#endif
