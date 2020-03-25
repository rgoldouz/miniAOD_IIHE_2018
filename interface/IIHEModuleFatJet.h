#ifndef UserCode_IIHETree_IIHEModuleFatJet_h
#define UserCode_IIHETree_IIHEModuleFatJet_h
#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
// class decleration

class IIHEModuleFatJet : public IIHEModule {
public:
  explicit IIHEModuleFatJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModuleFatJet(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModuleFatJet();

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

  edm::InputTag   FatJetLabel_;
  edm::InputTag   FatJetSmearedLabel_          ;
  edm::InputTag   FatJetSmearedJetResUpLabel_  ;
  edm::InputTag   FatJetSmearedJetResDownLabel_;

  edm::EDGetTokenT<edm::View<pat::Jet> >  FatJetToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> >  FatJetSmearedToken_          ;
  edm::EDGetTokenT<edm::View<pat::Jet> >  FatJetSmearedJetResUpToken_  ;
  edm::EDGetTokenT<edm::View<pat::Jet> >  FatJetSmearedJetResDownToken_;

  float ETThreshold_;
  bool  isMC_ ;
};
#endif
