#ifndef UserCode_IIHETree_IIHEModuleJet_h
#define UserCode_IIHETree_IIHEModuleJet_h
#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
// class decleration

class IIHEModuleJet : public IIHEModule {
public:
  explicit IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModuleJet(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModuleJet();
  
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
  edm::InputTag pfJetLabel_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetToken_;
  edm::InputTag pfJetLabelSmeared_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmeared_;
  edm::InputTag pfJetLabelSmearedJetResUp_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmearedJetResUp_;
  edm::InputTag pfJetLabelSmearedJetResDown_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmearedJetResDown_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenPrecor_;

  float ETThreshold_ ;
  bool isMC_;
};
#endif
