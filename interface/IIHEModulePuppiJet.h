#ifndef UserCode_IIHETree_IIHEModulePuppiJet_h
#define UserCode_IIHETree_IIHEModulePuppiJet_h
#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
// class decleration

class IIHEModulePuppiJet : public IIHEModule {
public:
  explicit IIHEModulePuppiJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModulePuppiJet(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModulePuppiJet();

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
  edm::InputTag pfPuppiJetLabel_                  ;
  edm::InputTag pfPuppiJetLabelPrecor_            ;
  edm::InputTag pfPuppiJetLabelSmeared_           ;
  edm::InputTag pfPuppiJetLabelSmearedJetResUp_   ;
  edm::InputTag pfPuppiJetLabelSmearedJetResDown_ ;

  edm::EDGetTokenT<edm::View<pat::Jet> > pfPuppiJetToken_                  ;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfPuppiJetTokenPrecor_            ;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfPuppiJetTokenSmeared_           ;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfPuppiJetTokenSmearedJetResUp_   ;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfPuppiJetTokenSmearedJetResDown_ ;


  float ETThreshold_ ;
  bool  isMC_        ;
};
#endif
