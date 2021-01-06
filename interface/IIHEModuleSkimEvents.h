#ifndef UserCode_IIHETree_IIHEModuleSkimEvents_h
#define UserCode_IIHETree_IIHEModuleSkimEvents_h

#include "UserCode/IIHETree/interface/IIHEModule.h"

// class decleration
class IIHEModuleSkimEvents : public IIHEModule {
private:
  int nAcceptAll_  = 0 ;
  
  double ptThresholdEle_;
  double ptThresholdmu_;
  double ptThresholdTau_;
  double ptThresholdPh_;
  double ptThresholdFatJet_;

  bool EleAccept_  ;
  bool EleMuAccept_;
  bool EleTauAccept_;
  bool MuAccept_;
  bool MuTauAccept_;
  bool TauAccept_;
  bool PhAccept_;
  bool FatJetAccept_;
  bool FatJetMuAccept_;
  bool FatJetPhAccept_;

  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonCollectionToken_;
  edm::EDGetTokenT<edm::View<pat::Tau> > tauCollectionToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollectionToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> >  FatJetToken_;

  edm::InputTag      electronCollectionLabel_ ;
  edm::InputTag      muonCollectionLabel_ ;
  edm::InputTag      tauCollectionLabel_ ;
  edm::InputTag      photonCollectionLabel_ ;
  edm::InputTag      FatJetLabel_;

public:
  explicit IIHEModuleSkimEvents(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModuleSkimEvents(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModuleSkimEvents();
  
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
};
#endif
