#include "UserCode/IIHETree/interface/IIHEModuleSkimEvents.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleSkimEvents::IIHEModuleSkimEvents(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  ptThresholdEle_      = iConfig.getUntrackedParameter<double>("electronPtThreshold") ;
  ptThresholdmu_       = iConfig.getUntrackedParameter<double>("muonPtThreshold") ;
  ptThresholdTau_      = iConfig.getUntrackedParameter<double>("tauPtTThreshold") ;
  ptThresholdPh_       = iConfig.getUntrackedParameter<double>("photonPtThreshold") ;
  ptThresholdFatJet_       = iConfig.getUntrackedParameter<double>("fatjetPtThreshold") ;

  EleAccept_       = iConfig.getUntrackedParameter<bool>("skimEvents_Ele"     ) ;
  EleMuAccept_     = iConfig.getUntrackedParameter<bool>("skimEvents_EleMu"   ) ;
  EleTauAccept_    = iConfig.getUntrackedParameter<bool>("skimEvents_EleTau"  ) ;
  MuAccept_        = iConfig.getUntrackedParameter<bool>("skimEvents_Mu"      ) ;
  MuTauAccept_     = iConfig.getUntrackedParameter<bool>("skimEvents_MuTau"   ) ;
  TauAccept_       = iConfig.getUntrackedParameter<bool>("skimEvents_Tau"     ) ;
  PhAccept_        = iConfig.getUntrackedParameter<bool>("skimEvents_Ph"      ) ;
  FatJetAccept_    = iConfig.getUntrackedParameter<bool>("skimEvents_FatJet"      ) ;
  FatJetMuAccept_  = iConfig.getUntrackedParameter<bool>("skimEvents_FatJetMu"      ) ;
  FatJetPhAccept_  = iConfig.getUntrackedParameter<bool>("skimEvents_FatJetPh"      ) ;

  electronCollectionLabel_     = iConfig.getParameter<edm::InputTag>("electronCollection"      ) ;
  muonCollectionLabel_         = iConfig.getParameter<edm::InputTag>("muonCollection"          ) ;
  tauCollectionLabel_          = iConfig.getParameter<edm::InputTag>("tauCollection");
  photonCollectionLabel_       = iConfig.getParameter<edm::InputTag>("photonCollection"        ) ;
  FatJetLabel_                  = iConfig.getParameter<edm::InputTag>("DeepAK8JetCollection"                 );

  electronCollectionToken_ =  iC.consumes<View<pat::Electron> > (electronCollectionLabel_);
  muonCollectionToken_     =  iC.consumes<View<pat::Muon> > (muonCollectionLabel_);
  tauCollectionToken_      =  iC.consumes<View<pat::Tau>> (tauCollectionLabel_);
  photonCollectionToken_   =  iC.consumes<View<pat::Photon> > (photonCollectionLabel_);
  FatJetToken_                  = iC.consumes<View<pat::Jet> >       ( FatJetLabel_                      );
}
IIHEModuleSkimEvents::~IIHEModuleSkimEvents(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleSkimEvents::beginJob(){
  nAcceptAll_ = 0 ;
}

// ------------ method called to for each event  ------------
void IIHEModuleSkimEvents::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Electron> > electronCollection_;
  iEvent.getByToken( electronCollectionToken_, electronCollection_) ;

  edm::Handle<edm::View<pat::Muon> > muonCollection_;
  iEvent.getByToken( muonCollectionToken_, muonCollection_) ;

  Handle<View<pat::Tau> > tauCollection_ ;
  iEvent.getByToken( tauCollectionToken_, tauCollection_ );

  edm::Handle<edm::View<pat::Photon> > photonCollection_;
  iEvent.getByToken( photonCollectionToken_, photonCollection_) ;

  edm::Handle<edm::View<pat::Jet> > FatJetHandle_;
    iEvent.getByToken(FatJetToken_ ,FatJetHandle_);

  int nEl = 0 ;
  int nMu = 0 ;
  int nTau = 0 ;
  int nPh = 0 ;
  int nEmu = 0;
  int nETau = 0;
  int nMuTau = 0;
  int nFatJet  = 0;
  
  for( unsigned int i = 0 ; i < electronCollection_->size() ; i++ ) {
    Ptr<pat::Electron> gsfiter = electronCollection_->ptrAt( i );
    float pt = gsfiter->pt() ;
    float HEEP_ET  = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    if(pt>ptThresholdEle_ || HEEP_ET>ptThresholdEle_) nEl++ ;
  }

  for( unsigned int i = 0 ; i < muonCollection_->size() ; i++ ) {
    Ptr<pat::Muon> muiter = muonCollection_->ptrAt( i );
    float pt = muiter->pt() ;
    if(pt>ptThresholdmu_ && muiter->isGlobalMuon()) {
      nMu++ ;
    }
  }

  for ( unsigned int i = 0; i <tauCollection_->size(); ++i) {
    Ptr<pat::Tau> tauni = tauCollection_->ptrAt( i );
    float pt = tauni->pt();
    if(pt>ptThresholdTau_) nTau++ ;
  }

  for( unsigned int i = 0 ; i < photonCollection_->size() ; i++ ) {
    Ptr<pat::Photon> phiter = photonCollection_->ptrAt( i );
    if(phiter->pt() > ptThresholdPh_ && phiter->hadTowOverEm()<0.05) {
      nPh++ ;
    }
  }

  for ( unsigned int i = 0; i <FatJetHandle_->size(); ++i) {
    Ptr<pat::Jet> fatjet = FatJetHandle_->ptrAt(i);
    if(fatjet->pt() > ptThresholdFatJet_) nFatJet++;
  }

  nEmu = nEl + nMu;
  nETau = nEl + nTau;
  nMuTau = nMu + nTau;

  bool acceptThisEvent=false;
  if(EleAccept_      && nEl >0      ) acceptThisEvent=true; 
  if(EleMuAccept_    && nEmu >0     ) acceptThisEvent=true;
  if(EleTauAccept_   && nETau >0    ) acceptThisEvent=true; 
  if(MuAccept_       && nMu>0       ) acceptThisEvent=true;
  if(MuTauAccept_    && nMuTau>0    ) acceptThisEvent=true;
  if(TauAccept_      && nTau>0      ) acceptThisEvent=true;
  if(PhAccept_       && nPh>0       ) acceptThisEvent=true;
  if(FatJetAccept_   && nFatJet>0   ) acceptThisEvent=true;
  if(FatJetAccept_   && nFatJet>0   ) acceptThisEvent=true;
  if(FatJetMuAccept_   && nFatJet>1 && nMu>0) acceptThisEvent=true;
  if(FatJetPhAccept_   && nFatJet>0 && nPh>0) acceptThisEvent=true;

//cout<<"nFatJet="<<nFatJet<<",nMu="<<nMu<<",nPh="<<nPh<<endl;
//  acceptThisEvent = (nEl >= nEleAccept_ || nMu >= nMuAccept_ || nTau >= nTauAccept_ || nPh >= nPhAccept_ || (nEmu >= nEleMuAccept_ && nEl >0 &&  nMu>0) || (nETau>= nEleTauAccept_  && nEl >0 && nTau>0) || (nMuTau>= nMuTauAccept_  && nMu>0 && nTau>0) || nFatJet >= nFatJetAccept_);
  
  

  // Save the event if we see something we like
  if(acceptThisEvent){
    acceptEvent() ;
    nAcceptAll_++ ;
  }
  else{
    //rejectEvent() ;
  }
}
void IIHEModuleSkimEvents::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleSkimEvents::beginEvent(){}
void IIHEModuleSkimEvents::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleSkimEvents::endJob(){
  std::cout << std::endl << "IIHEModuleSkimEvents report:" << std::endl ;
  std::cout << "  nAcceptAll  = " << nAcceptAll_   << std::endl ;
}

DEFINE_FWK_MODULE(IIHEModuleSkimEvents);
