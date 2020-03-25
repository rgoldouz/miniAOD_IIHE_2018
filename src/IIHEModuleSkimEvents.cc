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

  nEleThreshold_       = iConfig.getUntrackedParameter<int>("skimEvents_nEle"     ) ;
  nEleMuThreshold_     = iConfig.getUntrackedParameter<int>("skimEvents_nEleMu"   ) ;
  nEleTauThreshold_    = iConfig.getUntrackedParameter<int>("skimEvents_nEleTau"  ) ;
  nMuThreshold_        = iConfig.getUntrackedParameter<int>("skimEvents_nMu"      ) ;
  nMuTauThreshold_     = iConfig.getUntrackedParameter<int>("skimEvents_nMuTau"   ) ;
  nTauThreshold_       = iConfig.getUntrackedParameter<int>("skimEvents_nTau"     ) ;
  nPhThreshold_        = iConfig.getUntrackedParameter<int>("skimEvents_nPh"      ) ;

  electronCollectionLabel_     = iConfig.getParameter<edm::InputTag>("electronCollection"      ) ;
  muonCollectionLabel_         = iConfig.getParameter<edm::InputTag>("muonCollection"          ) ;
  tauCollectionLabel_          = iConfig.getParameter<edm::InputTag>("tauCollection");
  photonCollectionLabel_       = iConfig.getParameter<edm::InputTag>("photonCollection"        ) ;

  electronCollectionToken_ =  iC.consumes<View<pat::Electron> > (electronCollectionLabel_);
  muonCollectionToken_     =  iC.consumes<View<pat::Muon> > (muonCollectionLabel_);
  tauCollectionToken_      =  iC.consumes<View<pat::Tau>> (tauCollectionLabel_);
  photonCollectionToken_   =  iC.consumes<View<pat::Photon> > (photonCollectionLabel_);
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

  int nEl = 0 ;
  int nMu = 0 ;
  int nTau = 0 ;
  int nPh = 0 ;
  int nEmu = 0;
  int nETau = 0;
  int nMuTau = 0;
  
  for( unsigned int i = 0 ; i < electronCollection_->size() ; i++ ) {
    Ptr<pat::Electron> gsfiter = electronCollection_->ptrAt( i );
    float pt = gsfiter->pt() ;
    float HEEP_ET  = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    if(pt>ptThresholdEle_ || HEEP_ET>ptThresholdEle_) nEl++ ;
  }

  for( unsigned int i = 0 ; i < muonCollection_->size() ; i++ ) {
    Ptr<pat::Muon> muiter = muonCollection_->ptrAt( i );
    float pt = muiter->pt() ;
    if(pt>ptThresholdmu_) nMu++ ;
  }

  for ( unsigned int i = 0; i <tauCollection_->size(); ++i) {
    Ptr<pat::Tau> tauni = tauCollection_->ptrAt( i );
    float pt = tauni->pt();
    if(pt>ptThresholdTau_) nTau++ ;
  }

  for( unsigned int i = 0 ; i < photonCollection_->size() ; i++ ) {
    Ptr<pat::Photon> phiter = photonCollection_->ptrAt( i );
    if(phiter->pt() > ptThresholdPh_) nPh++ ;
  }

  nEmu = nEl + nMu;
  nETau = nEl + nTau;
  nMuTau = nMu + nTau;

  bool acceptThisEvent = (nEl >= nEleThreshold_ || nMu >= nMuThreshold_ || nTau >= nTauThreshold_ || nPh >= nPhThreshold_ || (nEmu >= nEleMuThreshold_ && nEl >0 &&  nMu>0) || (nETau>= nEleTauThreshold_  && nEl >0 && nTau>0) || (nMuTau>= nMuTauThreshold_  && nMu>0 && nTau>0));
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
