#include "UserCode/IIHETree/interface/IIHEModulePhoton.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"


#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModulePhoton::IIHEModulePhoton(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  ETThreshold_ = iConfig.getUntrackedParameter<double>("photonPtThreshold") ;
  photonCollectionLabel_       = iConfig.getParameter<edm::InputTag>("photonCollection"        ) ;
  photonCollectionToken_ =  iC.consumes<View<pat::Photon> > (photonCollectionLabel_);
}
IIHEModulePhoton::~IIHEModulePhoton(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModulePhoton::beginJob(){
  addBranch("ph_n", kUInt) ;
  
  setBranchType(kVectorFloat) ;
  addBranch("ph_pt") ;
  addBranch("ph_eta") ;
  addBranch("ph_superCluster_eta") ;
  addBranch("ph_theta") ;
  addBranch("ph_phi") ;
  addBranch("ph_energy") ;
  addBranch("ph_mass") ;
  addBranch("ph_hadronicOverEm") ;
  addBranch("ph_hadronicDepth1OverEm") ;
  addBranch("ph_hadronicDepth2OverEm") ;
  addBranch("ph_hadTowOverEm") ;
  addBranch("ph_hadTowDepth1OverEm") ;
  addBranch("ph_hadTowDepth2OverEm") ; 
  addBranch("ph_e1x5") ;
  addBranch("ph_e2x5") ;
  addBranch("ph_e3x3") ;
  addBranch("ph_e5x5") ;
  addBranch("ph_maxEnergyXtal") ;
  addBranch("ph_sigmaEtaEta") ;
  addBranch("ph_sigmaIetaIeta") ;
  addBranch("ph_r1x5") ;
  addBranch("ph_r2x5") ;
  addBranch("ph_r9") ;
  addBranch("ph_hadronicOverEm") ;
  addBranch("ph_hadronicDepth1OverEm") ;
  addBranch("ph_hadronicDepth2OverEm") ;
  addBranch("ph_hadTowOverEm") ;
  addBranch("ph_hadTowDepth1OverEm") ;
  addBranch("ph_hadTowDepth2OverEm") ;
  addBranch("ph_full5x5_sigmaIetaIeta") ;

  addBranch("ph_phoChargedIsolation") ;
  addBranch("ph_phoNeutralHadronIsolation") ;
  addBranch("ph_phoPhotonIsolation") ;
  addBranch("ph_ecalEnergyErrPreCorr") ;
  addBranch("ph_ecalEnergyPostCorr") ;
  addBranch("ph_energyScaleEtDown") ;
  addBranch("ph_energyScaleEtUp") ;


  setBranchType(kVectorInt) ; 
  addBranch("ph_mvaPhoID_RunIIFall17_v2_wp80") ;
  addBranch("ph_mvaPhoID_RunIIFall17_v2_wp90") ;
  addBranch("ph_cutBasedPhotonID_Fall17_94X_V2_loose") ;
  addBranch("ph_cutBasedPhotonID_Fall17_94X_V2_medium") ;
  addBranch("ph_cutBasedPhotonID_Fall17_94X_V2_tight") ;
  addBranch("ph_hasPixelSeed") ;
  addBranch("ph_isPFlowPhoton") ;
  addBranch("ph_passElectronVeto") ;
  addBranch("ph_isEB") ;
  addBranch("ph_isEE") ;

}

// ------------ method called to for each event  ------------
void IIHEModulePhoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Photon> > photonCollection_;
  iEvent.getByToken( photonCollectionToken_, photonCollection_) ;

  unsigned int ph_n = 0 ;
  for( unsigned int i = 0 ; i < photonCollection_->size() ; i++ ) {
    Ptr<pat::Photon> phiter = photonCollection_->ptrAt( i );
    if(phiter->pt() < ETThreshold_) continue ;
    store("ph_pt"    , phiter->pt()) ;
    store("ph_eta"   , phiter->eta()) ;
    store("ph_superCluster_eta"   , phiter->superCluster()->eta()) ;
    store("ph_theta" , phiter->theta()) ;
    store("ph_phi"   , phiter->phi()) ;
    store("ph_energy", phiter->energy()) ;
    store("ph_mass"  , phiter->mass()) ;

    store("ph_e1x5"                          , phiter->e1x5()                            ) ;
    store("ph_e2x5"                          , phiter->e2x5()                            ) ;
    store("ph_e3x3"                          , phiter->e3x3()                            ) ;
    store("ph_e5x5"                          , phiter->e5x5()                            ) ;
    store("ph_maxEnergyXtal"                 , phiter->maxEnergyXtal()                   ) ;
    store("ph_sigmaEtaEta"                   , phiter->sigmaEtaEta()                     ) ;
    store("ph_sigmaIetaIeta"                 , phiter->sigmaIetaIeta()                   ) ;
    store("ph_r1x5"                          , phiter->r1x5()                            ) ;
    store("ph_r2x5"                          , phiter->r2x5()                            ) ;
    store("ph_r9"                            , phiter->r9()                              ) ;
    store("ph_full5x5_sigmaIetaIeta"         , phiter->full5x5_sigmaIetaIeta()           ) ;

    store("ph_hasPixelSeed"                  , int(phiter->hasPixelSeed())               ) ;
    store("ph_passElectronVeto"              , int(phiter->passElectronVeto())            ) ;
    store("ph_isPFlowPhoton"                 , int(phiter->isPFlowPhoton())              ) ;
    store("ph_hasConversionTracks"           , int(phiter->hasConversionTracks())        ) ;
    store("ph_isEB"                          , int(phiter->isEB())                       ) ;
    store("ph_isEE"                          , int(phiter->isEE())                       ) ;
    store("ph_mvaPhoID_RunIIFall17_v2_wp80"            , int(phiter->photonID("mvaPhoID-RunIIFall17-v2-wp80"))) ;
    store("ph_mvaPhoID_RunIIFall17_v2_wp90"            , int(phiter->photonID("mvaPhoID-RunIIFall17-v2-wp90"))) ;
    store("ph_cutBasedPhotonID_Fall17_94X_V2_loose"    , int(phiter->photonID("cutBasedPhotonID-Fall17-94X-V2-loose"))) ;
    store("ph_cutBasedPhotonID_Fall17_94X_V2_medium"   , int(phiter->photonID("cutBasedPhotonID-Fall17-94X-V2-medium"))) ;
    store("ph_cutBasedPhotonID_Fall17_94X_V2_tight"    , int(phiter->photonID("cutBasedPhotonID-Fall17-94X-V2-tight"))) ;

    store("ph_hadronicOverEm"                , phiter->hadronicOverEm()                  ) ;
    store("ph_hadronicDepth1OverEm"          , phiter->hadronicDepth1OverEm()            ) ;
    store("ph_hadronicDepth2OverEm"          , phiter->hadronicDepth2OverEm()            ) ;
    store("ph_hadTowOverEm"                  , phiter->hadTowOverEm()                    ) ;
    store("ph_hadTowDepth1OverEm"            , phiter->hadTowDepth1OverEm()              ) ;
    store("ph_hadTowDepth2OverEm"            , phiter->hadTowDepth2OverEm()              ) ;


    store("ph_ecalEnergyErrPreCorr"         , phiter->userFloat("ecalEnergyErrPreCorr")) ;
    store("ph_ecalEnergyPostCorr"         , phiter->userFloat("ecalEnergyPostCorr")) ;
    store("ph_energyScaleEtDown"         , phiter->userFloat("energyScaleEtDown")) ;
    store("ph_energyScaleEtUp"         , phiter->userFloat("energyScaleEtUp")) ;

    store("ph_phoChargedIsolation"         , phiter->userFloat("phoChargedIsolation"       )) ;
    store("ph_phoNeutralHadronIsolation"   , phiter->userFloat("phoNeutralHadronIsolation" )) ;
    store("ph_phoPhotonIsolation"          , phiter->userFloat("phoPhotonIsolation"        )) ;


    
    // Now apply truth matching.
    int index = MCTruth_matchEtaPhi_getIndex(phiter->eta(), phiter->phi()) ;
    if(index>=0){
      const MCTruthObject* MCTruth = MCTruth_getRecordByIndex(index) ;
      store("ph_mc_bestDR", deltaR(phiter->eta(), phiter->phi(), MCTruth->eta(), MCTruth->phi())) ;
      store("ph_mc_index" , index) ;
      store("ph_mc_ERatio", phiter->energy()/MCTruth->energy()) ;
    }
    else{
      store("ph_mc_bestDR", 999.0) ;
      store("ph_mc_index" ,    -1) ;
      store("ph_mc_ERatio", 999.0) ;
    }
    ph_n++ ;
  }
  store("ph_n", ph_n ) ;
}

void IIHEModulePhoton::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModulePhoton::beginEvent(){}
void IIHEModulePhoton::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModulePhoton::endJob(){}

DEFINE_FWK_MODULE(IIHEModulePhoton);
