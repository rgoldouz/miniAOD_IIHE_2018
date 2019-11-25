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
  addBranch("ph_theta") ;
  addBranch("ph_phi") ;
  addBranch("ph_energy") ;
  addBranch("ph_mass") ;
  
  setBranchType(kVectorFloat) ;
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
  
  setBranchType(kVectorInt) ;
  addBranch("ph_hasPixelSeed") ;
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
    store("ph_hasPixelSeed"                  , phiter->hasPixelSeed()                    ) ;

    store("ph_mipChi2"                       , phiter->mipChi2()                         ) ;
    store("ph_mipTotEnergy"                  , phiter->mipTotEnergy()                    ) ;
    store("ph_mipSlope"                      , phiter->mipSlope()                        ) ;
    store("ph_mipIntercept"                  , phiter->mipIntercept()                    ) ;
    store("ph_mipNhitCone"                   , phiter->mipNhitCone()                     ) ;
    store("ph_mipIsHalo"                     , phiter->mipIsHalo()                       ) ;


    store("ph_ecalRecHitSumEtConeDR04"       , phiter->ecalRecHitSumEtConeDR04()         ) ;
    store("ph_hcalTowerSumEtConeDR04"        , phiter->hcalTowerSumEtConeDR04()          ) ;
    store("ph_hcalDepth1TowerSumEtConeDR04"  , phiter->hcalDepth1TowerSumEtConeDR04()    ) ;
    store("ph_hcalDepth2TowerSumEtConeDR04"  , phiter->hcalDepth2TowerSumEtConeDR04()    ) ;
    store("ph_hcalTowerSumEtBcConeDR04"      , phiter->hcalTowerSumEtBcConeDR04()        ) ;
    store("ph_hcalDepth1TowerSumEtBcConeDR04", phiter->hcalDepth1TowerSumEtBcConeDR04()  ) ;
    store("ph_hcalDepth2TowerSumEtBcConeDR04", phiter->hcalDepth2TowerSumEtBcConeDR04()  ) ;
    store("ph_trkSumPtSolidConeDR04"         , phiter->trkSumPtSolidConeDR04()           ) ;
    store("ph_trkSumPtHollowConeDR04"        , phiter->trkSumPtHollowConeDR04()          ) ;
    store("ph_nTrkSolidConeDR04"             , phiter->nTrkSolidConeDR04()               ) ;
    store("ph_nTrkHollowConeDR04"            , phiter->nTrkHollowConeDR04()              ) ;

    store("ph_ecalRecHitSumEtConeDR03"       , phiter->ecalRecHitSumEtConeDR03()         ) ;
    store("ph_hcalTowerSumEtConeDR03"        , phiter->hcalTowerSumEtConeDR03()          ) ;
    store("ph_hcalDepth1TowerSumEtConeDR03"  , phiter->hcalDepth1TowerSumEtConeDR03()    ) ;
    store("ph_hcalDepth2TowerSumEtConeDR03"  , phiter->hcalDepth2TowerSumEtConeDR03()    ) ;
    store("ph_hcalTowerSumEtBcConeDR03"      , phiter->hcalTowerSumEtBcConeDR03()        ) ;
    store("ph_hcalDepth1TowerSumEtBcConeDR03", phiter->hcalDepth1TowerSumEtBcConeDR03()  ) ;
    store("ph_hcalDepth2TowerSumEtBcConeDR03", phiter->hcalDepth2TowerSumEtBcConeDR03()  ) ;
    store("ph_trkSumPtSolidConeDR03"         , phiter->trkSumPtSolidConeDR03()           ) ;
    store("ph_trkSumPtHollowConeDR03"        , phiter->trkSumPtHollowConeDR03()          ) ;
    store("ph_nTrkSolidConeDR03"             , phiter->nTrkSolidConeDR03()               ) ;
    store("ph_nTrkHollowConeDR03"            , phiter->nTrkHollowConeDR03()              ) ;

    store("ph_chargedHadronIso"               , phiter->chargedHadronIso()               ) ;
    store("ph_neutralHadronIso"               , phiter->neutralHadronIso()               ) ;
    store("ph_photonIso"                      , phiter->photonIso()                      ) ;
    
    store("ph_chargedHadronIsoWrongVtx"       , phiter->chargedHadronIsoWrongVtx()       ) ;
    store("ph_sumChargedParticlePt"           , phiter->sumChargedParticlePt()           ) ;
    store("ph_sumNeutralHadronEtHighThreshold", phiter->sumNeutralHadronEtHighThreshold()) ;
    store("ph_sumPhotonEtHighThreshold"       , phiter->sumPhotonEtHighThreshold()       ) ;
    store("ph_sumPUPt"                        , phiter->sumPUPt()                        ) ;

    store("ph_nClusterOutsideMustache"        , phiter->nClusterOutsideMustache()        ) ;
    store("ph_etOutsideMustache"              , phiter->etOutsideMustache()              ) ;
    store("ph_pfMVA"                          , phiter->pfMVA()                          ) ;
    
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
