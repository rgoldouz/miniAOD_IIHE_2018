#include "UserCode/IIHETree/interface/IIHEModuleGedGsfElectron.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include <iostream>
#include <TMath.h>
#include <vector>
#include <Math/VectorUtil.h>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleGedGsfElectron::IIHEModuleGedGsfElectron(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  ebReducedRecHitCollection_ = iC.consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_ = iC.consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_ = iC.consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("esReducedRecHitCollection"));
  beamSpotToken_      = iC.consumes<reco::BeamSpot>(iConfig.getParameter<InputTag>("beamSpot")) ;
  rhoTokenAll_ = iC.consumes<double> (iConfig.getParameter<edm::InputTag>("eventRho"));
  electronCollectionToken_     = iC.consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("electronCollection")) ;
  ETThreshold_ = iConfig.getUntrackedParameter<double>("electronPtThreshold") ;
  primaryVertexLabel_          = iConfig.getParameter<edm::InputTag>("primaryVertex") ;
  vtxToken_ = iC.consumes<View<reco::Vertex>>(primaryVertexLabel_);
}
IIHEModuleGedGsfElectron::~IIHEModuleGedGsfElectron(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleGedGsfElectron::beginJob(){
  addBranch("gsf_n", kUInt) ;
  addBranch("gsf_classification", kVectorInt) ;

  setBranchType(kVectorFloat) ;
  addBranch("gsf_ecalTrkEnergyPostCorr") ;
  addBranch("gsf_ecalEnergyPostCorr") ;
  addBranch("gsf_energy") ;
  addBranch("gsf_p") ;
  addBranch("gsf_pt") ;
  addBranch("gsf_et") ;


  setBranchType(kVectorFloat) ;
  addBranch("gsf_energy") ;
  addBranch("gsf_p") ;
  addBranch("gsf_pt") ;
  addBranch("gsf_et") ;
  addBranch("gsf_scE1x5") ;
  addBranch("gsf_scE5x5") ;
  addBranch("gsf_scE2x5Max") ;
  addBranch("gsf_full5x5_e5x5") ;
  addBranch("gsf_full5x5_e1x5") ;
  addBranch("gsf_full5x5_e2x5Max") ;
  addBranch("gsf_full5x5_sigmaIetaIeta");
  addBranch("gsf_full5x5_hcalOverEcal");

  addBranch("gsf_eta") ;
  addBranch("gsf_phi") ;
  addBranch("gsf_theta") ;
  addBranch("gsf_px") ;
  addBranch("gsf_py") ;
  addBranch("gsf_pz") ;
  addBranch("gsf_caloEnergy") ;
  addBranch("gsf_deltaEtaSuperClusterTrackAtVtx") ;
  addBranch("gsf_deltaPhiSuperClusterTrackAtVtx") ;
  addBranch("gsf_hadronicOverEm") ;
  addBranch("gsf_hcalDepth1OverEcal") ;
  addBranch("gsf_hcalDepth2OverEcal") ;
  addBranch("gsf_dr03TkSumPt") ;
  addBranch("gsf_heepTrkPtIso") ;
  addBranch("gsf_dr03EcalRecHitSumEt") ;
  addBranch("gsf_dr03HcalDepth1TowerSumEt") ;
  addBranch("gsf_dr03HcalDepth2TowerSumEt") ;
  addBranch("gsf_charge", kVectorInt) ;
  addBranch("gsf_sigmaIetaIeta") ;
  setBranchType(kVectorInt) ;
  addBranch("gsf_ecaldrivenSeed"   ) ;
  addBranch("gsf_trackerdrivenSeed") ;
  addBranch("gsf_isEB") ;
  addBranch("gsf_isEE") ;
  addBranch("gsf_passConversionVeto") ;
  addBranch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose");
  addBranch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium"); 
  addBranch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight");
  addBranch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto"); 
  addBranch("gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose");
  addBranch("gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium"); 
  addBranch("gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight");
  addBranch("gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto");
  addBranch("gsf_VID_heepElectronID_HEEPV70");
  addBranch("gsf_VID_mvaEleID_Fall17_iso_V1_wp80");
  addBranch("gsf_VID_mvaEleID_Fall17_iso_V1_wp90"); 
  addBranch("gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose");
  addBranch("gsf_VID_mvaEleID_Fall17_noIso_V1_wp80");
  addBranch("gsf_VID_mvaEleID_Fall17_noIso_V1_wp90");
  addBranch("gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose");
  addBranch("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80"); 
  addBranch("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90");
  addBranch("gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose");

  setBranchType(kVectorFloat) ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaPhiSeedClusterTrackAtCalo") ;
  addBranch("gsf_ecalEnergy") ;
  addBranch("gsf_eSuperClusterOverP") ;
  addBranch("gsf_dxy") ;
  addBranch("gsf_dxy_beamSpot") ;
  addBranch("gsf_dxy_firstPVtx") ;
  addBranch("gsf_dxyError") ;
  addBranch("gsf_dz") ;
  addBranch("gsf_dz_beamSpot") ;
  addBranch("gsf_dz_firstPVtx") ;
  addBranch("gsf_dzError") ;
  addBranch("gsf_vz") ;
  setBranchType(kVectorInt) ;
  addBranch("gsf_numberOfValidHits") ;
  addBranch("gsf_nLostInnerHits"   ) ;
  addBranch("gsf_nLostOuterHits"   ) ;
  addBranch("gsf_convFlags"        ) ;
  setBranchType(kVectorFloat) ;
  addBranch("gsf_convDist") ;
  addBranch("gsf_convDcot") ;
  addBranch("gsf_convRadius") ;
  addBranch("gsf_fBrem") ;
  addBranch("gsf_e1x5") ;
  addBranch("gsf_e2x5Max") ;
  addBranch("gsf_e5x5") ;
  addBranch("gsf_r9") ;
  addBranch("gsf_deltaPhiSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtVtx") ;
  addBranch("gsf_relIso") ;
  addBranch("gsf_effArea") ;
  addBranch("gsf_sumChargedHadronPt") ;
  addBranch("gsf_sumNeutralHadronEt") ;
  addBranch("gsf_sumPhotonEt") ;
  addBranch("gsf_ooEmooP") ;
  addBranch("gsf_eSuperClusterOverP") ;

  addBranch("gsf_hitsinfo", kVectorVectorInt) ;

  setBranchType(kVectorFloat) ;
  addBranch("gsf_pixelMatch_dPhi1") ;
  addBranch("gsf_pixelMatch_dPhi2") ;
  addBranch("gsf_pixelMatch_dRz1" ) ;
  addBranch("gsf_pixelMatch_dRz2" ) ;
  setBranchType(kVectorInt) ;
  addBranch("gsf_pixelMatch_subDetector1") ;
  addBranch("gsf_pixelMatch_subDetector2") ;
  addBranch("gsf_mc_bestDR", kVectorFloat) ;
  addBranch("gsf_mc_index" , kVectorInt  ) ;
  addBranch("gsf_mc_ERatio", kVectorFloat) ;

  setBranchType(kVectorFloat) ;
  addBranch("gsf_sc_energy") ;
  addBranch("gsf_sc_seed_eta") ;
  addBranch("gsf_sc_eta") ;
  addBranch("gsf_sc_etacorr") ;
  addBranch("gsf_sc_theta") ;
  addBranch("gsf_sc_thetacorr") ;
  addBranch("gsf_sc_et") ;
  addBranch("gsf_sc_phi") ;
  addBranch("gsf_sc_px") ;
  addBranch("gsf_sc_py") ;
  addBranch("gsf_sc_pz") ;
  addBranch("gsf_sc_x") ;
  addBranch("gsf_sc_y") ;
  addBranch("gsf_sc_z") ;
  addBranch("gsf_sc_phiWidth") ;
  addBranch("gsf_sc_etaWidth") ;
  addBranch("gsf_sc_seed_rawId", kVectorInt) ;
  addBranch("gsf_sc_seed_ieta", kVectorInt) ;
  addBranch("gsf_sc_seed_iphi", kVectorInt) ;

  setBranchType(kVectorInt) ;
  addBranch("gsf_sc_seed_kHasSwitchToGain6") ;
  addBranch("gsf_sc_seed_kHasSwitchToGain1") ;

  setBranchType(kVectorFloat) ;
  addBranch("gsf_swissCross") ;
  addBranch("gsf_sc_rawEnergy") ;
  addBranch("gsf_sc_preshowerEnergy") ;
  addBranch("gsf_sc_lazyTools_e2x5Right") ;
  addBranch("gsf_sc_lazyTools_e2x5Left") ;
  addBranch("gsf_sc_lazyTools_e2x5Top") ;
  addBranch("gsf_sc_lazyTools_e2x5Bottom") ;
  addBranch("gsf_sc_lazyTools_eMax") ;
  addBranch("gsf_sc_lazyTools_e2nd") ;
  addBranch("gsf_sc_lazyTools_eRight") ;
  addBranch("gsf_sc_lazyTools_eLeft") ;
  addBranch("gsf_sc_lazyTools_eTop") ;
  addBranch("gsf_sc_lazyTools_eBottom") ;
  addBranch("gsf_sc_lazyTools_e2x2") ;
  addBranch("gsf_sc_lazyTools_e3x3") ;
  addBranch("gsf_sc_lazyTools_e4x4") ;
  addBranch("gsf_sc_lazyTools_e5x5") ;
  addBranch("gsf_sc_lazyTools_e1x3") ;
  addBranch("gsf_sc_lazyTools_e3x1") ;
  addBranch("gsf_sc_lazyTools_e1x5") ;
  addBranch("gsf_sc_lazyTools_e5x1") ;
  addBranch("gsf_sc_lazyTools_eshitsixix") ;
  addBranch("gsf_sc_lazyTools_eshitsiyiy") ;
  addBranch("gsf_sc_lazyTools_eseffsixix") ;
  addBranch("gsf_sc_lazyTools_eseffsiyiy") ;
  addBranch("gsf_sc_lazyTools_eseffsirir") ;
  addBranch("gsf_sc_lazyTools_BasicClusterSeedTime") ;

  // Saturation information
  addBranch("EHits_isSaturated", kInt) ;
  setBranchType(kVectorInt) ;
  addBranch("EBHits_rawId"   ) ;
  addBranch("EBHits_iRechit" ) ;
  addBranch("EBHits_ieta"    ) ;
  addBranch("EBHits_iphi"    ) ;
  addBranch("EBHits_RecoFlag") ;
  addBranch("EBHits_kSaturated"           ) ;
  addBranch("EBHits_kLeadingEdgeRecovered") ;
  addBranch("EBHits_kNeighboursRecovered" ) ;
  addBranch("EBHits_kWeird"               ) ;
  addBranch("EBHits_energy", kVectorFloat) ;

  setBranchType(kVectorInt) ;
  addBranch("EEHits_rawId"   ) ;
  addBranch("EEHits_iRechit" ) ;
  addBranch("EEHits_ieta"    ) ;
  addBranch("EEHits_iphi"    ) ;
  addBranch("EEHits_RecoFlag") ;
  addBranch("EEHits_kSaturated"           ) ;
  addBranch("EEHits_kLeadingEdgeRecovered") ;
  addBranch("EEHits_kNeighboursRecovered" ) ;
  addBranch("EEHits_kWeird"               ) ;
  addBranch("EEHits_energy", kVectorFloat) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleGedGsfElectron::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Delegate default electron collection name to IIHEAnalysis class
 
  Handle<EcalRecHitCollection> EBHits;
  Handle<EcalRecHitCollection> EEHits;
  Handle<EcalRecHitCollection> ESHits;
  iEvent.getByToken(ebReducedRecHitCollection_, EBHits) ;
  iEvent.getByToken(eeReducedRecHitCollection_, EEHits) ;
  iEvent.getByToken(esReducedRecHitCollection_, ESHits) ;
   
  const EcalRecHitCollection* theBarrelEcalRecHits = EBHits.product () ;
  const EcalRecHitCollection* theEndcapEcalRecHits = EEHits.product () ;

  edm::Handle<edm::View<pat::Electron>>     electronCollection_ ;
  iEvent.getByToken( electronCollectionToken_ , electronCollection_) ;

  edm::Handle<reco::BeamSpot> beamspotHandle_ ;
  iEvent.getByToken(beamSpotToken_, beamspotHandle_) ;

  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord>().get(pGeometry) ;
  CaloGeometry* geometry = (CaloGeometry*) pGeometry.product() ;
  const CaloSubdetectorGeometry* geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower) ;
  CaloSubdetectorTopology* topology_ES = (geometryES) ? new EcalPreshowerTopology(geometry) : 0 ;

  edm::Handle<View<reco::Vertex> > pvCollection_ ;
  iEvent.getByToken( vtxToken_ , pvCollection_);
  edm::Ptr<reco::Vertex> firstpvertex = pvCollection_->ptrAt( 0 );

  float pv_z = firstpvertex->z() ; 
  EcalClusterLazyTools lazytool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_) ;

  edm::Handle<double> rhoHandle ;
  iEvent.getByToken(rhoTokenAll_, rhoHandle) ;
  double rho = *rhoHandle ;
 
  unsigned int gsf_n = 0 ;
  unsigned int gsfref = -1 ;

  MiniAODHelper electronHelper;
  electronHelper.SetRho(rho);
  electronHelper.SetVertex(pvCollection_->at(0));

  double EffArea = 9999.;

  for( unsigned int i = 0 ; i < electronCollection_->size() ; i++ ) {
    Ptr<pat::Electron> gsfiter = electronCollection_->ptrAt( i );

    gsfref++;
    float ET = gsfiter->caloEnergy()*sin(2.*atan(exp(-1.*gsfiter->eta()))) ;
    if(ET<ETThreshold_ && gsfiter->pt()<ETThreshold_) continue ;
    gsf_n++ ;

    float sc_energy = gsfiter->superCluster()->energy();
    float sc_et     = sc_energy*sin(2.*atan(exp(-1.*gsfiter->superCluster()->eta()))) ;
    float etaCorr = etacorr( gsfiter->superCluster()->eta(), pv_z, gsfiter->superCluster()->position().z()) ;
    double Eta = abs(gsfiter->superCluster()->eta());
    if (Eta >= 0. && Eta < 1.0) EffArea = 0.1703;
    else if (Eta >= 1.0 && Eta < 1.479) EffArea = 0.1715;
    else if (Eta >= 1.479 && Eta < 2.0) EffArea = 0.1213;
    else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.1230;
    else if (Eta >= 2.2 && Eta < 2.3) EffArea = 0.1635;
    else if (Eta >= 2.3 && Eta < 2.4) EffArea = 0.1937;
    else if (Eta >= 2.4 && Eta < 5) EffArea = 0.2393;
    int gsf_nLostInnerHits = gsfiter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) ;
    int gsf_nLostOuterHits = gsfiter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) ;
    store("gsf_energy"                        , gsfiter->energy()                        ) ;
    store("gsf_ecalTrkEnergyPostCorr"         , gsfiter->userFloat("ecalTrkEnergyPostCorr")) ;
    store("gsf_ecalEnergyPostCorr"            , gsfiter->userFloat("ecalEnergyPostCorr")) ;
    store("gsf_p"                             , gsfiter->p()                             ) ;
    store("gsf_pt"                            , gsfiter->pt()                            ) ;
    store("gsf_et"                            , gsfiter->et()                            ) ;
    store("gsf_classification"                , gsfiter->classification()                ) ;
    store("gsf_scE1x5"                        , gsfiter->scE1x5()                        ) ;
    store("gsf_scE5x5"                        , gsfiter->scE5x5()                        ) ;
    store("gsf_scE2x5Max"                     , gsfiter->scE2x5Max()                     ) ;
    store("gsf_eta"                           , gsfiter->eta()                           ) ;
    store("gsf_phi"                           , gsfiter->phi()                           ) ;
    store("gsf_theta"                         , gsfiter->theta()                         ) ;
    store("gsf_px"                            , gsfiter->px()                            ) ;
    store("gsf_py"                            , gsfiter->py()                            ) ;
    store("gsf_pz"                            , gsfiter->pz()                            ) ;
    store("gsf_caloEnergy"                    , gsfiter->caloEnergy()                    ) ;
    store("gsf_deltaEtaSuperClusterTrackAtVtx", gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
    store("gsf_deltaPhiSuperClusterTrackAtVtx", gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
    store("gsf_hadronicOverEm"                , gsfiter->hadronicOverEm()                ) ;
    store("gsf_hcalDepth1OverEcal"            , gsfiter->hcalDepth1OverEcal()            ) ;
    store("gsf_hcalDepth2OverEcal"            , gsfiter->hcalDepth2OverEcal()            ) ;
    store("gsf_dr03TkSumPt"                   , gsfiter->dr03TkSumPt()                   ) ;
    store("gsf_heepTrkPtIso"              , gsfiter->userFloat("heepTrkPtIso"        ) );
    store("gsf_relIso"                        , electronHelper.GetElectronRelIso(electronCollection_->at(i), coneSize::R03, corrType::rhoEA, effAreaType::spring16)) ;
    store("gsf_effArea"                       , EffArea                                  ) ;
    store("gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose"   ,int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-loose")));
    store("gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium"  ,int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-medium")));
    store("gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight"   ,int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-tight")));
    store("gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto"    ,int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-veto")));
    store("gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose" ,int(gsfiter->electronID("cutBasedElectronID-Summer16-80X-V1-loose")));
    store("gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium",int(gsfiter->electronID("cutBasedElectronID-Summer16-80X-V1-medium")));
    store("gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight" ,int(gsfiter->electronID("cutBasedElectronID-Summer16-80X-V1-tight")));
    store("gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto"  ,int(gsfiter->electronID("cutBasedElectronID-Summer16-80X-V1-veto")));
    store("gsf_VID_heepElectronID_HEEPV70"                   ,int(gsfiter->electronID("heepElectronID-HEEPV70")));
    store("gsf_VID_mvaEleID_Fall17_iso_V1_wp80"              ,int(gsfiter->electronID("mvaEleID-Fall17-iso-V1-wp80")));
    store("gsf_VID_mvaEleID_Fall17_iso_V1_wp90"              ,int(gsfiter->electronID("mvaEleID-Fall17-iso-V1-wp90")));
    store("gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose"           ,int(gsfiter->electronID("mvaEleID-Fall17-iso-V1-wpLoose")));
    store("gsf_VID_mvaEleID_Fall17_noIso_V1_wp80"            ,int(gsfiter->electronID("mvaEleID-Fall17-noIso-V1-wp80")));
    store("gsf_VID_mvaEleID_Fall17_noIso_V1_wp90"            ,int(gsfiter->electronID("mvaEleID-Fall17-noIso-V1-wp90")));
    store("gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose"         ,int(gsfiter->electronID("mvaEleID-Fall17-noIso-V1-wpLoose")));
    store("gsf_VID_mvaEleIDSpring16_GeneralPurpose_V1_wp80" ,int(gsfiter->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp80")));
    store("gsf_VID_mvaEleIDSpring16_GeneralPurpose_V1_wp90" ,int(gsfiter->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp90")));
    store("gsf_VID_mvaEleIDSpring16_HZZ_V1_wpLoose"         ,int(gsfiter->electronID("mvaEleID-Spring16-HZZ-V1-wpLoose")));
    store("gsf_dr03EcalRecHitSumEt"           , gsfiter->dr03EcalRecHitSumEt()           ) ;
    store("gsf_dr03HcalDepth1TowerSumEt"      , gsfiter->dr03HcalDepth1TowerSumEt()      ) ;
    store("gsf_dr03HcalDepth2TowerSumEt"      , gsfiter->dr03HcalDepth2TowerSumEt()      ) ;
    store("gsf_charge"                        , gsfiter->charge()                        ) ;
    store("gsf_sigmaIetaIeta"                 , gsfiter->sigmaIetaIeta()                 ) ;
    store("gsf_ecaldrivenSeed"                , int(gsfiter->ecalDrivenSeed())                ) ;
    store("gsf_trackerdrivenSeed"             , int(gsfiter->trackerDrivenSeed())             ) ;
    store("gsf_isEB"                          , int(gsfiter->isEB())                          ) ;
    store("gsf_isEE"                          , int(gsfiter->isEE())                          ) ;
    store("gsf_passConversionVeto"            , gsfiter->passConversionVeto()            ) ;
    store("gsf_deltaEtaSeedClusterTrackAtCalo", gsfiter->deltaEtaSeedClusterTrackAtCalo()) ;
    store("gsf_deltaPhiSeedClusterTrackAtCalo", gsfiter->deltaPhiSeedClusterTrackAtCalo()) ;
    store("gsf_ecalEnergy"                    , gsfiter->ecalEnergy()                    ) ;
    store("gsf_eSuperClusterOverP"            , gsfiter->eSuperClusterOverP()            ) ;
    store("gsf_dxy"                           , gsfiter->gsfTrack()->dxy()               ) ;
    store("gsf_dxy_beamSpot"                  , gsfiter->gsfTrack()->dxy(beamspotHandle_->position())      ) ;
    store("gsf_dxy_firstPVtx"                 , gsfiter->gsfTrack()->dxy(firstpvertex->position())  ) ;
    store("gsf_dxyError"                      , gsfiter->gsfTrack()->dxyError()          ) ;
    store("gsf_dz"                            , gsfiter->gsfTrack()->dz()                ) ;
    store("gsf_dz_beamSpot"                   , gsfiter->gsfTrack()->dz(beamspotHandle_->position())       ) ;
    store("gsf_dz_firstPVtx"                  , gsfiter->gsfTrack()->dz(firstpvertex->position())   ) ;
    store("gsf_dzError"                       , gsfiter->gsfTrack()->dzError()           ) ; 
    store("gsf_vz"                            , gsfiter->gsfTrack()->vz()                ) ;
    store("gsf_numberOfValidHits"             , gsfiter->gsfTrack()->numberOfValidHits() ) ;
    store("gsf_nLostInnerHits"                , gsf_nLostInnerHits                       ) ;
    store("gsf_nLostOuterHits"                , gsf_nLostOuterHits                       ) ;
    store("gsf_convFlags"                     , gsfiter->convFlags()                     ) ;
    store("gsf_convDist"                      , gsfiter->convDist()                      ) ;
    store("gsf_convDcot"                      , gsfiter->convDcot()                      ) ;
    store("gsf_convRadius"                    , gsfiter->convRadius()                    ) ;
    store("gsf_fBrem"                         , gsfiter->fbrem()                         ) ;
    store("gsf_e1x5"                          , gsfiter->e1x5()                          ) ;
    store("gsf_r9"                            , gsfiter->r9()                            ) ;
    store("gsf_e2x5Max"                       , gsfiter->e2x5Max()                       ) ;
    store("gsf_e5x5"                          , gsfiter->e5x5()                          ) ;
    store("gsf_deltaPhiSeedClusterTrackAtCalo", gsfiter->deltaPhiSeedClusterTrackAtCalo()) ;
    store("gsf_deltaEtaSeedClusterTrackAtCalo", gsfiter->deltaEtaSeedClusterTrackAtCalo()) ;
    store("gsf_deltaEtaSeedClusterTrackAtVtx" , gsfiter->deltaEtaSeedClusterTrackAtVtx() ) ;
    store("gsf_full5x5_e5x5"                  ,gsfiter->full5x5_e5x5()) ;
    store("gsf_full5x5_e1x5"                  ,gsfiter->full5x5_e1x5()) ;
    store("gsf_full5x5_e2x5Max"               ,gsfiter->full5x5_e2x5Max()) ;
    store("gsf_full5x5_sigmaIetaIeta"         ,gsfiter->full5x5_sigmaIetaIeta()) ;
    store("gsf_full5x5_hcalOverEcal"          ,gsfiter->full5x5_hcalOverEcal());
    store("gsf_sumChargedHadronPt"            ,gsfiter->pfIsolationVariables().sumChargedHadronPt) ;
    store("gsf_sumNeutralHadronEt"            ,gsfiter->pfIsolationVariables().sumNeutralHadronEt) ;
    store("gsf_sumPhotonEt"                   ,gsfiter->pfIsolationVariables().sumPhotonEt) ;
    store("gsf_ooEmooP"                       ,fabs(1.0/gsfiter->ecalEnergy() - gsfiter->eSuperClusterOverP()/gsfiter->ecalEnergy() )) ;
    store("gsf_eSuperClusterOverP"            ,gsfiter->eSuperClusterOverP()) ;

    store("gsf_sc_eta"        , gsfiter->superCluster()->eta()                    ) ;
    store("gsf_sc_etacorr"    , etaCorr                                           ) ;
    store("gsf_sc_theta"      , 2.*atan(exp(-1.*gsfiter->superCluster()->eta()))  ) ;
    store("gsf_sc_thetacorr"  , 2.*atan(exp(-1.*etaCorr))                         ) ;
    store("gsf_sc_phi"        , gsfiter->superCluster()->phi()                    ) ;
    store("gsf_sc_energy"     , sc_energy                                         ) ;
    store("gsf_sc_et"         , sc_et                                             ) ;
    store("gsf_sc_px"         , sc_et*cos(gsfiter->superCluster()->phi())         ) ;
    store("gsf_sc_py"         , sc_et*sin(gsfiter->superCluster()->phi())         ) ;
    store("gsf_sc_pz"         , sc_energy*tanh(gsfiter->superCluster()->eta())    ) ;
    store("gsf_sc_x"          , gsfiter->superCluster()->position().x()           ) ;
    store("gsf_sc_y"          , gsfiter->superCluster()->position().y()           ) ;
    store("gsf_sc_z"          , gsfiter->superCluster()->position().z()           ) ;
    store("gsf_sc_phiWidth"   , gsfiter->superCluster()->phiWidth()               ) ;
    store("gsf_sc_etaWidth"   , gsfiter->superCluster()->etaWidth()               ) ;
    store("gsf_sc_seed_eta"   ,gsfiter->superCluster()->seed()->eta()             ) ;
    store("gsf_sc_seed_rawId" , gsfiter->superCluster()->seed()->seed().rawId()   );
    const EcalRecHitCollection *recHits = (gsfiter->isEB()) ? lazytool.getEcalEBRecHitCollection() : lazytool.getEcalEERecHitCollection();
    EcalRecHitCollection::const_iterator seedRecHit = recHits->find(gsfiter->superCluster()->seed()->seed()) ;
    store("gsf_sc_seed_kHasSwitchToGain6" , seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain6) ) ;
    store("gsf_sc_seed_kHasSwitchToGain1" , seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain1)  ) ;


    const std::vector<std::pair<DetId,float> > & hits= gsfiter->superCluster()->hitsAndFractions();
    if (gsfiter->isEB()){
      EBDetId EBscID = EBDetId(gsfiter->superCluster()->seed()->seed().rawId());
      store("gsf_sc_seed_ieta" , EBscID.ieta()   );
      store("gsf_sc_seed_iphi" , EBscID.iphi()   );
      std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits, theBarrelEcalRecHits); 
      store("gsf_swissCross"  , EcalTools::swissCross(id.first,*theBarrelEcalRecHits,0.));
    }
    else if (gsfiter->isEE()){
      EEDetId EEscID = EEDetId(gsfiter->superCluster()->seed()->seed().rawId());
      store("gsf_sc_seed_ieta" , EEscID.iy()   );
      store("gsf_sc_seed_iphi" , EEscID.ix()   );
      std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits, theEndcapEcalRecHits); 
      store("gsf_swissCross"  , EcalTools::swissCross(id.first,*theEndcapEcalRecHits,0.));
    }

    store("gsf_sc_rawEnergy"          , gsfiter->superCluster()->rawEnergy()      ) ;
    store("gsf_sc_preshowerEnergy"          , gsfiter->superCluster()->preshowerEnergy()) ;

    reco::SuperClusterRef    cl_ref = gsfiter->superCluster() ;
    const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed() ;

    store("gsf_sc_lazyTools_e2x5Right"   , lazytool.e2x5Right (*seed)            );
    store("gsf_sc_lazyTools_e2x5Left"    , lazytool.e2x5Left (*seed)             );
    store("gsf_sc_lazyTools_e2x5Top"     , lazytool.e2x5Top (*seed)              );
    store("gsf_sc_lazyTools_e2x5Bottom"  , lazytool.e2x5Bottom (*seed)           );
    store("gsf_sc_lazyTools_eMax"        , lazytool.eMax (*seed)                 );
    store("gsf_sc_lazyTools_e2nd"        , lazytool.e2nd (*seed)                 );
    store("gsf_sc_lazyTools_eRight"      , lazytool.eRight (*seed)               );
    store("gsf_sc_lazyTools_eLeft"       , lazytool.eLeft (*seed)                );
    store("gsf_sc_lazyTools_eTop"        , lazytool.eTop (*seed)                 );
    store("gsf_sc_lazyTools_eBottom"     , lazytool.eBottom (*seed)              );
    store("gsf_sc_lazyTools_e2x2"        , lazytool.e2x2 (*seed)                 );
    store("gsf_sc_lazyTools_e3x3"        , lazytool.e3x3 (*seed)                 );
    store("gsf_sc_lazyTools_e4x4"        , lazytool.e4x4 (*seed)                 );
    store("gsf_sc_lazyTools_e5x5"        , lazytool.e5x5 (*seed)                 );
    store("gsf_sc_lazyTools_e1x5"        , lazytool.e1x5 (*seed)                 );
    store("gsf_sc_lazyTools_e5x1"        , lazytool.e5x1 (*seed)                 );
    store("gsf_sc_lazyTools_e1x3"        , lazytool.e1x3 (*seed)                 );
    store("gsf_sc_lazyTools_e3x1"        , lazytool.e3x1 (*seed)                 );
    store("gsf_sc_lazyTools_BasicClusterSeedTime"        , lazytool.BasicClusterSeedTime (*seed)  );
    double x = gsfiter->superCluster()->x() ;
    double y = gsfiter->superCluster()->y() ;
    double z = gsfiter->superCluster()->z() ;
    store("gsf_sc_lazyTools_eshitsixix", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_ES, 0, 1)) ;
    store("gsf_sc_lazyTools_eshitsiyiy", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_ES, 0, 2)) ;
    store("gsf_sc_lazyTools_eseffsixix", lazytool.eseffsixix(*cl_ref)) ;
    store("gsf_sc_lazyTools_eseffsiyiy", lazytool.eseffsiyiy(*cl_ref)) ;
    store("gsf_sc_lazyTools_eseffsirir", lazytool.eseffsirir(*cl_ref)) ;

    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();
    int nbtrackhits = kfHitPattern.numberOfAllHits(reco::HitPattern::TRACK_HITS) ;
    std::vector<int> gsf_hitsinfo ;
    for(int hititer=0 ; hititer<25 ; hititer++){
      
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(reco::HitPattern::TRACK_HITS, hititer) : 0 ;
      
      gsf_hitsinfo.push_back(myhitbin) ;
    }
    store("gsf_hitsinfo", gsf_hitsinfo) ;
    
    store("gsf_pixelMatch_dPhi1"       , gsfiter->pixelMatchDPhi1()       ) ;
    store("gsf_pixelMatch_dPhi2"       , gsfiter->pixelMatchDPhi2()       ) ;
    store("gsf_pixelMatch_dRz1"        , gsfiter->pixelMatchDRz1()        ) ;
    store("gsf_pixelMatch_dRz2"        , gsfiter->pixelMatchDRz2()        ) ;
    store("gsf_pixelMatch_subDetector1", gsfiter->pixelMatchSubdetector1()) ;
    store("gsf_pixelMatch_subDetector2", gsfiter->pixelMatchSubdetector2()) ;
    
    // Now apply truth matching.
    int index = MCTruth_matchEtaPhi_getIndex(gsfiter->eta(), gsfiter->phi()) ;
    if(index>=0){
      const MCTruthObject* MCTruth = MCTruth_getRecordByIndex(index) ;
      store("gsf_mc_bestDR", deltaR(gsfiter->eta(), gsfiter->phi(), MCTruth->eta(), MCTruth->phi())) ;
      store("gsf_mc_index" , index) ;
      store("gsf_mc_ERatio", gsfiter->energy()/MCTruth->energy()) ;
    }
    else{
      store("gsf_mc_bestDR", 999.0) ;
      store("gsf_mc_index" ,    -1) ;
      store("gsf_mc_ERatio", 999.0) ;
    }
 }
  store("gsf_n", gsf_n) ;


  int nEBRecHits = 0 ;
  bool isSaturated = false;
  for(EcalRecHitCollection::const_iterator EBIt = theBarrelEcalRecHits->begin() ; EBIt!=theBarrelEcalRecHits->end() ; ++EBIt){
    if((*EBIt).checkFlag(EcalRecHit::kSaturated)) isSaturated = true;
   if( (*EBIt).energy() < 200.0 ) continue ;
    nEBRecHits++ ;
    EBDetId elementId = EBIt->id() ;
    store("EBHits_rawId"   , elementId.rawId()) ;
    store("EBHits_iRechit" , nEBRecHits) ;
    store("EBHits_energy"  , (*EBIt).energy() ) ;
    store("EBHits_ieta"    , elementId.ieta() ) ;
    store("EBHits_iphi"    , elementId.iphi() ) ;
    store("EBHits_RecoFlag", (*EBIt).recoFlag() ) ;

    store("EBHits_kSaturated"           , int((*EBIt).checkFlag(EcalRecHit::kSaturated           ))) ;
    store("EBHits_kLeadingEdgeRecovered", int((*EBIt).checkFlag(EcalRecHit::kLeadingEdgeRecovered))) ;
    store("EBHits_kNeighboursRecovered" , int((*EBIt).checkFlag(EcalRecHit::kNeighboursRecovered ))) ;
    store("EBHits_kWeird"               , int((*EBIt).checkFlag(EcalRecHit::kWeird               ))) ;
  }


  int nEERecHits = 0 ;
  for(EcalRecHitCollection::const_iterator EEIt = theEndcapEcalRecHits->begin() ; EEIt!=theEndcapEcalRecHits->end() ; ++EEIt){
    if((*EEIt).checkFlag(EcalRecHit::kSaturated)) isSaturated = true;
    if( (*EEIt).energy() < 200.0 ) continue ;
    nEERecHits++ ;
    EBDetId elementId = EEIt->id() ;
    store("EEHits_rawId"   , elementId.rawId()) ;
    store("EEHits_iRechit" , nEBRecHits) ;
    store("EEHits_energy"  , (*EEIt).energy() ) ;
    store("EEHits_ieta"    , elementId.ieta() ) ;
    store("EEHits_iphi"    , elementId.iphi() ) ;
    store("EEHits_RecoFlag", (*EEIt).recoFlag() ) ;

    store("EEHits_kSaturated"           , int((*EEIt).checkFlag(EcalRecHit::kSaturated           ))) ;
    store("EEHits_kLeadingEdgeRecovered", int((*EEIt).checkFlag(EcalRecHit::kLeadingEdgeRecovered))) ;
    store("EEHits_kNeighboursRecovered" , int((*EEIt).checkFlag(EcalRecHit::kNeighboursRecovered ))) ;
    store("EEHits_kWeird"               , int((*EEIt).checkFlag(EcalRecHit::kWeird               ))) ;
  }
  store("EHits_isSaturated", int(isSaturated));

}

void IIHEModuleGedGsfElectron::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleGedGsfElectron::beginEvent(){}
void IIHEModuleGedGsfElectron::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleGedGsfElectron::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleGedGsfElectron);
