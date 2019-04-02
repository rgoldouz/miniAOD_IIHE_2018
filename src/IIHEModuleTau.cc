#include "UserCode/IIHETree/interface/IIHEModuleTau.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;
using namespace tau;
IIHEModuleTau::IIHEModuleTau(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
	ETThreshold_ = iConfig.getUntrackedParameter<double>("tauPtTThreshold" ) ;
	tauCollectionLabel_     = iConfig.getParameter<edm::InputTag>("tauCollection");
	tauCollectionToken_     = iC.consumes<View<pat::Tau>> (tauCollectionLabel_);
	primaryVertexLabel_          = iConfig.getParameter<edm::InputTag>("primaryVertex") ;
	vtxToken_ = iC.consumes<View<reco::Vertex>>(primaryVertexLabel_);
}
IIHEModuleTau::~IIHEModuleTau(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTau::beginJob(){

	addBranch("tau_n", kUInt);

	setBranchType(kVectorFloat);
	addBranch("tau_px");
	addBranch("tau_py");
	addBranch("tau_pz");
	addBranch("tau_pt");
	addBranch("tau_eta");
	addBranch("tau_theta");
	addBranch("tau_phi");
	addBranch("tau_energy");
	addBranch("tau_mass");
	addBranch("tau_dxy");
	addBranch("tau_dxy_error");
	addBranch("tau_ptLeadChargedCand");
	addBranch("tau_decayModeFinding");
	addBranch("tau_decayModeFindingNewDMs");
	addBranch("tau_againstMuonLoose3");
	addBranch("tau_againstMuonTight3");
	addBranch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
	addBranch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits");
	addBranch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits");
	addBranch("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits");
	addBranch("tau_byIsolationMVArun2v1DBoldDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1DBoldDMwLT");
	addBranch("tau_byIsolationMVArun2v1DBnewDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1DBnewDMwLT");
	addBranch("tau_byIsolationMVArun2v1PWoldDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1PWoldDMwLT");
	addBranch("tau_byIsolationMVArun2v1PWnewDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1PWnewDMwLT");
	addBranch("tau_byIsolationMVArun2v1DBdR03oldDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT");
	addBranch("tau_byIsolationMVArun2v1PWdR03oldDMwLTraw");
	addBranch("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT");
	addBranch("tau_againstElectronMVA6Raw");
	addBranch("tau_againstElectronMVA6category");
	addBranch("tau_againstElectronVLooseMVA6");
	addBranch("tau_againstElectronLooseMVA6");
	addBranch("tau_againstElectronMediumMVA6");
	addBranch("tau_againstElectronTightMVA6");
	addBranch("tau_againstElectronVTightMVA6");
	addBranch("tau_mc_bestDR");
	addBranch("tau_mc_ERatio");

	/*
	   addBranch("tau_byVVLooseIsolationMVArun2v1DBoldDMwLTNew");
	   addBranch("tau_byVLooseIsolationMVArun2v1DBoldDMwLTNew");
	   addBranch("tau_byLooseIsolationMVArun2v1DBoldDMwLTNew");;
	   addBranch("tau_byMediumIsolationMVArun2v1DBoldDMwLTNew") ;
	   addBranch("tau_byTightIsolationMVArun2v1DBoldDMwLTNew");
	   addBranch("tau_byVTightIsolationMVArun2v1DBoldDMwLTNew");
	   addBranch("tau_byVVTightIsolationMVArun2v1DBoldDMwLTNew");
	   addBranch("tau_byIsolationMVArun2v1DBoldDMwLTrawNew");
	   addBranch("tau_againstElectronMVA6RawNew");
	   */

	addBranch("tau_chargedIsoPtSum");
	addBranch("tau_neutralIsoPtSum");
	addBranch("tau_puCorrPtSum");
	addBranch("tau_footprintCorrection");
	addBranch("tau_neutralIsoPtSumWeight");
	addBranch("tau_photonPtSumOutsideSignalCone");
	addBranch("tau_byPhotonPtSumOutsideSignalCone");
	addBranch("tau_footprintCorrectiondR03");
	addBranch("tau_chargedIsoPtSumdR03");
	addBranch("tau_neutralIsoPtSumWeightdR03");
	addBranch("tau_neutralIsoPtSumdR03");
	addBranch("tau_photonPtSumOutsideSignalConedR03");

	// by aman
	addBranch("tau_PFChargedHadIso");
	addBranch("tau_PFNeutralHadIso");
	addBranch("tau_PFPhotonIso");
	addBranch("tau_leadChargedParticlePt");
	addBranch("tau_trackRefPt");            
	addBranch("tau_lead_dxy");
	addBranch("tau_lead_dz");
	addBranch("tau_dxy_Sig");
	addBranch("tau_flightLengthSig");
	addBranch("tau_ip3d");
	addBranch("tau_ip3d_Sig");
	addBranch("tau_decayDistX");
	addBranch("tau_decayDistY");
	addBranch("tau_decayDistZ");
	addBranch("tau_decayDistMag");
	addBranch("tau_nPhoton");
	addBranch("tau_ptWeightedDetaStrip");
	addBranch("tau_ptWeightedDphiStrip");
	addBranch("tau_ptWeightedDrSignal");
	addBranch("tau_ptWeightedDrIsolation");
	addBranch("tau_leadingTrackChi2");
	addBranch("tau_eRatio");
	addBranch("tau_gjAngleDiff");;


	setBranchType(kVectorUInt);
	addBranch("tau_numberOfIsolationChargedHadrCands");
	addBranch("tau_numberOfSignalChargedHadrCands");
	addBranch("tau_numNeutralHadronsSignalCone");
	addBranch("tau_numPhotonsSignalCone");
	addBranch("tau_numParticlesSignalCone");
	addBranch("tau_numChargedParticlesIsoCone");
	addBranch("tau_numNeutralHadronsIsoCone");
	addBranch("tau_numPhotonsIsoCone");
	addBranch("tau_numParticlesIsoCone");

	setBranchType(kVectorInt);
	addBranch("tau_mc_index");
	addBranch("tau_decayMode");
	addBranch("tau_charge");

	setBranchType(kVectorInt);
	addBranch("tau_isPFTau");
	addBranch("tau_hasSecondaryVertex");
	addBranch("tau_leadChargedHadrAvailable");


	setBranchType(kVectorFloat);
	addBranch("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017");
	addBranch("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017");
	addBranch("tau_byIsolationMVArun2017v2DBnewDMwLTraw2017");
	addBranch("tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byTightIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017");
	addBranch("tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017");

	addBranch("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017");
	addBranch("tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");
	addBranch("tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017");

	addBranch("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017");
	addBranch("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017");
	addBranch("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017");

	addBranch("tau_byIsolationMVArun2v1DBnewDMwLTraw2016");
	addBranch("tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016");
	addBranch("tau_byLooseIsolationMVArun2v1DBnewDMwLT2016");
	addBranch("tau_byMediumIsolationMVArun2v1DBnewDMwLT2016");
	addBranch("tau_byTightIsolationMVArun2v1DBnewDMwLT2016");
	addBranch("tau_byVTightIsolationMVArun2v1DBnewDMwLT2016");
	addBranch("tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016");

	addBranch("tau_byIsolationMVArun2v1DBoldDMwLTraw2016");
	addBranch("tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016");
	addBranch("tau_byLooseIsolationMVArun2v1DBoldDMwLT2016");
	addBranch("tau_byMediumIsolationMVArun2v1DBoldDMwLT2016");
	addBranch("tau_byTightIsolationMVArun2v1DBoldDMwLT2016");
	addBranch("tau_byVTightIsolationMVArun2v1DBoldDMwLT2016");
	addBranch("tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016");


}

// ------------ method called to for each event  ------------
void IIHEModuleTau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	Handle<View<pat::Tau> > tauCollection_ ;
	iEvent.getByToken( tauCollectionToken_, tauCollection_ );

	edm::Handle<View<reco::Vertex> > pvCollection_ ;
	iEvent.getByToken( vtxToken_ , pvCollection_);
	edm::Ptr<reco::Vertex> firstpvertex = pvCollection_->ptrAt( 0 );      

	unsigned int tau_n = 0 ;
	for ( unsigned int i = 0; i <tauCollection_->size(); ++i) {
		Ptr<pat::Tau> tauni = tauCollection_->ptrAt( i );
		if(tauni->pt() < ETThreshold_) continue ;
		tau_n++ ;
		store("tau_px"    , tauni->px()) ;
		store("tau_py"    , tauni->py()) ;
		store("tau_pz"    , tauni->pz()) ;
		store("tau_pt"    , tauni->pt()) ;
		store("tau_eta"   , tauni->eta()) ;
		store("tau_theta" , tauni->theta()) ;
		store("tau_phi"   , tauni->phi()) ;
		store("tau_energy", tauni->energy()) ;
		store("tau_mass"  , tauni->mass()) ;
		store("tau_dxy"   , tauni->dxy()) ;
		store("tau_dxy_error"         , tauni->dxy_error()) ;
		store("tau_ptLeadChargedCand" , tauni->ptLeadChargedCand()) ;
		store("tau_decayModeFinding"                           , tauni->tauID("decayModeFinding") );
		store("tau_decayModeFindingNewDMs"                     , tauni->tauID("decayModeFindingNewDMs") );
		store("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits" , tauni->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) ;
		store("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", tauni->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") );
		store("tau_byTightCombinedIsolationDeltaBetaCorr3Hits" , tauni->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") );
		store("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits"   , tauni->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") );
		store("tau_byIsolationMVArun2v1DBoldDMwLTraw"          , tauni->tauID("byIsolationMVArun2v1DBoldDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1DBoldDMwLT"       , tauni->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byLooseIsolationMVArun2v1DBoldDMwLT"        , tauni->tauID("byLooseIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byMediumIsolationMVArun2v1DBoldDMwLT"       , tauni->tauID("byMediumIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byTightIsolationMVArun2v1DBoldDMwLT"        , tauni->tauID("byTightIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byVTightIsolationMVArun2v1DBoldDMwLT"       , tauni->tauID("byVTightIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1DBoldDMwLT"      , tauni->tauID("byVVTightIsolationMVArun2v1DBoldDMwLT") );
		store("tau_byIsolationMVArun2v1DBnewDMwLTraw"          , tauni->tauID("byIsolationMVArun2v1DBnewDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1DBnewDMwLT"       , tauni->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byLooseIsolationMVArun2v1DBnewDMwLT"        , tauni->tauID("byLooseIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byMediumIsolationMVArun2v1DBnewDMwLT"       , tauni->tauID("byMediumIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byTightIsolationMVArun2v1DBnewDMwLT"        , tauni->tauID("byTightIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byVTightIsolationMVArun2v1DBnewDMwLT"       , tauni->tauID("byVTightIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1DBnewDMwLT"      , tauni->tauID("byVVTightIsolationMVArun2v1DBnewDMwLT") );
		store("tau_byIsolationMVArun2v1PWoldDMwLTraw"          , tauni->tauID("byIsolationMVArun2v1PWoldDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1PWoldDMwLT"       , tauni->tauID("byVLooseIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byLooseIsolationMVArun2v1PWoldDMwLT"        , tauni->tauID("byLooseIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byMediumIsolationMVArun2v1PWoldDMwLT"       , tauni->tauID("byMediumIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byTightIsolationMVArun2v1PWoldDMwLT"        , tauni->tauID("byTightIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byVTightIsolationMVArun2v1PWoldDMwLT"       , tauni->tauID("byVTightIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1PWoldDMwLT"      , tauni->tauID("byVVTightIsolationMVArun2v1PWoldDMwLT") );
		store("tau_byIsolationMVArun2v1PWnewDMwLTraw"          , tauni->tauID("byIsolationMVArun2v1PWnewDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1PWnewDMwLT"       , tauni->tauID("byVLooseIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byLooseIsolationMVArun2v1PWnewDMwLT"        , tauni->tauID("byLooseIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byMediumIsolationMVArun2v1PWnewDMwLT"       , tauni->tauID("byMediumIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byTightIsolationMVArun2v1PWnewDMwLT"        , tauni->tauID("byTightIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byVTightIsolationMVArun2v1PWnewDMwLT"       , tauni->tauID("byVTightIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1PWnewDMwLT"      , tauni->tauID("byVVTightIsolationMVArun2v1PWnewDMwLT") );
		store("tau_byIsolationMVArun2v1DBdR03oldDMwLTraw"      , tauni->tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT"   , tauni->tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT"    , tauni->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT"   , tauni->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT"    , tauni->tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT"   , tauni->tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT"  , tauni->tauID("byVVTightIsolationMVArun2v1DBdR03oldDMwLT") );
		store("tau_byIsolationMVArun2v1PWdR03oldDMwLTraw"      , tauni->tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw") );
		store("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT"   , tauni->tauID("byVLooseIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT"    , tauni->tauID("byLooseIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT"   , tauni->tauID("byMediumIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT"    , tauni->tauID("byTightIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT"   , tauni->tauID("byVTightIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT"  , tauni->tauID("byVVTightIsolationMVArun2v1PWdR03oldDMwLT") );
		store("tau_againstMuonLoose3"          , tauni->tauID("againstMuonLoose3") );
		store("tau_againstMuonTight3"          , tauni->tauID("againstMuonTight3") );
		store("tau_againstElectronMVA6Raw"     , tauni->tauID("againstElectronMVA6Raw") );
		store("tau_againstElectronMVA6category", tauni->tauID("againstElectronMVA6category") );
		store("tau_againstElectronVLooseMVA6"  , tauni->tauID("againstElectronVLooseMVA6") );
		store("tau_againstElectronLooseMVA6"   , tauni->tauID("againstElectronLooseMVA6") );
		store("tau_againstElectronMediumMVA6"  , tauni->tauID("againstElectronMediumMVA6") );
		store("tau_againstElectronTightMVA6"   , tauni->tauID("againstElectronTightMVA6") );
		store("tau_againstElectronVTightMVA6"  , tauni->tauID("againstElectronVTightMVA6") );
		//		std::cout<<"tauID:"<< tauni->tauID("byIsolationMVArun2017v1DBoldDMwLTraw2017") << std::endl;
		/*		store("tau_byVVLooseIsolationMVArun2v1DBoldDMwLTNew"    , tauni->tauID("byVVLooseIsolationMVArun2v1DBoldDMwLTNew") );
				store("tau_byVLooseIsolationMVArun2v1DBoldDMwLTNew"    , tauni->tauID("byVLooseIsolationMVArun2v1DBoldDMwLTNew") );
				store("tau_byLooseIsolationMVArun2v1DBoldDMwLTNew"     , tauni->tauID("byLooseIsolationMVArun2v1DBoldDMwLTNew"));
				store("tau_byMediumIsolationMVArun2v1DBoldDMwLTNew"     , tauni->tauID("byMediumIsolationMVArun2v1DBoldDMwLTNew"));
				store("tau_byTightIsolationMVArun2v1DBoldDMwLTNew"     , tauni->tauID("byTightIsolationMVArun2v1DBoldDMwLTNew"));
				store("tau_byVTightIsolationMVArun2v1DBoldDMwLTNew"    , tauni->tauID("byVTightIsolationMVArun2v1DBoldDMwLTNew"));
				store("tau_byVVTightIsolationMVArun2v1DBoldDMwLTNew"   , tauni->tauID("byVVTightIsolationMVArun2v1DBoldDMwLTNew"));
				store("tau_byIsolationMVArun2v1DBoldDMwLTrawNew"       , tauni->tauID("byIsolationMVArun2v1DBoldDMwLTrawNew"));
				store("tau_againstElectronMVA6RawNew"  , tauni->tauID("againstElectronMVA6RawNew"));
				*/		store("tau_chargedIsoPtSum" , tauni->tauID("chargedIsoPtSum"));
		store("tau_neutralIsoPtSum", tauni->tauID("neutralIsoPtSum"));
		store("tau_puCorrPtSum", tauni->tauID("puCorrPtSum"));
		store("tau_footprintCorrection", tauni->tauID("footprintCorrection"));
		store("tau_neutralIsoPtSumWeight", tauni->tauID("neutralIsoPtSumWeight"));
		store("tau_photonPtSumOutsideSignalCone", tauni->tauID("photonPtSumOutsideSignalCone"));
		store("tau_byPhotonPtSumOutsideSignalCone", tauni->tauID("byPhotonPtSumOutsideSignalCone"));
		store("tau_footprintCorrectiondR03", tauni->tauID("footprintCorrectiondR03"));
		store("tau_chargedIsoPtSumdR03", tauni->tauID("chargedIsoPtSumdR03"));
		store("tau_neutralIsoPtSumWeightdR03", tauni->tauID("neutralIsoPtSumWeightdR03"));
		store("tau_neutralIsoPtSumdR03", tauni->tauID("neutralIsoPtSumdR03"));
		store("tau_photonPtSumOutsideSignalConedR03", tauni->tauID("photonPtSumOutsideSignalConedR03"));

		//////

		store("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017", tauni->tauID("byIsolationMVArun2017v1DBoldDMwLTraw2017"));
		store("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byVVLooseIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byVLooseIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byLooseIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byMediumIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byTightIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byVTightIsolationMVArun2017v1DBoldDMwLT2017"));
		store("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017", tauni->tauID("byVVTightIsolationMVArun2017v1DBoldDMwLT2017"));


		store("tau_byIsolationMVArun2017v2DBnewDMwLTraw2017", tauni->tauID("byIsolationMVArun2017v2DBnewDMwLTraw2017"));
		store("tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byVVLooseIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byVLooseIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byLooseIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byMediumIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byTightIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byTightIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byVTightIsolationMVArun2017v2DBnewDMwLT2017"));
		store("tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017", tauni->tauID("byVVTightIsolationMVArun2017v2DBnewDMwLT2017"));

		store("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", tauni->tauID("byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017"));
		store("tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		store("tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", tauni->tauID("byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));

		store("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017", tauni->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
		store("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
		store("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017", tauni->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));

		store("tau_byIsolationMVArun2v1DBnewDMwLTraw2016", tauni->tauID("byIsolationMVArun2v1DBnewDMwLTraw2016"));
		store("tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT2016"));
		store("tau_byLooseIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byLooseIsolationMVArun2v1DBnewDMwLT2016"));
		store("tau_byMediumIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byMediumIsolationMVArun2v1DBnewDMwLT2016"));
		store("tau_byTightIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byTightIsolationMVArun2v1DBnewDMwLT2016"));
		store("tau_byVTightIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byVTightIsolationMVArun2v1DBnewDMwLT2016"));
		store("tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016", tauni->tauID("byVVTightIsolationMVArun2v1DBnewDMwLT2016"));

		store("tau_byIsolationMVArun2v1DBoldDMwLTraw2016", tauni->tauID("byIsolationMVArun2v1DBoldDMwLTraw2016"));
		store("tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT2016"));
		store("tau_byLooseIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byLooseIsolationMVArun2v1DBoldDMwLT2016"));
		store("tau_byMediumIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byMediumIsolationMVArun2v1DBoldDMwLT2016"));
		store("tau_byTightIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byTightIsolationMVArun2v1DBoldDMwLT2016"));
		store("tau_byVTightIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byVTightIsolationMVArun2v1DBoldDMwLT2016"));
		store("tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016", tauni->tauID("byVVTightIsolationMVArun2v1DBoldDMwLT2016"));

		store("tau_decayMode"         , tauni->decayMode()) ;
		store("tau_charge"            , tauni->charge()) ;
		store("tau_isPFTau"           , int(tauni->isPFTau())) ;
		store("tau_hasSecondaryVertex", int(tauni->hasSecondaryVertex())) ;

		store("tau_PFChargedHadIso", tauni->chargedHadronIso());
		store("tau_PFNeutralHadIso", tauni->neutralHadronIso());
		store("tau_PFPhotonIso", tauni->photonIso());

		store("tau_numberOfIsolationChargedHadrCands" , tauni->isolationChargedHadrCands().size());
		store("tau_numberOfSignalChargedHadrCands"  , tauni->signalChargedHadrCands().size());

		store("tau_numNeutralHadronsSignalCone" , tauni->signalNeutrHadrCands().size());
		store("tau_numPhotonsSignalCone", tauni->signalGammaCands().size());
		store("tau_numParticlesSignalCone", tauni->signalCands().size());
		store("tau_numChargedParticlesIsoCone", tauni->isolationChargedHadrCands().size());
		store("tau_numNeutralHadronsIsoCone", tauni->isolationNeutrHadrCands().size());
		store("tau_numPhotonsIsoCone" ,tauni->isolationGammaCands().size());
		store("tau_numParticlesIsoCone", tauni->isolationCands().size());
		store("tau_leadChargedParticlePt",tauni->leadCand()->pt());
		store("tau_trackRefPt", (tauni->leadChargedHadrCand().isNonnull() ? tauni->leadChargedHadrCand()->pt() : 0.));


		// new variables for TauPOG study
		bool leadChargedHadrAvailable = false;
		float leadingTrackChi2  = -999.;
		float gjAngleDiff = -999;

		if ( tauni->leadChargedHadrCand().isNonnull() ) {
			leadChargedHadrAvailable = true; 
			leadingTrackChi2 = tauni->leadingTrackNormChi2(); 
		}


		int tauDecayMode = tauni->decayMode();            
		float decayDistX = tauni->flightLength().x();
		float decayDistY = tauni->flightLength().y();
		float decayDistZ = tauni->flightLength().z();
		float decayDistMag = std::sqrt(decayDistX*decayDistX + decayDistY*decayDistY + decayDistZ*decayDistZ);

		float nPhoton = n_photons_total(*tauni);
		float ptWeightedDetaStrip = pt_weighted_deta_strip(*tauni, tauDecayMode);
		float ptWeightedDphiStrip = pt_weighted_dphi_strip(*tauni, tauDecayMode);
		float ptWeightedDrSignal = pt_weighted_dr_signal(*tauni, tauDecayMode);
		float ptWeightedDrIsolation = pt_weighted_dr_iso(*tauni, tauDecayMode);

		float eRatio = eratio(*tauni);

		// Difference between measured and maximally allowed Gottfried-Jackson angle
		if ( tauDecayMode == 10 ) {
			double mTau = 1.77682;
			double mAOne = tauni->p4().M();
			double pAOneMag = tauni->p();
			double argumentThetaGJmax = (std::pow(mTau,2) - std::pow(mAOne,2) ) / ( 2 * mTau * pAOneMag );
			double argumentThetaGJmeasured = ( tauni->p4().px() * decayDistX + tauni->p4().py() * decayDistY + tauni->p4().pz() * decayDistZ ) / ( pAOneMag * decayDistMag );
			if ( std::abs(argumentThetaGJmax) <= 1. && std::abs(argumentThetaGJmeasured) <= 1. ) {
				double thetaGJmax = std::asin( argumentThetaGJmax );
				double thetaGJmeasured = std::acos( argumentThetaGJmeasured );
				gjAngleDiff = thetaGJmeasured - thetaGJmax;
			}
		}


		store("tau_dxy_Sig", tauni->dxy_Sig());
		store("tau_flightLengthSig", tauni->flightLengthSig());
		store("tau_ip3d", tauni->ip3d());
		store("tau_ip3d_Sig", tauni->ip3d_Sig());

		store("tau_leadChargedHadrAvailable", int(leadChargedHadrAvailable));
		store("tau_decayDistX", decayDistX);
		store("tau_decayDistY", decayDistY);;
		store("tau_decayDistZ",  decayDistZ);
		store("tau_decayDistMag", decayDistMag);
		store("tau_nPhoton", nPhoton);
		store("tau_ptWeightedDetaStrip", ptWeightedDetaStrip);
		store("tau_ptWeightedDphiStrip", ptWeightedDphiStrip);
		store("tau_ptWeightedDrSignal", ptWeightedDrSignal);
		store("tau_ptWeightedDrIsolation", ptWeightedDrIsolation);
		store("tau_leadingTrackChi2", leadingTrackChi2);
		store("tau_eRatio", eRatio);
		store("tau_gjAngleDiff", gjAngleDiff);;

		float dxy = 999.;
		float dz  = 999.;
		if (pvCollection_->size()>0) {

			pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tauni->leadChargedHadrCand().get());
			dz=packedLeadTauCand->dz();
			dxy=packedLeadTauCand->dxy();
		}

		store("tau_lead_dxy", dxy);
		store("tau_lead_dz", dz);

		// Now apply truth matching.
		int index = MCTruth_matchEtaPhi_getIndex(tauni->eta(), tauni->phi()) ;
		if(index>=0){
			const MCTruthObject* MCTruth = MCTruth_getRecordByIndex(index) ;
			store("tau_mc_bestDR", deltaR(tauni->eta(), tauni->phi(), MCTruth->eta(), MCTruth->phi())) ;
			store("tau_mc_index" , index) ;
			store("tau_mc_ERatio", tauni->energy()/MCTruth->energy()) ;
		}
		else{
			store("tau_mc_bestDR", 999.0) ;
			store("tau_mc_index" ,    -1) ;
			store("tau_mc_ERatio", 999.0) ;
		}

	}
	store("tau_n", tau_n) ;
}
void IIHEModuleTau::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleTau::beginEvent(){}
void IIHEModuleTau::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTau::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleTau);
