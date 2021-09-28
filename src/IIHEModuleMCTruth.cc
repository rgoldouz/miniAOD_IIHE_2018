#include "UserCode/IIHETree/interface/IIHEModuleMCTruth.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleMCTruth::IIHEModuleMCTruth(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
	pt_threshold_            = iConfig.getUntrackedParameter<double>("MCTruth_ptThreshold"            ) ;
	m_threshold_             = iConfig.getUntrackedParameter<double>("MCTruth_mThreshold"             ) ;
	DeltaROverlapThreshold_  = iConfig.getUntrackedParameter<double>("MCTruth_DeltaROverlapThreshold" ) ;
	puInfoSrc_               = iConfig.getParameter<edm::InputTag>("PileUpSummaryInfo") ;
	generatorLabel_ = iC.consumes<GenEventInfoProduct> (iConfig.getParameter<InputTag>("generatorLabel"));
	lheEventLabel_ = iC.consumes<LHEEventProduct> (iConfig.getParameter<InputTag>("LHELabel"));
	puCollection_ = iC.consumes<vector<PileupSummaryInfo> > (puInfoSrc_);
	genParticlesCollection_ = iC.consumes<vector<reco::GenParticle> > (iConfig.getParameter<InputTag>("genParticleSrc"));
	genJetsSrc_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>( "genJetsCollection" ));
        genJetsSrcAK8_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>( "genJetsCollectionAK8" ));
}
IIHEModuleMCTruth::~IIHEModuleMCTruth(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMCTruth::beginJob(){
	std::vector<int> MCPdgIdsToSave ;
	MCPdgIdsToSave.push_back(11) ; // Electron
	MCPdgIdsToSave.push_back(12) ; // Neutrino electron
	MCPdgIdsToSave.push_back(13) ; // Muon
	MCPdgIdsToSave.push_back(14) ; // Neutrino muon
	MCPdgIdsToSave.push_back(15) ; // Tau
	MCPdgIdsToSave.push_back(16) ; // Neutrino tau
	MCPdgIdsToSave.push_back(23) ; // Z boson
	MCPdgIdsToSave.push_back(24) ; // W boson
	MCPdgIdsToSave.push_back(25) ; // BEH boson
	MCPdgIdsToSave.push_back(32) ; // Z'  boson
	MCPdgIdsToSave.push_back(33) ; // Z'' boson
	MCPdgIdsToSave.push_back(34) ; // W'  boson
        MCPdgIdsToSave.push_back(37) ;
	MCPdgIdsToSave.push_back(1000016) ; // tau-sneutrino 
        MCPdgIdsToSave.push_back(600) ; // Excited top
	addToMCTruthWhitelist(MCPdgIdsToSave) ;
	addBranch("mc_n", kUInt) ;
        addBranch("mc_nMEPartons", kInt) ;
        addBranch("mc_nMEPartonsFiltered", kInt) ;
        addBranch("mc_DJRValues", kVectorFloat) ;
	addBranch("mc_weight", kFloat) ;
	addBranch("mc_w_sign", kFloat) ;
	addBranch("mc_id_first", kInt) ;
	addBranch("mc_id_second", kInt) ;
	addBranch("mc_x_first", kFloat) ;
	addBranch("mc_x_second", kFloat) ;
	addBranch("mc_xPDF_first", kFloat) ;
	addBranch("mc_xPDF_second", kFloat) ;
	addBranch("mc_scalePDF", kFloat) ;
	setBranchType(kVectorInt) ;
	addBranch("mc_index") ;
	addBranch("mc_pdgId") ;
	addBranch("mc_charge") ;
	addBranch("mc_status") ;
	addBranch("mc_status_flags");
	// new additions
//	addBranch("mc_status_tau_flags");
//	addBranch("mc_tau_charge");
//	addBranch("mc_tau_pdgId");
//	addBranch("mc_tau_decay");
//	addBranch("mc_tau_had_status");
//	addBranch("mc_tau_had_charge");
//	addBranch("mc_tau_had_pdgId");

	setBranchType(kVectorFloat) ;
	addBranch("mc_mass") ;
	addBranch("mc_px") ;
	addBranch("mc_py") ;
	addBranch("mc_pz") ;
	addBranch("mc_pt") ;
	addBranch("mc_eta") ;
	addBranch("mc_phi") ;
	addBranch("mc_energy") ;
	// new additions

//	addBranch("mc_tau_pt");
//	addBranch("mc_tau_eta");
//	addBranch("mc_tau_phi");
//	addBranch("mc_tau_energy");
//
//	addBranch("mc_tau_had_pt");
//	addBranch("mc_tau_had_eta");
//	addBranch("mc_tau_had_phi");
//	addBranch("mc_tau_had_energy");
	////

	setBranchType(kVectorUInt) ;
	addBranch("mc_numberOfDaughters") ;
	addBranch("mc_numberOfMothers"  ) ;
	setBranchType(kVectorVectorInt) ;
	addBranch("mc_mother_index") ;
	addBranch("mc_mother_pdgId") ;
	setBranchType(kVectorVectorFloat) ;
	addBranch("mc_mother_px"    ) ;
	addBranch("mc_mother_py"    ) ;
	addBranch("mc_mother_pz"    ) ;
	addBranch("mc_mother_pt"    ) ;
	addBranch("mc_mother_eta"   ) ;
	addBranch("mc_mother_phi"   ) ;
	addBranch("mc_mother_energy") ;
	addBranch("mc_mother_mass"  ) ;

	setBranchType(kInt) ;
	addBranch("mc_trueNumInteractions") ;
	addBranch("mc_PU_NumInteractions" ) ;

	addValueToMetaTree("MCTruth_ptThreshold"           , pt_threshold_          ) ;
	addValueToMetaTree("MCTruth_mThreshold"            , m_threshold_           ) ;
	addValueToMetaTree("MCTruth_DeltaROverlapThreshold", DeltaROverlapThreshold_) ;

	setBranchType(kVectorFloat) ;
	addBranch("genjet_pt") ;
	addBranch("genjet_eta") ;
	addBranch("genjet_phi") ;
	addBranch("genjet_energy") ;
        addBranch("genjetAK8_pt") ;
        addBranch("genjetAK8_eta") ;
        addBranch("genjetAK8_phi") ;
        addBranch("genjetAK8_energy") ;
	setBranchType(kVectorFloat) ;
	addBranch("gen_weight_sys");

	nEventsWeighted_ = 0.0 ;
	pileupDist_ = new TH1F("pileupDist","pileup distribution",120,0,120);
}

// ------------ method called to for each event  ------------
void IIHEModuleMCTruth::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	edm::Handle<std::vector<reco::GenJet> > genJets;
	iEvent.getByToken(genJetsSrc_, genJets);

	for (size_t j = 0; j < genJets->size();++j){
		store("genjet_pt", genJets->at(j).pt()) ;
		store("genjet_eta", genJets->at(j).eta()) ;
		store("genjet_phi", genJets->at(j).phi()) ;
		store("genjet_energy", genJets->at(j).energy()) ;
	}

        edm::Handle<std::vector<reco::GenJet> > genJetsAK8;
        iEvent.getByToken(genJetsSrcAK8_, genJetsAK8);

        for (size_t j = 0; j < genJetsAK8->size();++j){
                store("genjetAK8_pt", genJetsAK8->at(j).pt()) ;
                store("genjetAK8_eta", genJetsAK8->at(j).eta()) ;
                store("genjetAK8_phi", genJetsAK8->at(j).phi()) ;
                store("genjetAK8_energy", genJetsAK8->at(j).energy()) ;
        }

	edm::Handle<GenEventInfoProduct> genEventInfoHandle;
	iEvent.getByToken(generatorLabel_, genEventInfoHandle);
        store("mc_nMEPartons"            ,genEventInfoHandle->nMEPartons());
        store("mc_nMEPartonsFiltered"    ,genEventInfoHandle->nMEPartonsFiltered());
        for (auto djr_val: genEventInfoHandle->DJRValues()) {
          store("mc_DJRValues"            , djr_val);
        }

	float weight = genEventInfoHandle->weight() ;
	float w_sign = (weight>=0) ? 1 : -1 ;
	store("mc_weight"                  ,weight);
	store("mc_w_sign"             , w_sign) ;
	nEventsWeighted_ += w_sign ;

	store("mc_id_first" , genEventInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
	store("mc_id_second" , genEventInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
	store("mc_x_first", genEventInfoHandle->pdf()->x.first);     // x value of parton #1
	store("mc_x_second" , genEventInfoHandle->pdf()->x.second);    // x value of parton #2
	store("mc_xPDF_first" , genEventInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
	store("mc_xPDF_second" , genEventInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
	store("mc_scalePDF" , genEventInfoHandle->pdf()->scalePDF);    // scale of the hard interaction

	/*
	   First two weights (weightID= 0 and 1) correspond to central ME weight value and replica.
	   The remaining 12 values (weightIDs = 2 to 13) correspond to the PS weights in the following order (ISR up, FSR up, ISR down, FSR down) x 3 sets, i.e.:
	   2 = isrRedHi isr:muRfac=0.707, 3 = fsrRedHi fsr:muRfac=0.707, 4 = isrRedLo isr:muRfac=1.414, 5 = fsrRedLo fsr:muRfac=1.414, 
	   6 = isrDefHi isr:muRfac=0.5, 7 = fsrDefHi fsr:muRfac=0.5,  8 = isrDefLo isr:muRfac=2.0,   9 = fsrDefLo fsr:muRfac=2.0, 
	   10 = isrConHi isr:muRfac=0.25, 11 = fsrConHi fsr:muRfac=0.25, 12 = isrConLo isr:muRfac=4.0, 13 = fsrConLo fsr:muRfac=4.0 
	   */

	if(!(genEventInfoHandle->weights().empty())) {
		for (unsigned i = 0; i < genEventInfoHandle->weights().size(); ++i) {
			if(sumofgenWeights_.size() != genEventInfoHandle->weights().size()) {
				sumofgenWeights_.push_back(0);
			}
			store("gen_weight_sys",(float) genEventInfoHandle->weights().at(i));
			sumofgenWeights_[i] += (float) genEventInfoHandle->weights().at(i);
		}
	} else {
		store("gen_weight_sys", 1);
	}
	// Fill pile-up related informations
	// --------------------------------
	edm::Handle<std::vector< PileupSummaryInfo > >  puInfo ;
	iEvent.getByToken(puCollection_, puInfo) ;
	int trueNumInteractions = -1 ;
	int PU_NumInteractions  = -1 ;
	if(puInfo.isValid()){
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for(PVI = puInfo->begin() ; PVI != puInfo->end() ; ++PVI){
			int BX = PVI->getBunchCrossing() ;
			if(BX==0){ // "0" is the in-time crossing, negative values are the early crossings, positive are late
				pileupDist_->Fill(PVI->getTrueNumInteractions());
				trueNumInteractions = PVI->getTrueNumInteractions() ;
				PU_NumInteractions  = PVI->getPU_NumInteractions() ;
			}
		}
	}

	Handle<GenParticleCollection> pGenParticles ;
	iEvent.getByToken(genParticlesCollection_, pGenParticles) ;
	GenParticleCollection genParticles(pGenParticles->begin(),pGenParticles->end()) ;
	// tau variables

	/*unsigned int NGenSel = genParticles.size();
	for (unsigned int iGen = 0; iGen < NGenSel; iGen++)
	{
		const reco::Candidate* genP = &genParticles.at(iGen);
		const GenParticle* genPClone = (GenParticle*) &genParticles.at(iGen);
		if (fabs(genP->pdgId()) != 15) continue;
		store("mc_tau_pt", genP->pt());
		store("mc_tau_eta", genP->eta());
		store("mc_tau_phi", genP->phi());
		store("mc_tau_energy",genP->energy());
		store("mc_tau_charge",genP->charge());
		store("mc_tau_pdgId",genP->pdgId());
		store("mc_tau_decay",(GetTauDecay(genP)));
		if(GetTauDecay(genP) == 2) {
			store("mc_tau_had_pt",(GetTauHad(genP).pt()));
			store("mc_tau_had_eta",(GetTauHad(genP).eta()));
			store("mc_tau_had_phi",(GetTauHad(genP).phi()));
			store("mc_tau_had_energy",(GetTauHad(genP).energy()));
			store("mc_tau_had_status",(GetTauHad(genP).status()));
			store("mc_tau_had_charge",(GetTauHad(genP).charge()));
			store("mc_tau_had_pdgId",(GetTauHad(genP).pdgId()));

		}
		int flags = 0;
		if (genPClone->statusFlags().isPrompt())                  flags |= (1 << 0);
		if (genPClone->statusFlags().isDecayedLeptonHadron())     flags |= (1 << 1);
		if (genPClone->statusFlags().isTauDecayProduct())         flags |= (1 << 2);
		if (genPClone->statusFlags().isPromptTauDecayProduct())   flags |= (1 << 3);
		if (genPClone->statusFlags().isDirectTauDecayProduct())   flags |= (1 << 4);
		if (genPClone->statusFlags().isDirectPromptTauDecayProduct())       flags |= (1 << 5);
		if (genPClone->statusFlags().isDirectHadronDecayProduct())          flags |= (1 << 6);
		if (genPClone->statusFlags().isHardProcess())             flags |= (1 << 7);
		if (genPClone->statusFlags().fromHardProcess())           flags |= (1 << 8);
		if (genPClone->statusFlags().isHardProcessTauDecayProduct())       flags |= (1 << 9);
		if (genPClone->statusFlags().isDirectHardProcessTauDecayProduct()) flags |= (1 << 10);
		if (genPClone->statusFlags().fromHardProcessBeforeFSR())  flags |= (1 << 11);
		if (genPClone->statusFlags().isLastCopyBeforeFSR())       flags |= (1 << 12);

		store("mc_status_tau_flags" , flags);
	}
*/
	// These variables are used to match up mothers to daughters at the end.
	int counter = 0 ;

	MCTruthRecord_.clear() ;

	MCTruthObject* MCTruth ;
	//  const Candidate* parent ;
	const Candidate* child ;

	for(GenParticleCollection::const_iterator mc_iter = genParticles.begin() ; mc_iter!=genParticles.end() ; ++mc_iter){
		int pdgId = mc_iter->pdgId() ;
		float pt  = mc_iter->pt()    ;
                int status = mc_iter->status();
		// First check the whitelist.
		bool whitelist_accept = false ;
		for(unsigned int i=0 ; i<whitelist_.size() ; ++i){
			if(abs(pdgId)==abs(whitelist_.at(i))){
				whitelist_accept = true ;
				break ;
			}
		}


		// Remove objects with zero pT.
		bool nonZeroPt_accept = (pt>1e-3) ;

		// Now check the thresholds.
		bool thresholds_accept = false ;
		float pt_threshold = (pdgId==21 || abs(pdgId)<5) ? 10 : pt_threshold_ ;
		if(             pt>pt_threshold ) thresholds_accept = true  ;
		if(mc_iter->mass()> m_threshold_) thresholds_accept = true  ;

		// Now combine them all.
		bool accept = (whitelist_accept && thresholds_accept && nonZeroPt_accept) ;
                if(status == 3 || (status > 20 && status < 30)) accept=true; //keep matrix element summary
                if(pdgId == 22 && status == 1 && (pt > 10 || mc_iter->isPromptFinalState())) accept=true; // keep gamma above 10 GeV (or all prompt) 
                if(mc_iter->isHardProcess() ||  mc_iter->fromHardProcessDecayed()  || mc_iter->fromHardProcessFinalState() || (mc_iter->statusFlags().fromHardProcess() && mc_iter->statusFlags().isLastCopy())) accept=true; //keep event summary based on status flags
//                if(status > 70 && status < 80 && pt > 15) accept=true; //keep high pt partons right before hadronization
		if(false==accept) continue ;
		// Now go up the ancestry until we find the real parent
		//    parent = mc_iter->mother() ;
		child  = mc_iter->clone()  ;
		//    while(parent->pdgId()==pdgId){
		//      child  = parent ;
		//      parent = parent->mother() ;
		//    }

		// Create a truth record instance.
		MCTruth = new MCTruthObject((reco::Candidate*)&*mc_iter) ;

		// Add all the mothers
		for(unsigned int mother_iter=0 ; mother_iter<child->numberOfMothers() ; ++mother_iter){
			MCTruth->addMother(child->mother(mother_iter)) ;
		}

		// Finally check to see if this overlaps with an existing truth particle.
		bool overlap = false ;
		for(unsigned int i=0 ; i<MCTruthRecord_.size() ; ++i){
			const reco::Candidate* comp = MCTruthRecord_.at(i)->getCandidate() ;
			float DR = deltaR(comp->eta(),comp->phi(),MCTruth->getCandidate()->eta(),MCTruth->getCandidate()->phi()) ;
			if(DR<DeltaROverlapThreshold_){
				overlap = true ;
				break ;
			}
		}
		if(true==overlap) continue ;

		// Then push back the MC truth information
		MCTruthRecord_.push_back(MCTruth) ;
		counter++ ;
	}
	for(unsigned int i=0 ; i<MCTruthRecord_.size() ; ++i){
                MCTruthObject* ob = MCTruthRecord_.at(i) ;
		std::vector<int  > mc_mother_index ;
		std::vector<int  > mc_mother_pdgId ;
		std::vector<float> mc_mother_px ;
		std::vector<float> mc_mother_py ;
		std::vector<float> mc_mother_pz ;
		std::vector<float> mc_mother_pt ;
		std::vector<float> mc_mother_eta ;
		std::vector<float> mc_mother_phi ;
		std::vector<float> mc_mother_energy ;
		std::vector<float> mc_mother_mass ;
		for(unsigned int j=0 ; j<ob->nMothers() ; ++j){
//                        if (j>2) continue;
			const reco::Candidate* mother = ob->getMother(j) ;
			if(mother){
				int mother_index_tmp = ob->matchMother(MCTruthRecord_, j) ;
				mc_mother_index .push_back(mother_index_tmp) ;
				mc_mother_pdgId .push_back(mother->pdgId() ) ;
				mc_mother_px    .push_back(mother->px()    ) ;
				mc_mother_py    .push_back(mother->py()    ) ;
				mc_mother_pz    .push_back(mother->pz()    ) ;
				mc_mother_pt    .push_back(mother->pt()    ) ;
				mc_mother_eta   .push_back(mother->eta()   ) ;
				mc_mother_phi   .push_back(mother->phi()   ) ;
				mc_mother_energy.push_back(mother->energy()) ;
				mc_mother_mass  .push_back(mother->mass()  ) ;
			}
		}
		if(mc_mother_index.size()==0) mc_mother_index.push_back(0) ;
		store("mc_mother_index" , mc_mother_index ) ;
		store("mc_mother_pdgId" , mc_mother_pdgId ) ;
		store("mc_mother_px"    , mc_mother_px    ) ;
		store("mc_mother_py"    , mc_mother_py    ) ;
		store("mc_mother_pz"    , mc_mother_pz    ) ;
		store("mc_mother_pt"    , mc_mother_pt    ) ;
		store("mc_mother_eta"   , mc_mother_eta   ) ;
		store("mc_mother_phi"   , mc_mother_phi   ) ;
		store("mc_mother_energy", mc_mother_energy) ;
		store("mc_mother_mass"  , mc_mother_mass  ) ;

		store("mc_numberOfDaughters", (unsigned int)(ob->getCandidate()->numberOfDaughters())) ;
		store("mc_numberOfMothers"  , (unsigned int)(ob->nMothers())) ;

		store("mc_px"     , ob->getCandidate()->px()    ) ;
		store("mc_py"     , ob->getCandidate()->py()    ) ;
		store("mc_pz"     , ob->getCandidate()->pz()    ) ;
		store("mc_pt"     , ob->getCandidate()->pt()    ) ;
		store("mc_eta"    , ob->getCandidate()->eta()   ) ;
		store("mc_phi"    , ob->getCandidate()->phi()   ) ;
		store("mc_energy" , ob->getCandidate()->energy()) ;
		store("mc_mass"   , ob->getCandidate()->mass()  ) ;

		store("mc_index"  , i ) ;
		store("mc_pdgId"  , ob->getCandidate()->pdgId() ) ;
		store("mc_charge" , ob->getCandidate()->charge()) ;
		store("mc_status" , ob->getCandidate()->status()) ;
		// storing flags of a gen particle
		int flags = 0;
		const reco::GenParticle* p = (GenParticle*) ob->getCandidate();

		if (p->statusFlags().isPrompt())                  flags |= (1 << 0);
		if (p->statusFlags().isDecayedLeptonHadron())     flags |= (1 << 1);
		if (p->statusFlags().isTauDecayProduct())         flags |= (1 << 2);
		if (p->statusFlags().isPromptTauDecayProduct())   flags |= (1 << 3);
		if (p->statusFlags().isDirectTauDecayProduct())   flags |= (1 << 4);
		if (p->statusFlags().isDirectPromptTauDecayProduct())       flags |= (1 << 5);
		if (p->statusFlags().isDirectHadronDecayProduct())          flags |= (1 << 6);
		if (p->statusFlags().isHardProcess())             flags |= (1 << 7);
		if (p->statusFlags().fromHardProcess())           flags |= (1 << 8);
		if (p->statusFlags().isHardProcessTauDecayProduct())       flags |= (1 << 9);
		if (p->statusFlags().isDirectHardProcessTauDecayProduct()) flags |= (1 << 10);
		if (p->statusFlags().fromHardProcessBeforeFSR())  flags |= (1 << 11);
		if (p->statusFlags().isLastCopyBeforeFSR())       flags |= (1 << 12);

		store("mc_status_flags" , flags);
	}

	store("mc_trueNumInteractions", trueNumInteractions) ;
	store("mc_PU_NumInteractions" , PU_NumInteractions ) ;

	store("mc_n", (unsigned int)(MCTruthRecord_.size())) ;
}

int IIHEModuleMCTruth::matchEtaPhi_getIndex(float eta, float phi){
	float bestDR = 1e6 ;
	int bestIndex = -1 ;
	for(unsigned int i=0 ; i<MCTruthRecord_.size() ; ++i){
		MCTruthObject* MCTruth = MCTruthRecord_.at(i) ;
		float DR = deltaR(MCTruth->eta(),MCTruth->phi(),eta,phi) ;
		if(DR<bestDR){
			bestDR = DR ;
			bestIndex = i ;
		}
	}
	return bestIndex ;
}

const MCTruthObject* IIHEModuleMCTruth::matchEtaPhi(float eta, float phi){
	int index = matchEtaPhi_getIndex(eta, phi) ;
	return getRecordByIndex(index) ;
}

const MCTruthObject* IIHEModuleMCTruth::getRecordByIndex(int index){
	if(index<0) return 0 ;
	return MCTruthRecord_.at(index) ;
}

int IIHEModuleMCTruth::GetTauDecay (const reco::Candidate* part)
{
	if (abs(part->pdgId()) != 15) return -1; // only on taus
	int decay = -1;
	int nele = 0;
	int nmu = 0;
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter(iDau);
		int dauId = abs(Dau->pdgId());
		if (dauId == 11) nele++;
		if (dauId == 13) nmu++;
	}

	if (nmu == 1 && nele == 0) decay = 0;
	if (nmu == 0 && nele == 1) decay = 1;
	if (nmu == 0 && nele == 0) decay = 2;

	return decay; // -1 if strange things happen
}
int IIHEModuleMCTruth::GetTauDecay (const reco::GenParticle& part)
{
	if ( !(part.statusFlags().isLastCopy()) ) return -1; // only for last copies
	const reco::Candidate* p = &part;
	return GetTauDecay(p);
}

const reco::Candidate* IIHEModuleMCTruth::IsFromID (const reco::Candidate* part, int targetPDGId)
{
	if (abs(part->pdgId()) == targetPDGId){
		if(abs(part->pdgId()) == 5) return GetFirstCopy(part);
		else return part;
	}

	for (unsigned int i = 0; i < part->numberOfMothers(); i++)
	{
		const reco::Candidate* matchMoth = IsFromID(part->mother(i), targetPDGId);
		if ( matchMoth != NULL) return matchMoth;
	}
	return NULL;

}

const reco::Candidate* IIHEModuleMCTruth::GetFirstCopy (const reco::Candidate* part)
{
	int cloneInd = -1;
	int id = part->pdgId();
	for (unsigned int iMot = 0; iMot < part->numberOfMothers(); iMot++)
	{
		const reco::Candidate * Mot = part->mother( iMot );
		if (id == Mot->pdgId())
		{
			cloneInd = iMot;
			break;
		}
	}

	if (cloneInd == -1) return part;
	else return (GetFirstCopy (part->mother(cloneInd)));

}
int IIHEModuleMCTruth::GetIndexInOutput (const reco::Candidate* part, std::vector<const reco::Candidate *> cands)
{
	int index = -1;
	std::vector<const reco::Candidate *>::const_iterator found = find(cands.begin(), cands.end(), part);
	if(found != cands.end()) index = found - cands.begin();
	return index;

}

reco::GenParticle IIHEModuleMCTruth::GetTauHad (const reco::Candidate* part)
{
	if (abs(part->pdgId()) != 15)
	{
		reco::GenParticle fakeTauH = reco::GenParticle (0, reco::Candidate::LorentzVector(0.,0.,0.,0.), reco::Candidate::Point (0.,0.,0.), -999999, 0, true);
		std::cout << "Warning: building had tau from a particle with pdgId != 15 --> dummy entry returned" << std::endl;
		return fakeTauH;
	}

	reco::Candidate::LorentzVector p4Had (0,0,0,0);
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter( iDau );
		int dauId = abs(Dau->pdgId());
		if (dauId != 12 && dauId != 14 && dauId != 16) // no neutrinos
			p4Had += Dau->p4();
	}

	int sign = part->pdgId() / abs(part->pdgId());
	reco::GenParticle TauH = reco::GenParticle (part->charge(), p4Had, part->vertex(), sign*66615, part->status(), true);
	return TauH;
}

reco::GenParticle IIHEModuleMCTruth::GetTauHadNeutrals (const reco::Candidate* part)
{
	if (abs(part->pdgId()) != 15)
	{
		reco::GenParticle fakeTauH = reco::GenParticle (0, reco::Candidate::LorentzVector(0.,0.,0.,0.), reco::Candidate::Point (0.,0.,0.), -999999, 0, true);
		std::cout << "Warning: building had tau from a particle with pdgId != 15 --> dummy entry returned" << std::endl;
		return fakeTauH;
	}
	reco::Candidate::LorentzVector p4Had (0,0,0,0);
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter( iDau );
		int dauId = abs(Dau->pdgId());
		if (dauId != 12 && dauId != 14 && dauId != 16 && Dau->charge()==0) // no neutrinos
			p4Had += Dau->p4();
	}

	int sign = part->pdgId() / abs(part->pdgId());
	reco::GenParticle TauH = reco::GenParticle (part->charge(), p4Had, part->vertex(), sign*77715, part->status(), true);
	return TauH;
}

void IIHEModuleMCTruth::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMCTruth::beginEvent(){}
void IIHEModuleMCTruth::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMCTruth::endJob(){
	addValueToMetaTree("mc_nEventsWeighted", nEventsWeighted_) ;
	addFVValueToMetaTree("mc_sumofgenWeights", sumofgenWeights_) ;
	parent_->saveToFile(pileupDist_) ;
}

DEFINE_FWK_MODULE(IIHEModuleMCTruth);
