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
	addBranch("GenPart_pdgId") ;
	addBranch("GenPart_status") ;
	addBranch("GenPart_statusFlags");
        addBranch("GenPart_genPartIdxMother") ;

	setBranchType(kVectorFloat) ;
	addBranch("GenPart_mass") ;
	addBranch("GenPart_pt") ;
	addBranch("GenPart_eta") ;
	addBranch("GenPart_phi") ;
	addBranch("GenPart_energy") ;

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

	for(GenParticleCollection::const_iterator mc_iter = genParticles.begin() ; mc_iter!=genParticles.end() ; ++mc_iter){
                store("GenPart_pt"     , mc_iter->pt()    ) ;
                store("GenPart_eta"    , mc_iter->eta()   ) ;
                store("GenPart_phi"    , mc_iter->phi()   ) ;
                store("GenPart_energy" , mc_iter->energy()) ;
                store("GenPart_mass"   , mc_iter->mass()  ) ;
                store("GenPart_pdgId"  , mc_iter->pdgId() ) ;
                store("GenPart_status" , mc_iter->status()) ;
                if (mc_iter->numberOfMothers()>0) store("GenPart_genPartIdxMother" ,mc_iter->motherRef(0).key()) ;
                else store("GenPart_genPartIdxMother" ,-1);
		// storing flags of a gen particle
		int flags = 0;

		if (mc_iter->statusFlags().isPrompt())                  flags |= (1 << 0);
		if (mc_iter->statusFlags().isDecayedLeptonHadron())     flags |= (1 << 1);
		if (mc_iter->statusFlags().isTauDecayProduct())         flags |= (1 << 2);
		if (mc_iter->statusFlags().isPromptTauDecayProduct())   flags |= (1 << 3);
		if (mc_iter->statusFlags().isDirectTauDecayProduct())   flags |= (1 << 4);
		if (mc_iter->statusFlags().isDirectPromptTauDecayProduct())       flags |= (1 << 5);
		if (mc_iter->statusFlags().isDirectHadronDecayProduct())          flags |= (1 << 6);
		if (mc_iter->statusFlags().isHardProcess())             flags |= (1 << 7);
		if (mc_iter->statusFlags().fromHardProcess())           flags |= (1 << 8);
		if (mc_iter->statusFlags().isHardProcessTauDecayProduct())       flags |= (1 << 9);
		if (mc_iter->statusFlags().isDirectHardProcessTauDecayProduct()) flags |= (1 << 10);
		if (mc_iter->statusFlags().fromHardProcessBeforeFSR())  flags |= (1 << 11);
                if (mc_iter->statusFlags().isFirstCopy())  flags |= (1 << 12);
                if (mc_iter->statusFlags().isLastCopy())  flags |= (1 << 13);
		if (mc_iter->statusFlags().isLastCopyBeforeFSR())       flags |= (1 << 14);

		store("GenPart_statusFlags" , flags);
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
