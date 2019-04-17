#include "UserCode/IIHETree/interface/IIHEModuleJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;
//////////////////////////////////////////////////////////////////////////////////////////
////                                  Main IIHEJetModule
////////////////////////////////////////////////////////////////////////////////////////////

IIHEModuleJet::IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC):IIHEModule(iConfig){
  pfJetLabel_                  =  iConfig.getParameter<edm::InputTag>("JetCollection");
  pfJetToken_                  =  iC.consumes<View<pat::Jet> > (pfJetLabel_);
  pfJetLabelSmeared_           =  iConfig.getParameter<edm::InputTag>("JetCollectionSmeared");
  pfJetTokenSmeared_           =  iC.consumes<View<pat::Jet> > (pfJetLabelSmeared_);
  pfJetLabelSmearedJetResUp_   =  iConfig.getParameter<edm::InputTag>("JetCollectionSmearedJetResUp");
  pfJetTokenSmearedJetResUp_   =  iC.consumes<View<pat::Jet> > (pfJetLabelSmearedJetResUp_);
  pfJetLabelSmearedJetResDown_ =  iConfig.getParameter<edm::InputTag>("JetCollectionSmearedJetResDown");
  pfJetTokenSmearedJetResDown_ = iC.consumes<View<pat::Jet> > (pfJetLabelSmearedJetResDown_);
  pfJetTokenPrecor_                  =  iC.consumes<View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetCollectionPrecor"));

  looseBtagSFdownToken_      =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("looseBtagSFdown"));
  looseBtagSFnominalToken_   =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("looseBtagSFnominal"));
  looseBtagSFupToken_        =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("looseBtagSFup"));
  mediumBtagSFdownToken_     =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("mediumBtagSFdown"));
  mediumBtagSFnominalToken_  =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("mediumBtagSFnominal"));
  mediumBtagSFupToken_       =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("mediumBtagSFup"));
  tightBtagSFdownToken_      =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("tightBtagSFdown"));
  tightBtagSFnominalToken_   =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("tightBtagSFnominal"));
  tightBtagSFupToken_        =  iC.consumes<nanoaod::FlatTable>  (iConfig.getParameter<edm::InputTag>("tightBtagSFup"));

  ETThreshold_ = iConfig.getUntrackedParameter<double>("jetPtThreshold") ;
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC") ;

}
IIHEModuleJet::~IIHEModuleJet(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleJet::beginJob(){
  addBranch("jet_n", kUInt);

  setBranchType(kVectorFloat);
  addBranch("jet_px");
  addBranch("jet_py");
  addBranch("jet_pz");
  addBranch("jet_pt");
  addBranch("jet_eta");
  addBranch("jet_theta");
  addBranch("jet_phi");
  addBranch("jet_energy");
  addBranch("jet_mass");
  addBranch("jet_chargedEmEnergyFraction");
  addBranch("jet_neutralHadronEnergyFraction");
  addBranch("jet_neutralEmEnergyFraction");
  addBranch("jet_chargedHadronEnergyFraction");
  addBranch("jet_muonEnergyFraction");
  setBranchType(kVectorInt);
  addBranch("jet_chargedMultiplicity");
  addBranch("jet_neutralMultiplicity");
  addBranch("jet_partonFlavour");
  addBranch("jet_hadronFlavour");

  setBranchType(kVectorFloat);
  setBranchType(kVectorFloat);
  addBranch("jet_CSVv2");
  addBranch("jet_CvsL");
  addBranch("jet_CvsB");
  addBranch("jet_MVA2BJets");
  addBranch("jet_CvsB_DeepJet_charm_tagger");
  addBranch("jet_CvsL_DeepJet_charm_tagger");
  addBranch("jet_CvsB_DeepCSV_charm_tagger");
  addBranch("jet_CvsL_DeepCSV_charm_tagger");
  addBranch("jet_DeepJet");
  addBranch("jet_DeepCSV");
  setBranchType(kVectorInt);
  addBranch("jet_isJetIDLoose_2016");
  addBranch("jet_isJetIDTight_2016");
  addBranch("jet_isJetIDTightLepVeto_2016");
  addBranch("jet_isJetID_2017");
  addBranch("jet_isJetIDLepVeto_2017");
  addBranch("jet_isJetID_2018");
  addBranch("jet_isJetIDLepVeto_2018");
  if (isMC_){

  setBranchType(kVectorFloat);
  addBranch("jet_Smeared_pt");
  addBranch("jet_SmearedJetResUp_pt");
  addBranch("jet_SmearedJetResDown_pt");
  addBranch("jet_EnUp_pt");
  addBranch("jet_EnDown_pt");

  setBranchType(kFloat);
  addBranch("BtagSF_Deepcsv_loose");
  addBranch("BtagSF_Deepcsv_Up_loose");
  addBranch("BtagSF_Deepcsv_Down_loose");

  addBranch("BtagSF_Deepcsv_medium");
  addBranch("BtagSF_Deepcsv_Up_medium");
  addBranch("BtagSF_Deepcsv_Down_medium");

  addBranch("BtagSF_Deepcsv_tight");
  addBranch("BtagSF_Deepcsv_Up_tight");
  addBranch("BtagSF_Deepcsv_Down_tight");
  }

}

// ------------ method called to for each event  ------------
void IIHEModuleJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Jet> > pfJetHandlePrecor_;
  iEvent.getByToken(pfJetTokenPrecor_, pfJetHandlePrecor_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandle_;
  iEvent.getByToken(pfJetToken_, pfJetHandle_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmeared_;
  iEvent.getByToken(pfJetTokenSmeared_, pfJetHandleSmeared_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmearedJetResUp_;
  iEvent.getByToken(pfJetTokenSmearedJetResUp_, pfJetHandleSmearedJetResUp_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmearedJetResDown_;
  iEvent.getByToken(pfJetTokenSmearedJetResDown_, pfJetHandleSmearedJetResDown_);


  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  for ( unsigned int i = 0; i <pfJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfJetHandle_->ptrAt( i );
    if(pfjet->pt() < ETThreshold_) continue ;
    store("jet_px"    , pfjet->px()) ;
    store("jet_py"    , pfjet->py()) ;
    store("jet_pz"    , pfjet->pz()) ;
    store("jet_pt"    , pfjet->pt()) ;
    store("jet_pt_Precor"    , pfJetHandlePrecor_->ptrAt( i )->pt()) ;
    store("jet_eta"   , pfjet->eta()) ;
    store("jet_theta" , pfjet->theta()) ;
    store("jet_phi"   , pfjet->phi()) ;
    store("jet_energy", pfjet->energy()) ;
    store("jet_mass"  , pfjet->mass()) ;
    store("jet_neutralHadronEnergyFraction"         ,pfjet->neutralHadronEnergyFraction());
    store("jet_neutralEmEnergyFraction"             ,pfjet->neutralEmEnergyFraction());
    store("jet_chargedHadronEnergyFraction"         ,pfjet->chargedHadronEnergyFraction());
    store("jet_muonEnergyFraction"                  ,pfjet->muonEnergyFraction());
    store("jet_chargedEmEnergyFraction"             ,pfjet->chargedEmEnergyFraction());
    store("jet_chargedMultiplicity"                 ,pfjet->chargedMultiplicity());
    store("jet_neutralMultiplicity"                 ,pfjet->neutralMultiplicity());
    store("jet_partonFlavour"                       ,pfjet->partonFlavour());
    store("jet_hadronFlavour"                       ,pfjet->hadronFlavour());

    store("jet_CSVv2",pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    store("jet_CvsL",pfjet->bDiscriminator("pfCombinedCvsLJetTags"));
    store("jet_CvsB",pfjet->bDiscriminator("pfCombinedCvsBJetTags"));
    store("jet_MVA2BJets", pfjet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    store("jet_DeepCSV", pfjet->bDiscriminator("pfDeepCSVJetTags:probb") + pfjet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    store("jet_DeepJet", pfjet->bDiscriminator("pfDeepFlavourJetTags:probb") + pfjet->bDiscriminator("pfDeepFlavourJetTags:probbb") + pfjet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    store("jet_CvsL_DeepCSV_charm_tagger", pfjet->bDiscriminator("pfDeepCSVJetTags:probc")/(pfjet->bDiscriminator("pfDeepCSVJetTags:probc") + pfjet->bDiscriminator("pfDeepCSVJetTags:probudsg")));

    store("jet_CvsB_DeepCSV_charm_tagger", pfjet->bDiscriminator("pfDeepCSVJetTags:probc")/(pfjet->bDiscriminator("pfDeepCSVJetTags:probc")+pfjet->bDiscriminator("pfDeepCSVJetTags:probb")+pfjet->bDiscriminator("pfDeepCSVJetTags:probbb")));

    store("jet_CvsL_DeepJet_charm_tagger", pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")/(pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")+pfjet->bDiscriminator("pfDeepFlavourJetTags:probuds")+pfjet->bDiscriminator("pfDeepFlavourJetTags:probg")));

    store("jet_CvsB_DeepJet_charm_tagger", pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")/(pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")+pfjet->bDiscriminator("pfDeepFlavourJetTags:probb")+pfjet->bDiscriminator("pfDeepFlavourJetTags:probbb")+ pfjet->bDiscriminator("pfDeepFlavourJetTags:problepb")));

    float eta = pfjet->eta();
    float NHF = pfjet->neutralHadronEnergyFraction();
    float NEMF = pfjet->neutralEmEnergyFraction();
    float CHF = pfjet->chargedHadronEnergyFraction();
    float MUF = pfjet->muonEnergyFraction();
    float CEMF = pfjet->chargedEmEnergyFraction();
    int NumConst = pfjet->chargedMultiplicity()+pfjet->neutralMultiplicity();
    int NumNeutralParticles =pfjet->neutralMultiplicity();
    int CHM = pfjet->chargedMultiplicity(); 

    bool isJetIDLoose_2016;
    bool isJetIDTight_2016;
    bool isJetIDTightLepVeto_2016;
    bool isJetID_2017;
    bool isJetIDLepVeto_2017;
    bool isJetID_2018;
    bool isJetIDLepVeto_2018;

//****** 2016 IDs **** https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016
    if( fabs(eta) <= 2.7){
      isJetIDLoose_2016 = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTight_2016 = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTightLepVeto_2016 = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4));
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetIDLoose_2016 = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTight_2016 = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTightLepVeto_2016 = false;
      }
    else{
      isJetIDLoose_2016 = (NEMF<0.90 && NumNeutralParticles>10 );
      isJetIDTight_2016 = (NEMF<0.90 && NumNeutralParticles>10 ) ;
      isJetIDTightLepVeto_2016 = false;
    }
    store("jet_isJetIDLoose_2016"                 ,int(isJetIDLoose_2016));
    store("jet_isJetIDTight_2016"                 ,int(isJetIDTight_2016));
    store("jet_isJetIDTightLepVeto_2016"          ,int(isJetIDTightLepVeto_2016));

//****** 2017 IDs **** https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    if( fabs(eta)<=2.4){
      isJetID_2017 = (abs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2017 = (abs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.4 && abs(eta)<=2.7 ){
      isJetID_2017  = ( abs(eta)>2.4 && abs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2017  = ( abs(eta)>2.4 && abs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetID_2017 = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      isJetIDLepVeto_2017  = false;
      }
    else{
      isJetID_2017 = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 && abs(eta)>3.0 );
      isJetIDLepVeto_2017  = false;
    }
    store("jet_isJetIDLepVeto_2017"                 ,int(isJetIDLepVeto_2017));
    store("jet_isJetID_2017"                        ,int(isJetID_2017));


//****** 2018 IDs **** https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    if( fabs(eta)<=2.6){
      isJetID_2018 = (abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2018 = (abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.6 && abs(eta)<=2.7 ){
      isJetID_2018  = ( abs(eta)>2.6 && abs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 );
      isJetIDLepVeto_2018  = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetID_2018 = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      isJetIDLepVeto_2018  = false;
      }
    else{
      isJetID_2018 = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && abs(eta)>3.0 );
      isJetIDLepVeto_2018  = false;
    }
    store("jet_isJetIDLepVeto_2018"                 ,int(isJetIDLepVeto_2018));
    store("jet_isJetID_2018"                        ,int(isJetID_2018));



    if (isMC_){
      jecUnc.setJetEta(pfjet->eta());
      jecUnc.setJetPt(pfjet->pt());
      double unc = jecUnc.getUncertainty(true);
      store("jet_EnUp_pt", (1+unc)*pfjet->pt());
      store("jet_EnDown_pt", (1-unc)*pfjet->pt());
      store("jet_Smeared_pt",pfJetHandleSmeared_->at(i).pt());
      store("jet_SmearedJetResUp_pt",pfJetHandleSmearedJetResUp_->at(i).pt());
      store("jet_SmearedJetResDown_pt",pfJetHandleSmearedJetResDown_->at(i).pt());
    }

  }

  if (isMC_){
    edm::Handle<nanoaod::FlatTable>  looseBtagSFdownHandle_;
    iEvent.getByToken(looseBtagSFdownToken_, looseBtagSFdownHandle_);

    edm::Handle<nanoaod::FlatTable> looseBtagSFnominalHandle_;
    iEvent.getByToken(looseBtagSFnominalToken_, looseBtagSFnominalHandle_);

    edm::Handle<nanoaod::FlatTable>  looseBtagSFupHandle_;
    iEvent.getByToken(looseBtagSFupToken_, looseBtagSFupHandle_);

    edm::Handle<nanoaod::FlatTable> mediumBtagSFdownHandle_;
    iEvent.getByToken(mediumBtagSFdownToken_, mediumBtagSFdownHandle_);

    edm::Handle<nanoaod::FlatTable> mediumBtagSFnominalHandle_;
    iEvent.getByToken(mediumBtagSFnominalToken_, mediumBtagSFnominalHandle_);

    edm::Handle<nanoaod::FlatTable> mediumBtagSFupHandle_;
    iEvent.getByToken(mediumBtagSFupToken_, mediumBtagSFupHandle_);

    edm::Handle<nanoaod::FlatTable> tightBtagSFdownHandle_;
    iEvent.getByToken(tightBtagSFdownToken_, tightBtagSFdownHandle_);

    edm::Handle<nanoaod::FlatTable> tightBtagSFnominalHandle_;
    iEvent.getByToken(tightBtagSFnominalToken_, tightBtagSFnominalHandle_);

    edm::Handle<nanoaod::FlatTable> tightBtagSFupHandle_;
    iEvent.getByToken(tightBtagSFupToken_, tightBtagSFupHandle_);

    store("BtagSF_Deepcsv_loose"       ,looseBtagSFnominalHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Up_loose"    ,looseBtagSFupHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Down_loose"  ,looseBtagSFdownHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_medium"      ,mediumBtagSFnominalHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Up_medium"   ,mediumBtagSFupHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Down_medium" ,mediumBtagSFdownHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_tight"       ,tightBtagSFnominalHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Up_tight"    ,tightBtagSFupHandle_->getAnyValue(0,0));
    store("BtagSF_Deepcsv_Down_tight"  ,tightBtagSFdownHandle_->getAnyValue(0,0));
  }
}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleJet);
