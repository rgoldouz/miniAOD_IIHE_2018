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
  addBranch("jet_CSVv2");
  addBranch("jet_CvsL");
  addBranch("jet_CvsB");
  addBranch("jet_MVA2BJets");
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
  addBranch("jet_Smeared_energy");
  addBranch("jet_SmearedJetResUp_pt");
  addBranch("jet_SmearedJetResUp_energy");
  addBranch("jet_SmearedJetResDown_pt");
  addBranch("jet_SmearedJetResDown_energy");
  addBranch("jet_EnUp_pt");
  addBranch("jet_EnUp_energy");
  addBranch("jet_EnDown_pt");
  addBranch("jet_EnDown_energy");

  addBranch("jet_BtagSF_loose");
  addBranch("jet_BtagSFbcUp_loose");
  addBranch("jet_BtagSFbcDown_loose");
  addBranch("jet_BtagSFudsgUp_loose");
  addBranch("jet_BtagSFudsgDown_loose");

  addBranch("jet_BtagSF_medium");
  addBranch("jet_BtagSFbcUp_medium");
  addBranch("jet_BtagSFbcDown_medium");
  addBranch("jet_BtagSFudsgUp_medium");
  addBranch("jet_BtagSFudsgDown_medium");

  addBranch("jet_BtagSF_tight");
  addBranch("jet_BtagSFbcUp_tight");
  addBranch("jet_BtagSFbcDown_tight");
  addBranch("jet_BtagSFudsgUp_tight");
  addBranch("jet_BtagSFudsgDown_tight");
  }

}

// ------------ method called to for each event  ------------
void IIHEModuleJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  edm::Handle<edm::View<pat::Jet> > pfJetHandle_;
  iEvent.getByToken(pfJetToken_, pfJetHandle_);

  const string ctr = "central";
  const string vup = "up";
  const string vdown = "down";


  for ( unsigned int i = 0; i <pfJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfJetHandle_->ptrAt( i );

    if(pfjet->pt() < ETThreshold_) continue ;

    store("jet_px"    , pfjet->px()) ;
    store("jet_py"    , pfjet->py()) ;
    store("jet_pz"    , pfjet->pz()) ;
    store("jet_pt"    , pfjet->pt()) ;
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

//JetBTagWeight( edm::View<pat::Jet>&b, size_t ijet, const vector<BTagEntry::OperatingPoinconst string &bc_full_syst, const string &udsg_full_syst,const string &bc_full_syst, const string &udsg_full_syst,const string &bc_fast_syst, const string &udsg_fast_syst,bool do_deep_csv, bool do_by_proc, Runs runs)


   }

  }
}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleJet);
