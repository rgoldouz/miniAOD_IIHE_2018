#include "UserCode/IIHETree/interface/IIHEModulePuppiJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;
//////////////////////////////////////////////////////////////////////////////////////////
////                                  Main IIHEPuppiJetModule
////////////////////////////////////////////////////////////////////////////////////////////

IIHEModulePuppiJet::IIHEModulePuppiJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC):IIHEModule(iConfig){
  pfPuppiJetLabel_                  =  iConfig.getParameter<edm::InputTag>("PuppiJetCollection"                  );
  pfPuppiJetLabelPrecor_            =  iConfig.getParameter<edm::InputTag>("PuppiJetCollectionPrecor"            );
  pfPuppiJetLabelSmeared_           =  iConfig.getParameter<edm::InputTag>("PuppiJetCollectionSmeared"           );
  pfPuppiJetLabelSmearedJetResUp_   =  iConfig.getParameter<edm::InputTag>("PuppiJetCollectionSmearedJetResUp"   );
  pfPuppiJetLabelSmearedJetResDown_ =  iConfig.getParameter<edm::InputTag>("PuppiJetCollectionSmearedJetResDown" );

  pfPuppiJetToken_                  =  iC.consumes<View<pat::Jet> > (pfPuppiJetLabel_                  );
  pfPuppiJetTokenPrecor_            =  iC.consumes<View<pat::Jet> > (pfPuppiJetLabelPrecor_            );
  pfPuppiJetTokenSmeared_           =  iC.consumes<View<pat::Jet> > (pfPuppiJetLabelSmeared_           );
  pfPuppiJetTokenSmearedJetResUp_   =  iC.consumes<View<pat::Jet> > (pfPuppiJetLabelSmearedJetResUp_   );
  pfPuppiJetTokenSmearedJetResDown_ =  iC.consumes<View<pat::Jet> > (pfPuppiJetLabelSmearedJetResDown_ );

  ETThreshold_ = iConfig.getUntrackedParameter<double>("puppiJetPtThreshold" );
  isMC_        = iConfig.getUntrackedParameter<bool>(  "isMC"                );

}
IIHEModulePuppiJet::~IIHEModulePuppiJet(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModulePuppiJet::beginJob(){
  addBranch("puppiJet_n", kUInt);

  setBranchType(kVectorFloat);
  addBranch("puppiJet_px"                          );
  addBranch("puppiJet_py"                          );
  addBranch("puppiJet_pz"                          );
  addBranch("puppiJet_pt"                          );
  addBranch("puppiJet_eta"                         );
  addBranch("puppiJet_theta"                       );
  addBranch("puppiJet_phi"                         );
  addBranch("puppiJet_energy"                      );
  addBranch("puppiJet_mass"                        );
  addBranch("puppiJet_chargedEmEnergyFraction"     );
  addBranch("puppiJet_neutralHadronEnergyFraction" );
  addBranch("puppiJet_neutralEmEnergyFraction"     );
  addBranch("puppiJet_chargedHadronEnergyFraction" );
  addBranch("puppiJet_muonEnergyFraction"          );
  addBranch("puppiJet_CSVv2"                       );
  addBranch("puppiJet_CvsL"                        );
  addBranch("puppiJet_CvsB"                        );
  addBranch("puppiJet_MVA2BJets"                   );
  addBranch("puppiJet_CvsB_DeepJet_charm_tagger"   );
  addBranch("puppiJet_CvsL_DeepJet_charm_tagger"   );
  addBranch("puppiJet_CvsB_DeepCSV_charm_tagger"   );
  addBranch("puppiJet_CvsL_DeepCSV_charm_tagger"   );
  addBranch("puppiJet_DeepJet"                     );
  addBranch("puppiJet_DeepCSV"                     );

  setBranchType(kVectorInt);
  addBranch("puppiJet_chargedMultiplicity"         );
  addBranch("puppiJet_neutralMultiplicity"         );
  addBranch("puppiJet_partonFlavour"               );
  addBranch("puppiJet_hadronFlavour"               );
  addBranch("puppiJet_isJetIDLoose_2016"           );
  addBranch("puppiJet_isJetIDTight_2016"           );
  addBranch("puppiJet_isJetIDTightLepVeto_2016"    );
  addBranch("puppiJet_isJetID_2017"                );
  addBranch("puppiJet_isJetIDLepVeto_2017"         );
  addBranch("puppiJet_isJetID_2018"                );
  addBranch("puppiJet_isJetIDLepVeto_2018"         );

  if (isMC_) {
  setBranchType(kVectorFloat);
  addBranch("puppiJet_Smeared_pt"                  );
  addBranch("puppiJet_SmearedJetResUp_pt"          );
  addBranch("puppiJet_SmearedJetResDown_pt"        );
  addBranch("puppiJet_SmearedJetEnUp_pt"           );
  addBranch("puppiJet_SmearedJetEnDown_pt"         );
  }

}

// ------------ method called to for each event  ------------
void IIHEModulePuppiJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Jet> > pfPuppiJetHandlePrecor_;
  edm::Handle<edm::View<pat::Jet> > pfPuppiJetHandle_;
  edm::Handle<edm::View<pat::Jet> > pfPuppiJetHandleSmeared_;
  edm::Handle<edm::View<pat::Jet> > pfPuppiJetHandleSmearedJetResUp_;
  edm::Handle<edm::View<pat::Jet> > pfPuppiJetHandleSmearedJetResDown_;

  iEvent.getByToken(pfPuppiJetTokenPrecor_            ,pfPuppiJetHandlePrecor_            );
  iEvent.getByToken(pfPuppiJetToken_                  ,pfPuppiJetHandle_                  );
  iEvent.getByToken(pfPuppiJetTokenSmeared_           ,pfPuppiJetHandleSmeared_           );
  iEvent.getByToken(pfPuppiJetTokenSmearedJetResUp_   ,pfPuppiJetHandleSmearedJetResUp_   );
  iEvent.getByToken(pfPuppiJetTokenSmearedJetResDown_ ,pfPuppiJetHandleSmearedJetResDown_ );

  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFPuppi", JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  for ( unsigned int i = 0; i <pfPuppiJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfPuppiJetHandle_->ptrAt( i );
    if(pfjet->pt() < ETThreshold_) continue ;
    store("puppiJet_px"         ,pfjet->px()     );
    store("puppiJet_py"         ,pfjet->py()     );
    store("puppiJet_pz"         ,pfjet->pz()     );
    store("puppiJet_pt"         ,pfjet->pt()     );
    store("puppiJet_pt_Precor"  ,pfPuppiJetHandlePrecor_->ptrAt(i)->pt() );
    store("puppiJet_eta"        ,pfjet->eta()    );
    store("puppiJet_theta"      ,pfjet->theta()  );
    store("puppiJet_phi"        ,pfjet->phi()    );
    store("puppiJet_energy"     ,pfjet->energy() );
    store("puppiJet_mass"       ,pfjet->mass()   );
    store("puppiJet_neutralHadronEnergyFraction" ,pfjet->neutralHadronEnergyFraction() );
    store("puppiJet_neutralEmEnergyFraction"     ,pfjet->neutralEmEnergyFraction()     );
    store("puppiJet_chargedHadronEnergyFraction" ,pfjet->chargedHadronEnergyFraction() );
    store("puppiJet_muonEnergyFraction"          ,pfjet->muonEnergyFraction()          );
    store("puppiJet_chargedEmEnergyFraction"     ,pfjet->chargedEmEnergyFraction()     );
    store("puppiJet_chargedMultiplicity"         ,pfjet->chargedMultiplicity()         );
    store("puppiJet_neutralMultiplicity"         ,pfjet->neutralMultiplicity()         );
    store("puppiJet_partonFlavour"               ,pfjet->partonFlavour()               );
    store("puppiJet_hadronFlavour"               ,pfjet->hadronFlavour()               );

    store("puppiJet_CSVv2"     ,pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    store("puppiJet_CvsL"      ,pfjet->bDiscriminator("pfCombinedCvsLJetTags")   );
    store("puppiJet_CvsB"      ,pfjet->bDiscriminator("pfCombinedCvsBJetTags")   );
    store("puppiJet_MVA2BJets" ,pfjet->bDiscriminator("pfCombinedMVAV2BJetTags") );
    store("puppiJet_DeepCSV"   ,pfjet->bDiscriminator("pfDeepCSVJetTags:probb")
                               +pfjet->bDiscriminator("pfDeepCSVJetTags:probbb")
    );
    store("puppiJet_DeepJet" ,pfjet->bDiscriminator("pfDeepFlavourJetTags:probb")
                             +pfjet->bDiscriminator("pfDeepFlavourJetTags:probbb")
                             +pfjet->bDiscriminator("pfDeepFlavourJetTags:problepb")
    );
    store("puppiJet_CvsL_DeepCSV_charm_tagger",
      pfjet->bDiscriminator("pfDeepCSVJetTags:probc") / (
      pfjet->bDiscriminator("pfDeepCSVJetTags:probc")
     +pfjet->bDiscriminator("pfDeepCSVJetTags:probudsg") )
    );

    store("puppiJet_CvsB_DeepCSV_charm_tagger",
      pfjet->bDiscriminator("pfDeepCSVJetTags:probc") / (
      pfjet->bDiscriminator("pfDeepCSVJetTags:probc")
     +pfjet->bDiscriminator("pfDeepCSVJetTags:probb")
     +pfjet->bDiscriminator("pfDeepCSVJetTags:probbb")  )
    );

    store("puppiJet_CvsL_DeepJet_charm_tagger",
      pfjet->bDiscriminator("pfDeepFlavourJetTags:probc") / (
      pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")
     +pfjet->bDiscriminator("pfDeepFlavourJetTags:probuds")
     +pfjet->bDiscriminator("pfDeepFlavourJetTags:probg")   )
    );

    store("puppiJet_CvsB_DeepJet_charm_tagger",
      pfjet->bDiscriminator("pfDeepFlavourJetTags:probc") /  (
      pfjet->bDiscriminator("pfDeepFlavourJetTags:probc")
     +pfjet->bDiscriminator("pfDeepFlavourJetTags:probb")
     +pfjet->bDiscriminator("pfDeepFlavourJetTags:probbb")
     +pfjet->bDiscriminator("pfDeepFlavourJetTags:problepb") )
    );

    float eta               = pfjet->eta();
    float NHF               = pfjet->neutralHadronEnergyFraction();
    float NEMF              = pfjet->neutralEmEnergyFraction();
    float CHF               = pfjet->chargedHadronEnergyFraction();
    float MUF               = pfjet->muonEnergyFraction();
    float CEMF              = pfjet->chargedEmEnergyFraction();
    int NumConst            = pfjet->chargedMultiplicity()+pfjet->neutralMultiplicity();
    int NumNeutralParticles = pfjet->neutralMultiplicity();
    int CHM                 = pfjet->chargedMultiplicity();

    bool isJetIDLoose_2016        ;
    bool isJetIDTight_2016        ;
    bool isJetIDTightLepVeto_2016 ;
    bool isJetID_2017             ;
    bool isJetIDLepVeto_2017      ;
    bool isJetID_2018             ;
    bool isJetIDLepVeto_2018      ;

//****** 2016 IDs **** https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016
    if( fabs(eta) <= 2.7 ) {
      isJetIDLoose_2016 = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) );
      isJetIDTight_2016 = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) );
      isJetIDTightLepVeto_2016 = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) );
    }
    else if ( abs(eta)>2.7 && abs(eta)<=3.0 ) {
      isJetIDLoose_2016 = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTight_2016 = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTightLepVeto_2016 = false;
    }
    else {
      isJetIDLoose_2016 = (NEMF<0.90 && NumNeutralParticles>10 );
      isJetIDTight_2016 = (NEMF<0.90 && NumNeutralParticles>10 ) ;
      isJetIDTightLepVeto_2016 = false;
    }
    store("puppiJet_isJetIDLoose_2016"        ,int(isJetIDLoose_2016)        );
    store("puppiJet_isJetIDTight_2016"        ,int(isJetIDTight_2016)        );
    store("puppiJet_isJetIDTightLepVeto_2016" ,int(isJetIDTightLepVeto_2016) );

//****** 2017 IDs **** https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    if( fabs(eta)<=2.4 ) {
      isJetID_2017        = (abs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2017 = (abs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
    }
    else if ( abs(eta)>2.4 && abs(eta)<=2.7 ) {
      isJetID_2017         = ( abs(eta)>2.4 && abs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2017  = ( abs(eta)>2.4 && abs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
    }
    else if ( abs(eta)>2.7 && abs(eta)<=3.0 ) {
      isJetID_2017         = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      isJetIDLepVeto_2017  = false;
    }
    else {
      isJetID_2017         = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 && abs(eta)>3.0 );
      isJetIDLepVeto_2017  = false;
    }
    store("puppiJet_isJetIDLepVeto_2017" ,int(isJetIDLepVeto_2017) );
    store("puppiJet_isJetID_2017"        ,int(isJetID_2017)        );


//****** 2018 IDs **** https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    if( fabs(eta)<=2.6){
      isJetID_2018        = (abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      isJetIDLepVeto_2018 = (abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.6 && abs(eta)<=2.7 ){
      isJetID_2018         = ( abs(eta)>2.6 && abs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 );
      isJetIDLepVeto_2018  = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetID_2018         = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      isJetIDLepVeto_2018  = false;
      }
    else{
      isJetID_2018         = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && abs(eta)>3.0 );
      isJetIDLepVeto_2018  = false;
    }
    store("jet_isJetIDLepVeto_2018" ,int(isJetIDLepVeto_2018) );
    store("jet_isJetID_2018"        ,int(isJetID_2018)        );

    if (isMC_) {
      jecUnc.setJetEta(pfPuppiJetHandleSmeared_->at(i).eta() );
      jecUnc.setJetPt( pfPuppiJetHandleSmeared_->at(i).pt()  );
      double unc = jecUnc.getUncertainty(true);
      store("puppiJet_SmearedJetEnUp_pt"   ,(1+unc)*pfPuppiJetHandleSmeared_  ->at(i).pt() );
      store("puppiJet_SmearedJetEnDown_pt" ,(1-unc)*pfPuppiJetHandleSmeared_  ->at(i).pt() );
      store("puppiJet_Smeared_pt"          ,pfPuppiJetHandleSmeared_          ->at(i).pt() );
      store("puppiJet_SmearedJetResUp_pt"  ,pfPuppiJetHandleSmearedJetResUp_  ->at(i).pt() );
      store("puppiJet_SmearedJetResDown_pt",pfPuppiJetHandleSmearedJetResDown_->at(i).pt() );
    }

  }
}
void IIHEModulePuppiJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModulePuppiJet::beginEvent(){}
void IIHEModulePuppiJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModulePuppiJet::endJob(){
}

DEFINE_FWK_MODULE(IIHEModulePuppiJet);
