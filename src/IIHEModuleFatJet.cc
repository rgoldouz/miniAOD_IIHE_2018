#include "UserCode/IIHETree/interface/IIHEModuleFatJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;
//////////////////////////////////////////////////////////////////////////////////////////
////                                  Main IIHEFatJetModule
//////////////////////////////////////////////////////////////////////////////////////////

IIHEModuleFatJet::IIHEModuleFatJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC):IIHEModule(iConfig){

  FatJetLabel_                  = iConfig.getParameter<edm::InputTag>("DeepAK8JetCollection"                 );
  FatJetSmearedLabel_           = iConfig.getParameter<edm::InputTag>("DeepAK8JetCollectionSmeared"          );
  FatJetSmearedJetResUpLabel_   = iConfig.getParameter<edm::InputTag>("DeepAK8JetCollectionSmearedJetResUp"  );
  FatJetSmearedJetResDownLabel_ = iConfig.getParameter<edm::InputTag>("DeepAK8JetCollectionSmearedJetResDown");
  FatJetToken_                  = iC.consumes<View<pat::Jet> >       ( FatJetLabel_                      );
  FatJetSmearedToken_           = iC.consumes<View<pat::Jet> >       ( FatJetSmearedLabel_               );
  FatJetSmearedJetResUpToken_   = iC.consumes<View<pat::Jet> >       ( FatJetSmearedJetResUpLabel_       );
  FatJetSmearedJetResDownToken_ = iC.consumes<View<pat::Jet> >       ( FatJetSmearedJetResDownLabel_     );

  ETThreshold_ = iConfig.getUntrackedParameter<double>("fatjetPtThreshold");
  isMC_        = iConfig.getUntrackedParameter<bool>  ("isMC"             );

}
IIHEModuleFatJet::~IIHEModuleFatJet(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleFatJet::beginJob(){
  setBranchType(kVectorFloat);
  // kinematics
  addBranch("fatjet_px"    );
  addBranch("fatjet_py"    );
  addBranch("fatjet_pz"    );
  addBranch("fatjet_pt"    );
  addBranch("fatjet_eta"   );
  addBranch("fatjet_theta" );
  addBranch("fatjet_phi"   );
  addBranch("fatjet_energy");
  addBranch("fatjet_mass"  );
  addBranch("fatjet_area"  );
  addBranch("fatjet_pt_raw"    );
  addBranch("fatjet_mass_raw"  );
  addBranch("fatjet_energy_raw");
  if (isMC_) {
    setBranchType(kVectorFloat);
    addBranch("fatjet_Smeared_pt"           );
    addBranch("fatjet_SmearedJetResUp_pt"   );
    addBranch("fatjet_SmearedJetResDown_pt" );
    addBranch("fatjet_SmearedJetEnUp_pt"    );
    addBranch("fatjet_SmearedJetEnDown_pt"  );
  }

  // N-jettiness
  addBranch("fatjet_puppi_tau1" );
  addBranch("fatjet_puppi_tau2" );
  addBranch("fatjet_puppi_tau3" );
  addBranch("fatjet_puppi_tau4" );
  addBranch("fatjet_chs_tau1"   );
  addBranch("fatjet_chs_tau2"   );
  addBranch("fatjet_chs_tau3"   );
  addBranch("fatjet_chs_tau4"   );

  // Energy correlation function variables
  addBranch("fatjet_puppi_n2b1" );
  addBranch("fatjet_puppi_n3b1" );
  addBranch("fatjet_puppi_n2b2" );
  addBranch("fatjet_puppi_n3b2" );

  // Softdrop and pruned masses
  addBranch("fatjet_puppi_softdrop_mass");
  addBranch("fatjet_chs_softdrop_mass"  );
  addBranch("fatjet_chs_pruned_mass"    );

  // b Tagger
  addBranch("fatjet_DeepCSV_b"            );
  addBranch("fatjet_DeepCSV_bb"           );
  addBranch("fatjet_DeepCSV_c"            );
  addBranch("fatjet_DeepCSV_udsg"         );
  addBranch("fatjet_DeepCSV_bTag"         );

  addBranch("fatjet_DeepFlavour_b"        );
  addBranch("fatjet_DeepFlavour_bb"       );
  addBranch("fatjet_DeepFlavour_lepb"     );
  addBranch("fatjet_DeepFlavour_c"        );
  addBranch("fatjet_DeepFlavour_uds"      );
  addBranch("fatjet_DeepFlavour_g"        );
  addBranch("fatjet_DeepFlavour_bTag"     );
  addBranch("fatjet_DeepFlavour_bTagLepb" );

  addBranch("fatjet_bTag_CMVA"            );
  addBranch("fatjet_bTag_CSVV2"           );
  addBranch("fatjet_bTag_Hbb"             );

  addBranch("fatjet_DeepDouble_BvL_QCD"   );
  addBranch("fatjet_DeepDouble_BvL_Hbb"   );
  addBranch("fatjet_DeepDouble_CvL_QCD"   );
  addBranch("fatjet_DeepDouble_CvL_Hcc"   );
  addBranch("fatjet_DeepDouble_CvB_Hbb"   );
  addBranch("fatjet_DeepDouble_CvB_Hcc"   );
  addBranch("fatjet_DeepDoubleMD_BvL_QCD" );
  addBranch("fatjet_DeepDoubleMD_BvL_Hbb" );
  addBranch("fatjet_DeepDoubleMD_CvL_QCD" );
  addBranch("fatjet_DeepDoubleMD_CvL_Hcc" );
  addBranch("fatjet_DeepDoubleMD_CvB_Hbb" );
  addBranch("fatjet_DeepDoubleMD_CvB_Hcc" );

  // DeepAK8 Nominal Tagger
  // // Aggregated
  addBranch("fatjet_DeepBoosted_TvsQCD"    );
  addBranch("fatjet_DeepBoosted_WvsQCD"    );
  addBranch("fatjet_DeepBoosted_ZvsQCD"    );
  addBranch("fatjet_DeepBoosted_ZbbvsQCD"  );
  addBranch("fatjet_DeepBoosted_HbbvsQCD"  );
  addBranch("fatjet_DeepBoosted_H4qcsQCD"  );
  // // Raw
  addBranch("fatjet_DeepBoosted_Tbcq"      );
  addBranch("fatjet_DeepBoosted_Tbqq"      );
  addBranch("fatjet_DeepBoosted_Tbc"       );
  addBranch("fatjet_DeepBoosted_Tbq"       );
  addBranch("fatjet_DeepBoosted_Wcq"       );
  addBranch("fatjet_DeepBoosted_Wqq"       );
  addBranch("fatjet_DeepBoosted_Zbb"       );
  addBranch("fatjet_DeepBoosted_Zcc"       );
  addBranch("fatjet_DeepBoosted_Zqq"       );
  addBranch("fatjet_DeepBoosted_Hbb"       );
  addBranch("fatjet_DeepBoosted_Hcc"       );
  addBranch("fatjet_DeepBoosted_H4q"       );
  addBranch("fatjet_DeepBoosted_QCDbb"     );
  addBranch("fatjet_DeepBoosted_QCDcc"     );
  addBranch("fatjet_DeepBoosted_QCDb"      );
  addBranch("fatjet_DeepBoosted_QCDc"      );
  addBranch("fatjet_DeepBoosted_QCDothers" );
  // // Raw Combined
  addBranch("fatjet_DeepBoosted_T"         );
  addBranch("fatjet_DeepBoosted_W"         );
  addBranch("fatjet_DeepBoosted_Z"         );
  addBranch("fatjet_DeepBoosted_H"         );
  addBranch("fatjet_DeepBoosted_QCD"       );

  // DeepAK8 Mass-Decorrelated Tagger
  // // Aggregated
  addBranch("fatjet_DeepBoostedMD_TvsQCD"    );
  addBranch("fatjet_DeepBoostedMD_WvsQCD"    );
  addBranch("fatjet_DeepBoostedMD_ZvsQCD"    );
  addBranch("fatjet_DeepBoostedMD_ZbbvsQCD"  );
  addBranch("fatjet_DeepBoostedMD_HbbvsQCD"  );
  addBranch("fatjet_DeepBoostedMD_H4qvsQCD"  );
  addBranch("fatjet_DeepBoostedMD_ZHbbvsQCD" );
  addBranch("fatjet_DeepBoostedMD_ZHccvsQCD" );
  addBranch("fatjet_DeepBoostedMD_bbvsLight" );
  addBranch("fatjet_DeepBoostedMD_ccvsLight" );
  // // Raw
  addBranch("fatjet_DeepBoostedMD_Tbcq"      );
  addBranch("fatjet_DeepBoostedMD_Tbqq"      );
  addBranch("fatjet_DeepBoostedMD_Tbc"       );
  addBranch("fatjet_DeepBoostedMD_Tbq"       );
  addBranch("fatjet_DeepBoostedMD_Wcq"       );
  addBranch("fatjet_DeepBoostedMD_Wqq"       );
  addBranch("fatjet_DeepBoostedMD_Zbb"       );
  addBranch("fatjet_DeepBoostedMD_Zcc"       );
  addBranch("fatjet_DeepBoostedMD_Zqq"       );
  addBranch("fatjet_DeepBoostedMD_Hbb"       );
  addBranch("fatjet_DeepBoostedMD_Hcc"       );
  addBranch("fatjet_DeepBoostedMD_H4q"       );
  addBranch("fatjet_DeepBoostedMD_QCDbb"     );
  addBranch("fatjet_DeepBoostedMD_QCDcc"     );
  addBranch("fatjet_DeepBoostedMD_QCDb"      );
  addBranch("fatjet_DeepBoostedMD_QCDc"      );
  addBranch("fatjet_DeepBoostedMD_QCDothers" );
  // // Raw Combined
  addBranch("fatjet_DeepBoostedMD_T"         );
  addBranch("fatjet_DeepBoostedMD_W"         );
  addBranch("fatjet_DeepBoostedMD_Z"         );
  addBranch("fatjet_DeepBoostedMD_H"         );
  addBranch("fatjet_DeepBoostedMD_QCD"       );
}

// ------------ method called to for each event  ------------
void IIHEModuleFatJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // ------------------------------
  // Handles
    edm::Handle<edm::View<pat::Jet> >                FatJetHandle_                   ;
    edm::Handle<edm::View<pat::Jet> >                FatJetSmearedHandle_            ;
    edm::Handle<edm::View<pat::Jet> >                FatJetSmearedJetResUpHandle_    ;
    edm::Handle<edm::View<pat::Jet> >                FatJetSmearedJetResDownHandle_  ;
    iEvent.getByToken(FatJetToken_                  ,FatJetHandle_                  );
    iEvent.getByToken(FatJetSmearedToken_           ,FatJetSmearedHandle_           );
    iEvent.getByToken(FatJetSmearedJetResUpToken_   ,FatJetSmearedJetResUpHandle_   );
    iEvent.getByToken(FatJetSmearedJetResDownToken_ ,FatJetSmearedJetResDownHandle_ );

  // ------------------------------
  // Corrections
    ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get("AK8PFPuppi", JetCorParColl);
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    JetCorrectionUncertainty jecUnc(JetCorPar);

  // ------------------------------------------------------------
  // Loop over Fat Jets
  for ( unsigned int i = 0; i <FatJetHandle_->size(); ++i) {
    Ptr<pat::Jet> fatjet = FatJetHandle_->ptrAt(i);
    if(fatjet->pt() < ETThreshold_) continue ;
    store("fatjet_px"    ,fatjet->px()      );
    store("fatjet_py"    ,fatjet->py()      );
    store("fatjet_pz"    ,fatjet->pz()      );
    store("fatjet_pt"    ,fatjet->pt()      );
    store("fatjet_eta"   ,fatjet->eta()     );
    store("fatjet_theta" ,fatjet->theta()   );
    store("fatjet_phi"   ,fatjet->phi()     );
    store("fatjet_energy",fatjet->energy()  );
    store("fatjet_mass"  ,fatjet->mass()    );
    store("fatjet_area"  ,fatjet->jetArea() );

    // Uncorrected FatJets
    store("fatjet_pt_raw"    ,fatjet->pt()     * fatjet->jecFactor("Uncorrected") );
    store("fatjet_mass_raw"  ,fatjet->mass()   * fatjet->jecFactor("Uncorrected") );
    store("fatjet_energy_raw",fatjet->energy() * fatjet->jecFactor("Uncorrected") );

    // N-jettiness
    store("fatjet_puppi_tau1",fatjet->userFloat("NjettinessAK8Puppi:tau1") );
    store("fatjet_puppi_tau2",fatjet->userFloat("NjettinessAK8Puppi:tau2") );
    store("fatjet_puppi_tau3",fatjet->userFloat("NjettinessAK8Puppi:tau3") );
    store("fatjet_puppi_tau4",fatjet->userFloat("NjettinessAK8Puppi:tau4") );
    store("fatjet_chs_tau1"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") );
    store("fatjet_chs_tau2"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") );
    store("fatjet_chs_tau3"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") );
    store("fatjet_chs_tau4"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau4") );

    // Energy correlation function variables with PUPPI N2 (for W/Z/H-tagging) and N3 (for top-tagging) with beta=1 and beta=2
    store("fatjet_puppi_n2b1",fatjet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2") );
    store("fatjet_puppi_n3b1",fatjet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3") );
    store("fatjet_puppi_n2b2",fatjet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2") );
    store("fatjet_puppi_n3b2",fatjet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3") );

    // Softdrop and pruned Mass
    store("fatjet_puppi_softdrop_mass",fatjet->userFloat("ak8PFJetsPuppiSoftDropMass")                    );
    store("fatjet_chs_softdrop_mass"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") );
    store("fatjet_chs_pruned_mass"    ,fatjet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   );


    // ------------------------------
    // Flavor Tagging
    store("fatjet_bTag_CMVA"            ,fatjet->bDiscriminator("pfCombinedMVAV2BJetTags")                      );
    store("fatjet_bTag_CSVV2"           ,fatjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    store("fatjet_bTag_Hbb"             ,fatjet->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")    );

    store("fatjet_DeepCSV_b"            ,fatjet->bDiscriminator("pfDeepCSVJetTags:probb")         );
    store("fatjet_DeepCSV_bb"           ,fatjet->bDiscriminator("pfDeepCSVJetTags:probbb")        );
    store("fatjet_DeepCSV_c"            ,fatjet->bDiscriminator("pfDeepCSVJetTags:probc")         );
    store("fatjet_DeepCSV_udsg"         ,fatjet->bDiscriminator("pfDeepCSVJetTags:probudsg")      );
    store("fatjet_DeepCSV_bTag"         ,fatjet->bDiscriminator("pfDeepCSVJetTags:probb")
                                        +fatjet->bDiscriminator("pfDeepCSVJetTags:probbb")        );

    store("fatjet_DeepFlavour_b"        ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probb")     );
    store("fatjet_DeepFlavour_bb"       ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probbb")    );
    store("fatjet_DeepFlavour_lepb"     ,fatjet->bDiscriminator("pfDeepFlavourJetTags:problepb")  );
    store("fatjet_DeepFlavour_c"        ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probc")     );
    store("fatjet_DeepFlavour_uds"      ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probuds")   );
    store("fatjet_DeepFlavour_g"        ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probg")     );
    store("fatjet_DeepFlavour_bTag"     ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probb")
                                        +fatjet->bDiscriminator("pfDeepFlavourJetTags:probbb")    );
    store("fatjet_DeepFlavour_bTagLepb" ,fatjet->bDiscriminator("pfDeepFlavourJetTags:probb")
                                        +fatjet->bDiscriminator("pfDeepFlavourJetTags:probbb")
                                        +fatjet->bDiscriminator("pfDeepFlavourJetTags:problepb")  );

    store("fatjet_DeepDouble_BvL_QCD"   ,fatjet->bDiscriminator("pfDeepDoubleBvLJetTags:probQCD") );
    store("fatjet_DeepDouble_BvL_Hbb"   ,fatjet->bDiscriminator("pfDeepDoubleBvLJetTags:probHbb") );
    store("fatjet_DeepDouble_CvL_QCD"   ,fatjet->bDiscriminator("pfDeepDoubleCvLJetTags:probQCD") );
    store("fatjet_DeepDouble_CvL_Hcc"   ,fatjet->bDiscriminator("pfDeepDoubleCvLJetTags:probHcc") );
    store("fatjet_DeepDouble_CvB_Hbb"   ,fatjet->bDiscriminator("pfDeepDoubleCvBJetTags:probHbb") );
    store("fatjet_DeepDouble_CvB_Hcc"   ,fatjet->bDiscriminator("pfDeepDoubleCvBJetTags:probHcc") );
    store("fatjet_DeepDoubleMD_BvL_QCD" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probQCD") );
    store("fatjet_DeepDoubleMD_BvL_Hbb" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb") );
    store("fatjet_DeepDoubleMD_CvL_QCD" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probQCD") );
    store("fatjet_DeepDoubleMD_CvL_Hcc" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probHcc") );
    store("fatjet_DeepDoubleMD_CvB_Hbb" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHbb") );
    store("fatjet_DeepDoubleMD_CvB_Hcc" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHcc") );


    // ------------------------------
    // Nominal DeepAK8 Tagger

      // aggregated discriminants
      store("fatjet_DeepBoosted_TvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD")   );
      store("fatjet_DeepBoosted_WvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD")   );
      store("fatjet_DeepBoosted_ZvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD")   );
      store("fatjet_DeepBoosted_ZbbvsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD") );
      store("fatjet_DeepBoosted_HbbvsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD") );
      store("fatjet_DeepBoosted_H4qcsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD") );

      // raw discriminants https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoBTag/MXNet/python/pfDeepBoostedDiscriminatorsJetTags_cfi.py
      store("fatjet_DeepBoosted_Tbcq"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbcq")      );
      store("fatjet_DeepBoosted_Tbqq"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbqq")      );
      store("fatjet_DeepBoosted_Tbc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbc")       );
      store("fatjet_DeepBoosted_Tbq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbq")       );
      store("fatjet_DeepBoosted_Wcq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWcq")       );
      store("fatjet_DeepBoosted_Wqq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWqq")       );
      store("fatjet_DeepBoosted_Zbb"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZbb")       );
      store("fatjet_DeepBoosted_Zcc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZcc")       );
      store("fatjet_DeepBoosted_Zqq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZqq")       );
      store("fatjet_DeepBoosted_Hbb"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHbb")       );
      store("fatjet_DeepBoosted_Hcc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHcc")       );
      store("fatjet_DeepBoosted_H4q"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_DeepBoosted_QCDbb"    ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDbb")     );
      store("fatjet_DeepBoosted_QCDcc"    ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDcc")     );
      store("fatjet_DeepBoosted_QCDb"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDb")      );
      store("fatjet_DeepBoosted_QCDc"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDc")      );
      store("fatjet_DeepBoosted_QCDothers",fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDothers") );

      // raw discriminants -- combined
      store("fatjet_DeepBoosted_T"
        ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbcq")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbqq")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbc")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbq")
      );
      store("fatjet_DeepBoosted_W"
        ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWcq")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probWqq")
      );
      store("fatjet_DeepBoosted_Z"
        ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZbb")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probZcc")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probZqq")
      );
      store("fatjet_DeepBoosted_H"
        ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHbb")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probHcc")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probHqqqq")
      );
      store("fatjet_DeepBoosted_QCD"
        ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDbb")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDcc")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDb")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDc")
        +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDothers")
      );

    // ------------------------------
    // Mass-Decorrelated DeepAK8 Tagger

      // aggregated discriminants https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoBTag/MXNet/python/pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags_cfi.py
      store("fatjet_DeepBoostedMD_TvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD")    );
      store("fatjet_DeepBoostedMD_WvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD")    );
      store("fatjet_DeepBoostedMD_ZvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD")    );
      store("fatjet_DeepBoostedMD_ZbbvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD")  );
      store("fatjet_DeepBoostedMD_HbbvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD")  );
      store("fatjet_DeepBoostedMD_H4qvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD")  );
      store("fatjet_DeepBoostedMD_ZHbbvsQCD",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD") );
      store("fatjet_DeepBoostedMD_ZHccvsQCD",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD") );
      store("fatjet_DeepBoostedMD_bbvsLight",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight") );
      store("fatjet_DeepBoostedMD_ccvsLight",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight") );

      // raw discriminants
      store("fatjet_DeepBoostedMD_Tbcq"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")      );
      store("fatjet_DeepBoostedMD_Tbqq"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq")      );
      store("fatjet_DeepBoostedMD_Tbc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbc")       );
      store("fatjet_DeepBoostedMD_Tbq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbq")       );
      store("fatjet_DeepBoostedMD_Wcq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")       );
      store("fatjet_DeepBoostedMD_Wqq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq")       );
      store("fatjet_DeepBoostedMD_Zbb"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")       );
      store("fatjet_DeepBoostedMD_Zcc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")       );
      store("fatjet_DeepBoostedMD_Zqq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq")       );
      store("fatjet_DeepBoostedMD_Hbb"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")       );
      store("fatjet_DeepBoostedMD_Hcc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")       );
      store("fatjet_DeepBoostedMD_H4q"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_DeepBoostedMD_QCDbb"    ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")     );
      store("fatjet_DeepBoostedMD_QCDcc"    ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")     );
      store("fatjet_DeepBoostedMD_QCDb"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")      );
      store("fatjet_DeepBoostedMD_QCDc"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")      );
      store("fatjet_DeepBoostedMD_QCDothers",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers") );

      // raw discriminants -- combined
      store("fatjet_DeepBoostedMD_T"
        ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbc")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbq")
      );
      store("fatjet_DeepBoostedMD_W"
        ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq")
      );
      store("fatjet_DeepBoostedMD_Z"
        ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq")
      );
      store("fatjet_DeepBoostedMD_H"
        ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq")
      );
      store("fatjet_DeepBoostedMD_QCD"
        ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")
        +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers")
      );


    if (isMC_) {
      // ------------------------------
      // Smearing
        jecUnc.setJetEta(FatJetSmearedHandle_->at(i).eta() );
        jecUnc.setJetPt(FatJetSmearedHandle_ ->at(i).pt()  );
        double unc = jecUnc.getUncertainty(true);
        store("fatjet_SmearedJetEnUp_pt"   ,(1+unc)*FatJetSmearedHandle_  ->at(i).pt() );
        store("fatjet_SmearedJetEnDown_pt" ,(1-unc)*FatJetSmearedHandle_  ->at(i).pt() );
        store("fatjet_Smeared_pt"          ,FatJetSmearedHandle_          ->at(i).pt() );
        store("fatjet_SmearedJetResUp_pt"  ,FatJetSmearedJetResUpHandle_  ->at(i).pt() );
        store("fatjet_SmearedJetResDown_pt",FatJetSmearedJetResDownHandle_->at(i).pt() );

    }

  }

}

void IIHEModuleFatJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleFatJet::beginEvent(){}
void IIHEModuleFatJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleFatJet::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleFatJet);
