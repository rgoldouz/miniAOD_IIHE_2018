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
  // puppi, chs
  addBranch("fatjet_puppi_tau1" );
  addBranch("fatjet_puppi_tau2" );
  addBranch("fatjet_puppi_tau3" );
  addBranch("fatjet_puppi_tau4" );
  addBranch("fatjet_chs_tau1"   );
  addBranch("fatjet_chs_tau2"   );
  addBranch("fatjet_chs_tau3"   );
  addBranch("fatjet_chs_tau4"   );
  addBranch("fatjet_puppi_n2b1" );
  addBranch("fatjet_puppi_n3b1" );
  addBranch("fatjet_puppi_n2b2" );
  addBranch("fatjet_puppi_n3b2" );
  addBranch("fatjet_puppi_softdrop_mass");
  addBranch("fatjet_chs_softdrop_mass"  );
  addBranch("fatjet_chs_pruned_mass"    );
  // b Tagger
  addBranch("fatjet_bTag_CMVA"     );
  addBranch("fatjet_bTag_DeepB"    );
  addBranch("fatjet_bTag_CSVV2"    );
  addBranch("fatjet_bTag_Hbb"      );
  addBranch("fatjet_bTag_DDBvL"    ); // btagDDBvL_noMD
  addBranch("fatjet_bTag_DDCvL"    ); // btagDDCvL_noMD
  addBranch("fatjet_bTag_DDCvB"    ); // btagDDCvB_noMD
  addBranch("fatjet_bTagMD_DDBvL"  ); // btagDDBvL
  addBranch("fatjet_bTagMD_DDCvL"  ); // btagDDCvL
  addBranch("fatjet_bTagMD_DDCvB"  ); // btagDDCvB

  // DeepAK8 Nominal Tagger
  // // Aggregated
  addBranch("fatjet_deepTag_TvsQCD"    );
  addBranch("fatjet_deepTag_WvsQCD"    );
  addBranch("fatjet_deepTag_ZvsQCD"    );
  addBranch("fatjet_deepTag_ZbbvsQCD"  );
  addBranch("fatjet_deepTag_HbbvsQCD"  );
  addBranch("fatjet_deepTag_H4qcsQCD"  );
  // // Raw
  addBranch("fatjet_deepTag_Tbcq"      );
  addBranch("fatjet_deepTag_Tbqq"      );
  addBranch("fatjet_deepTag_Tbc"       );
  addBranch("fatjet_deepTag_Tbq"       );
  addBranch("fatjet_deepTag_Wcq"       );
  addBranch("fatjet_deepTag_Wqq"       );
  addBranch("fatjet_deepTag_Zbb"       );
  addBranch("fatjet_deepTag_Zcc"       );
  addBranch("fatjet_deepTag_Zqq"       );
  addBranch("fatjet_deepTag_Hbb"       );
  addBranch("fatjet_deepTag_Hcc"       );
  addBranch("fatjet_deepTag_H4q"       );
  addBranch("fatjet_deepTag_QCDbb"     );
  addBranch("fatjet_deepTag_QCDcc"     );
  addBranch("fatjet_deepTag_QCDb"      );
  addBranch("fatjet_deepTag_QCDc"      );
  addBranch("fatjet_deepTag_QCDothers" );
  // // Raw Combined
  addBranch("fatjet_deepTag_T"         );
  addBranch("fatjet_deepTag_W"         );
  addBranch("fatjet_deepTag_Z"         );
  addBranch("fatjet_deepTag_H"         );
  addBranch("fatjet_deepTag_QCD"       );

  // DeepAK8 Mass-Decorrelated Tagger
  // // Aggregated
  addBranch("fatjet_deepTagMD_TvsQCD"    );
  addBranch("fatjet_deepTagMD_WvsQCD"    );
  addBranch("fatjet_deepTagMD_ZvsQCD"    );
  addBranch("fatjet_deepTagMD_ZbbvsQCD"  );
  addBranch("fatjet_deepTagMD_HbbvsQCD"  );
  addBranch("fatjet_deepTagMD_H4qvsQCD"  );
  addBranch("fatjet_deepTagMD_ZHbbvsQCD" );
  addBranch("fatjet_deepTagMD_ZHccvsQCD" );
  addBranch("fatjet_deepTagMD_bbvsLight" );
  addBranch("fatjet_deepTagMD_ccvsLight" );
  // // Raw
  addBranch("fatjet_deepTagMD_Tbcq"      );
  addBranch("fatjet_deepTagMD_Tbqq"      );
  addBranch("fatjet_deepTagMD_Tbc"       );
  addBranch("fatjet_deepTagMD_Tbq"       );
  addBranch("fatjet_deepTagMD_Wcq"       );
  addBranch("fatjet_deepTagMD_Wqq"       );
  addBranch("fatjet_deepTagMD_Zbb"       );
  addBranch("fatjet_deepTagMD_Zcc"       );
  addBranch("fatjet_deepTagMD_Zqq"       );
  addBranch("fatjet_deepTagMD_Hbb"       );
  addBranch("fatjet_deepTagMD_Hcc"       );
  addBranch("fatjet_deepTagMD_H4q"       );
  addBranch("fatjet_deepTagMD_QCDbb"     );
  addBranch("fatjet_deepTagMD_QCDcc"     );
  addBranch("fatjet_deepTagMD_QCDb"      );
  addBranch("fatjet_deepTagMD_QCDc"      );
  addBranch("fatjet_deepTagMD_QCDothers" );
  // // Raw Combined
  addBranch("fatjet_deepTagMD_T"         );
  addBranch("fatjet_deepTagMD_W"         );
  addBranch("fatjet_deepTagMD_Z"         );
  addBranch("fatjet_deepTagMD_H"         );
  addBranch("fatjet_deepTagMD_QCD"       );
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

    // N-subjettiness
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

    // Mass
    store("fatjet_puppi_softdrop_mass",fatjet->userFloat("ak8PFJetsPuppiSoftDropMass")                    );
    store("fatjet_chs_softdrop_mass"  ,fatjet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") );
    store("fatjet_chs_pruned_mass"    ,fatjet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   );


    // ------------------------------
    // Flavor Tagging
    store("fatjet_bTag_CMVA"    ,fatjet->bDiscriminator("pfCombinedMVAV2BJetTags") );
    store("fatjet_bTag_DeepB"   ,fatjet->bDiscriminator("pfDeepCSVJetTags:probb")
                                +fatjet->bDiscriminator("pfDeepCSVJetTags:probbb") );
    store("fatjet_bTag_CSVV2"   ,fatjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    store("fatjet_bTag_Hbb"     ,fatjet->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")    );
    store("fatjet_bTag_DDBvL"   ,fatjet->bDiscriminator("pfDeepDoubleBvLJetTags:probHbb") ); // btagDDBvL_noMD
    store("fatjet_bTag_DDCvL"   ,fatjet->bDiscriminator("pfDeepDoubleCvLJetTags:probHcc") ); // btagDDCvL_noMD
    store("fatjet_bTag_DDCvB"   ,fatjet->bDiscriminator("pfDeepDoubleCvBJetTags:probHcc") ); // btagDDCvB_noMD
    store("fatjet_bTagMD_DDBvL" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb") ); // btagDDBvL
    store("fatjet_bTagMD_DDCvL" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probHcc") ); // btagDDCvL
    store("fatjet_bTagMD_DDCvB" ,fatjet->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHcc") ); // btagDDCvB

    // ------------------------------
    // Nominal DeepAK8 Tagger

      // aggregated discriminants
      store("fatjet_deepTag_TvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD")   );
      store("fatjet_deepTag_WvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD")   );
      store("fatjet_deepTag_ZvsQCD"   ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD")   );
      store("fatjet_deepTag_ZbbvsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD") );
      store("fatjet_deepTag_HbbvsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD") );
      store("fatjet_deepTag_H4qcsQCD" ,fatjet->bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD") );

      // raw discriminants https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoBTag/MXNet/python/pfDeepBoostedDiscriminatorsJetTags_cfi.py
      store("fatjet_deepTag_Tbcq"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbcq")      );
      store("fatjet_deepTag_Tbqq"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbqq")      );
      store("fatjet_deepTag_Tbc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbc")       );
      store("fatjet_deepTag_Tbq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbq")       );
      store("fatjet_deepTag_Wcq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWcq")       );
      store("fatjet_deepTag_Wqq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWqq")       );
      store("fatjet_deepTag_Zbb"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZbb")       );
      store("fatjet_deepTag_Zcc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZcc")       );
      store("fatjet_deepTag_Zqq"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZqq")       );
      store("fatjet_deepTag_Hbb"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHbb")       );
      store("fatjet_deepTag_Hcc"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHcc")       );
      store("fatjet_deepTag_H4q"      ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_deepTag_QCDbb"    ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDbb")     );
      store("fatjet_deepTag_QCDcc"    ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDcc")     );
      store("fatjet_deepTag_QCDb"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDb")      );
      store("fatjet_deepTag_QCDc"     ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDc")      );
      store("fatjet_deepTag_QCDothers",fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDothers") );

      // raw discriminants -- combined
      store("fatjet_deepTag_T"  ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbcq")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbqq")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbc")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probTbq")       );
      store("fatjet_deepTag_W"  ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probWcq")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probWqq")       );
      store("fatjet_deepTag_Z"  ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probZbb")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probZcc")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probZqq")       );
      store("fatjet_deepTag_H"  ,fatjet->bDiscriminator("pfDeepBoostedJetTags:probHbb")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probHcc")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_deepTag_QCD",fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDbb")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDcc")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDb")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDc")
                                +fatjet->bDiscriminator("pfDeepBoostedJetTags:probQCDothers") );

    // ------------------------------
    // Mass-Decorrelated DeepAK8 Tagger

      // aggregated discriminants https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoBTag/MXNet/python/pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags_cfi.py
      store("fatjet_deepTagMD_TvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD")    );
      store("fatjet_deepTagMD_WvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD")    );
      store("fatjet_deepTagMD_ZvsQCD"   ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD")    );
      store("fatjet_deepTagMD_ZbbvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD")  );
      store("fatjet_deepTagMD_HbbvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD")  );
      store("fatjet_deepTagMD_H4qvsQCD" ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD")  );
      store("fatjet_deepTagMD_ZHbbvsQCD",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD") );
      store("fatjet_deepTagMD_ZHccvsQCD",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD") );
      store("fatjet_deepTagMD_bbvsLight",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight") );
      store("fatjet_deepTagMD_ccvsLight",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight") );

      // raw discriminants
      store("fatjet_deepTagMD_Tbcq"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")      );
      store("fatjet_deepTagMD_Tbqq"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq")      );
      store("fatjet_deepTagMD_Tbc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbc")       );
      store("fatjet_deepTagMD_Tbq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbq")       );
      store("fatjet_deepTagMD_Wcq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")       );
      store("fatjet_deepTagMD_Wqq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq")       );
      store("fatjet_deepTagMD_Zbb"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")       );
      store("fatjet_deepTagMD_Zcc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")       );
      store("fatjet_deepTagMD_Zqq"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq")       );
      store("fatjet_deepTagMD_Hbb"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")       );
      store("fatjet_deepTagMD_Hcc"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")       );
      store("fatjet_deepTagMD_H4q"      ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_deepTagMD_QCDbb"    ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")     );
      store("fatjet_deepTagMD_QCDcc"    ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")     );
      store("fatjet_deepTagMD_QCDb"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")      );
      store("fatjet_deepTagMD_QCDc"     ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")      );
      store("fatjet_deepTagMD_QCDothers",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers") );

      // raw discriminants -- combined
      store("fatjet_deepTagMD_T"  ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbc")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbq")       );
      store("fatjet_deepTagMD_W"  ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq")       );
      store("fatjet_deepTagMD_Z"  ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq")       );
      store("fatjet_deepTagMD_H"  ,fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq")     );
      store("fatjet_deepTagMD_QCD",fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")
                                  +fatjet->bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers") );


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
