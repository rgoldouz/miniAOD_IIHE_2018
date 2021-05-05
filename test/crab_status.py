import os

mc2017_samples = {}
mc2017_samples["2017_GJetsHT100To200"]=[    ["/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_GJetsHT200To400"]=[    ["/GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_GJetsHT400To600"]=[    ["/GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_GJetsHT600ToInf"]=[    ["/GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
        
mc2017_samples["2017_QCDHT100to200_1"]=[    ["/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT100to200_2"]=[    ["/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT200to300"]=[    ["/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT300to500"]=[    ["/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT500to700"]=[    ["/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT700to1000"]=[    ["/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT1000to1500"]=[    ["/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT1500to2000"]=[    ["/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_QCDHT2000toInf"]=[    ["/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
        
mc2017_samples["2017_DYJetsToLLM50"]=[    ["/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_DYJetsToLLM10to50"]=[    ["/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_WJetsToLNu"]=[    ["/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"],    "mc2017"]
        
mc2017_samples["2017_WGPtG40to130"]=[    ["/WGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_WGPtG130"]=[    ["/WGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
       
mc2017_samples["2017_tt"]=[    ["/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_ttG"]=[    ["/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]

mc2017_samples["2017_ttToSemiLeptonic"]=[    ["/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_STtchtop"]=[    ["/ST_t-channel_top_4f_InclusiveDecays_TuneCP5up_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_STtchatop"]=[    ["/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_STtwchtop"]=[    ["/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]
mc2017_samples["2017_STtwchatop"]=[    ["/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"],    "mc2017"]

data2017_samples = {}
data2017_samples["2017_B_SinglePhoton"]=[    ["/SinglePhoton/Run2017B-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_C_SinglePhoton"]=[    ["/SinglePhoton/Run2017C-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_D_SinglePhoton"]=[    ["/SinglePhoton/Run2017D-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_E_SinglePhoton"]=[    ["/SinglePhoton/Run2017E-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_F_SinglePhoton"]=[    ["/SinglePhoton/Run2017F-31Mar2018-v1/MINIAOD"],    "data2017"]

data2017_samples["2017_B_SingleMuon"]=[    ["/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_C_SingleMuon"]=[    ["/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_D_SingleMuon"]=[    ["/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_E_SingleMuon"]=[    ["/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_F_SingleMuon"]=[    ["/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"],    "data2017"]

#for key, value in mc2017_samples.items():
#    os.system("crab status -d crab_ExitedTopSamplesMCJan2021/crab_" + key)
#    os.system("crab kill -d crab_ExitedTopSamplesMCJan2021/crab_" + key)
for key, value in data2017_samples.items():
    os.system("crab status -d crab_ExitedTopSamplesDataJan2021/crab_" + key)
#    os.system("crab resubmit -d crab_ExitedTopSamplesDataJan2021/crab_" + key)

