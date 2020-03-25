##########################
#                        #
#   tW                   #
#                        #
##########################
#dataset = {
#"SingleMuon_runB":      ["/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"],
#"SingleMuon_runC":      ["/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"],
#"SingleMuon_runD":      ["/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"],
#"SingleMuon_runE":      ["/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"],
#"SingleMuon_runF":      ["/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"],
#
#"SingleElectron_runB":  ["/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD"],
#"SingleElectron_runC":  ["/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD"],
#"SingleElectron_runD":  ["/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD"],
#"SingleElectron_runE":  ["/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD"],
#"SingleElectron_runF":  ["/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD"],
#
#"DoubleEG_runB":        ["/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD"],
#"DoubleEG_runC":        ["/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD"],
#"DoubleEG_runD":        ["/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD"],
#"DoubleEG_runE":        ["/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD"],
#"DoubleEG_runF":        ["/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD"],
#
#"DoubleMuon_runB":      ["/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD"],
#"DoubleMuon_runC":      ["/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD"],
#"DoubleMuon_runD":      ["/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD"],
#"DoubleMuon_runE":      ["/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD"],
#"DoubleMuon_runF":      ["/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD"],
#
#"MuonEG_runB":  ["/MuonEG/Run2017B-31Mar2018-v1/MINIAOD"],
#"MuonEG_runC":  ["/MuonEG/Run2017C-31Mar2018-v1/MINIAOD"],
#"MuonEG_runD":  ["/MuonEG/Run2017D-31Mar2018-v1/MINIAOD"],
#"MuonEG_runE":  ["/MuonEG/Run2017E-31Mar2018-v1/MINIAOD"],
#"MuonEG_runF":  ["/MuonEG/Run2017F-31Mar2018-v1/MINIAOD"],
#
#"SinglePhoton_runB":    ["/SinglePhoton/Run2017B-31Mar2018-v1/MINIAOD"],
#"SinglePhoton_runC":    ["/SinglePhoton/Run2017C-31Mar2018-v1/MINIAOD"],
#"SinglePhoton_runD":    ["/SinglePhoton/Run2017D-31Mar2018-v1/MINIAOD"],
#"SinglePhoton_runE":    ["/SinglePhoton/Run2017E-31Mar2018-v1/MINIAOD"],
#"SinglePhoton_runF":    ["/SinglePhoton/Run2017F-31Mar2018-v1/MINIAOD"],
#}

dataset = {
#"TT_Mtt_700to1000":      ["/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"],
#"TT_Mtt_1000toInf":      ["/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"],
#"TT":      ["/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"]
#"TT_MG":      ["/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"]
"TT_TuneCP5":      ["/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"]
}

nfiles = {
}

filesPerJob = {
}

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    def submit(config):
        res = crabCommand('submit', config = config)

    from CRABClient.UserUtilities import config
    config = config()

    name = 'TOPptSamples'
    config.General.workArea = 'crab_'+name
    config.General.transferLogs = False
    config.General.transferOutputs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'IIHE.py'
    config.Data.inputDBS = 'global'
#    config.Data.splitting = 'LumiBased'
    config.Data.splitting = 'FileBased'
#    config.Data.lumiMask = '/user/xgao/IIHENtuples/2017_TOP_prefire/CMSSW_10_2_13/src/UserCode/IIHETree/test/crab_20190516_data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    config.Site.storageSite = 'T2_BE_IIHE'
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pyCfgParams = ['DataProcessing=mc2016']
    config.JobType.inputFiles   = ['data']
#    config.Data.splitting = 'Automatic'
    config.Data.unitsPerJob = 1
    config.Data.publication = False
    config.Site.storageSite = 'T2_BE_IIHE'

    for sample in dataset:
        config.General.requestName = sample
        config.Data.inputDataset = dataset[sample][0]
        config.Data.outLFNDirBase = '/store/user/rgoldouz/' + name
#        config.Data.outputDatasetTag = sample
        submit(config)


