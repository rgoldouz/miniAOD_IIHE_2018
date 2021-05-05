data2017_samples = {}
#data2017_samples["2017_B_SinglePhoton"]=[    ["/SinglePhoton/Run2017B-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_C_SinglePhoton"]=[    ["/SinglePhoton/Run2017C-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_D_SinglePhoton"]=[    ["/SinglePhoton/Run2017D-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_E_SinglePhoton"]=[    ["/SinglePhoton/Run2017E-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_F_SinglePhoton"]=[    ["/SinglePhoton/Run2017F-31Mar2018-v1/MINIAOD"],    "data2017"]

data2017_samples["2017_B_SingleMuon"]=[    ["/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_C_SingleMuon"]=[    ["/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_D_SingleMuon"]=[    ["/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_E_SingleMuon"]=[    ["/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"],    "data2017"]
data2017_samples["2017_F_SingleMuon"]=[    ["/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"],    "data2017"]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    def submit(config):
        res = crabCommand('submit', config = config)

    from CRABClient.UserUtilities import config
    config = config()

    name = 'ExitedTopSamplesDataJan2021'
    config.General.workArea = 'crab_'+name
    config.General.transferLogs = False
    config.General.transferOutputs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'IIHE.py'
    config.Data.inputDBS = 'global'
#    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'FileBased'
    config.Data.lumiMask = '/afs/cern.ch/work/r/rgoldouz/IIHETree/miniAOD/slc7/CMSSW_10_2_13/src/UserCode/IIHETree/test/data/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
    config.Site.storageSite = 'T3_US_NotreDame'
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pyCfgParams = ['DataProcessing=data2017']
    config.JobType.inputFiles   = ['data']
    config.Data.splitting = 'Automatic'
#    config.Data.unitsPerJob = 40
    config.Data.publication = False

    for key, value in data2017_samples.items():
        config.General.requestName = key
        config.Data.inputDataset = value[0][0]
        config.Data.outLFNDirBase = '/store/user/rgoldouz/' + name
#        config.Data.outputDatasetTag = sample
        submit(config)


