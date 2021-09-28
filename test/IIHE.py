# In order to run the code for MC/DATA on lxplus
#[cmsRun IIHE.py grid=False DataProcessing="mc2016" DataFormat="mc" dataset="RunIISummer16MiniAODv2" sample="ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1" address="MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000" file="0E4CD60A-0FC3-E611-BCB5-0CC47A7C3420.root"  ]
#root://eoscms//cms/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root

# In order to run the code for DATA on lxplus
#[cmsRun IIHE.py grid=False DataProcessing="rerecodata" DataFormat="data" dataset="Run2016D" sample="DoubleEG" address="MINIAOD/03Feb2017-v1/100000" file="CC7F39AC-F5EA-E611-A125-0CC47A13CDB0.root"  ]
#root://eoscms//cms/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6E918C9-87A6-E511-B3D3-0CC47A4D76B2.root

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy
import os
from os import environ


options = opts.VarParsing ("analysis")
options.register("sample",
                 "ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Sample to analyze")
options.register("address",
                 "MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "address of sample in eos like: MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000")
options.register("file",
                 "0E4CD60A-0FC3-E611-BCB5-0CC47A7C3420.root",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "file to analyze")
options.register("DataProcessing",
                 "mc2017",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Data processing types. Options are:mc2016,rerecodata,promptdata")
options.register("dataset",
                 "RunIISummer16MiniAODv2",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "datasets to analyze: RunIISummer16MiniAODv2, Run2016B-H")
options.register("DataFormat",
                 "mc",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Data format. Options are:mc,data")
options.register('grid',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'If you run on grid or localy on eos')
options.parseArguments()

###########################################################################################
#                                      Global tags                                       #
##########################################################################################
globalTag = "80"

if options.DataProcessing == "data2016":
    globalTag = "94X_dataRun2_v10"
if options.DataProcessing == "mc2016":
    globalTag = "94X_mcRun2_asymptotic_v3"

if options.DataProcessing == "data2017":
    globalTag = "102X_dataRun2_v13"
#    globalTag = "94X_dataRun2_v11"
if options.DataProcessing == "mc2017":
    globalTag = "94X_mc2017_realistic_v17"

if options.DataProcessing == "data2018ABC":
    globalTag = "102X_dataRun2_v12"
if options.DataProcessing == "data2018D":
    globalTag = "102X_dataRun2_Prompt_v15"
if options.DataProcessing == "mc2018":
    globalTag = "102X_upgrade2018_realistic_v20"

if globalTag == "80":
    print '*****ERROR**** you are using wrong global tag so lets stop running'
    quit()
##########################################################################################
#                                  Start the sequences                                   #
##########################################################################################
process = cms.Process("IIHEAnalysis")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag
##process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
##########################################################################################
#                                         Files                                          #
##########################################################################################
if options.DataFormat == "mc":
  path = "root://eoscms//eos/cms/store/"+ options.DataFormat + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

if options.DataFormat == "data":
  path = "root://eoscms//eos/cms/store/"+ options.DataFormat + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
#    lumisToProcess = cms.untracked.VLuminosityBlockRange('319756:1567','319337:56'),
#    lumisToProcess = cms.untracked.VLuminosityBlockRange('1:427', '1:428'),
#    eventsToProcess = cms.untracked.VEventRange('297483:182:250326304')
#    skipEvents=cms.untracked.uint32(8000)
)
#process.source.fileNames.append( "process.source.fileNames.append( "/store/mc/RunIIFall17MiniAODv2/TprimeBToBW_M-1200_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/30000/CA1B02AB-4433-E911-8B8E-FA163E27F7AC.root")
#process.source.fileNames.append( "/store/mc/RunIIFall17MiniAODv2/ZprimeToTTJet_M1000_TuneCP2_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/3C2F88E6-572E-E911-80DD-0242AC130002.root")
#process.source.fileNames.append( "file:3C2F88E6-572E-E911-80DD-0242AC130002.root")
#process.source.fileNames.append( "/store/mc/RunIIFall17MiniAODv2/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/120000/5047A56F-B9EB-E811-8CAE-D067E5F91B8A.root")
#process.source.fileNames.append("/store/mc/RunIIFall17MiniAODv2/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/1EF6C9F2-2D5D-E811-9D2A-1C6A7A263523.root")
process.source.fileNames.append("/store/mc/RunIIFall17MiniAODv2/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/1EF6C9F2-2D5D-E811-9D2A-1C6A7A263523.root")
#/store/mc/RunIIFall17MiniAODv2/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/110000/885ECDE7-10B9-E811-B764-549F35AD8C0A.root")

#process.source.fileNames.append("file:pickevents.root")

filename_out = "outfile.root"
if options.DataFormat == "mc" and not options.grid:
#  filename_out = "file:/tmp/output_%s" % (options.sample + "_" + options.file)
  filename_out = "outfile_MC.root"
if options.DataFormat == "data" and not options.grid:
  filename_out = "outfile_Data.root"

#filename_out = "outfile.root"
#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string(filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename_out) )


##########################################################################################
#                                   NEW tau id disc.                                     #
##########################################################################################
from UserCode.IIHETree.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
    debug=True,
    toKeep = ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()


##########################################################################################
#                                     Jet Collections                                    #
##########################################################################################
process.load("UserCode.IIHETree.JetProcessing")
from UserCode.IIHETree.JetProcessing import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

##############################
# AK4 CHS Jets with Corrections
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
    pvSource  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource  = cms.InputTag('slimmedSecondaryVertices'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
)
process.AK4CHSJetSequence = cms.Sequence(
    process.patJetCorrFactorsUpdatedJEC *
    process.updatedPatJetsUpdatedJEC
)

##############################
# AK4 Puppi Jets with Corrections
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsPuppi'),
    pvSource  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource  = cms.InputTag('slimmedSecondaryVertices'),
   labelName = 'PuppiJEC',
   jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
#   outputModules = 'noOutput'
)
process.AK4PuppiJetSequence = cms.Sequence(
    process.patJetCorrFactorsPuppiJEC *
    process.updatedPatJetsPuppiJEC
)

###############################
#AK8 Puppi Jets with Corrections
from JMEAnalysis.JetToolbox.jetToolbox_cff import *
from RecoBTag.MXNet.pfDeepBoostedJet_cff   import _pfDeepBoostedJetTagsAll

runOnMC = False
if (options.DataFormat == "mc"):
    runOnMC = True
tagDiscriminatorsDeepAK8 = [
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    'pfDeepCSVJetTags:probb',
    'pfDeepCSVJetTags:probbb',
    'pfDeepCSVJetTags:probc',
    'pfDeepCSVJetTags:probudsg',
    'pfDeepFlavourJetTags:probb',
    'pfDeepFlavourJetTags:probbb',
    'pfDeepFlavourJetTags:problepb',
    'pfDeepFlavourJetTags:probc',
    'pfDeepFlavourJetTags:probuds',
    'pfDeepFlavourJetTags:probg',
    'pfDeepDoubleBvLJetTags:probQCD',
    'pfDeepDoubleBvLJetTags:probHbb',
    'pfDeepDoubleCvLJetTags:probQCD',
    'pfDeepDoubleCvLJetTags:probHcc',
    'pfDeepDoubleCvBJetTags:probHbb',
    'pfDeepDoubleCvBJetTags:probHcc',
    'pfMassIndependentDeepDoubleBvLJetTags:probQCD',
    'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
    'pfMassIndependentDeepDoubleCvLJetTags:probQCD',
    'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
    'pfMassIndependentDeepDoubleCvBJetTags:probHbb',
    'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
]
tagDiscriminatorsDeepAK8+=_pfDeepBoostedJetTagsAll

bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    'pfMassIndependentDeepDoubleBvLJetTags:probQCD',
    'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
    'pfDeepCSVJetTags:probb',
    'pfDeepCSVJetTags:probbb',
]
subjetBTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfDeepCSVJetTags:probb',
    'pfDeepCSVJetTags:probbb',
]


updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsAK8'),
    pvSource  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource  = cms.InputTag('slimmedSecondaryVertices'),
    rParam    = 0.8,
    jetCorrections = ('AK8PFPuppi', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual']), 'None'),
    btagDiscriminators = tagDiscriminatorsDeepAK8,
    postfix        = 'DeepAK8',
    printWarning   = False,
)



jetVer = ''
if '2016' in options.DataProcessing:
    jetVer = 'WINTER16'
if '2017' in options.DataProcessing:
    jetVer = 'WINTER17PUPPI'
if '2018' in options.DataProcessing:
    jetVer = 'SUMMER18PUPPI'

process.mylooseJetIdAK8 = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
        version = cms.string(jetVer),
        quality = cms.string('LOOSE'),
    ),
    src = cms.InputTag("selectedUpdatedPatJetsDeepAK8")
)

process.mytightJetIdAK8 = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
        version = cms.string(jetVer),
        quality = cms.string('TIGHT'),
    ),
    src = cms.InputTag("selectedUpdatedPatJetsDeepAK8")
#    src = cms.InputTag("selectedUpdatedPatJetsDeepAK8")
)

process.mytightJetIdLepVetoAK8 = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
        version = cms.string(jetVer),
        quality = cms.string('TIGHTLEPVETO'),
    ),
    src = cms.InputTag("selectedUpdatedPatJetsDeepAK8")
)

process.selectedUpdatedPatJetsDeepAK8WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
     src = cms.InputTag("selectedUpdatedPatJetsDeepAK8"),
     userFloats = cms.PSet(),
     userInts = cms.PSet(
     tightId = cms.InputTag("mytightJetIdAK8"),
     tightIdLepVeto = cms.InputTag("mytightJetIdLepVetoAK8"),
     looseId = cms.InputTag("mylooseJetIdAK8")
     ),
)

process.p10 = cms.Path(
process.mylooseJetIdAK8 *
process.mytightJetIdAK8 *
process.mytightJetIdLepVetoAK8 *
process.selectedUpdatedPatJetsDeepAK8WithUserData
)

###########################################################################################
##                                       Jet Smearing                                     #
###########################################################################################
process.SmearedJetsSequence = cms.Sequence(
    mySmearedJets     *
    mySmearedJetsUp   *
    mySmearedJetsDown
)
process.SmearedPuppiJetsSequence = cms.Sequence(
    mySmearedPuppiJets     *
    mySmearedPuppiJetsUp   *
    mySmearedPuppiJetsDown
)
process.SmearedDeepAK8Sequence = cms.Sequence(
    mySmearedDeepAK8Jets     *
    mySmearedDeepAK8JetsUp   *
    mySmearedDeepAK8JetsDown
)
##########################################################################################
#                                   2018 electron scale smearing                         #
##########################################################################################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

if "2016" in options.DataProcessing:
    setupEgammaPostRecoSeq(process, era='2016-Legacy')
if "2017" in options.DataProcessing:
    setupEgammaPostRecoSeq(process, era='2017-Nov17ReReco')
if "2018" in options.DataProcessing:
    setupEgammaPostRecoSeq(process, era='2018-Prompt')


##########################################################################################
#                                   2018 electron L1 prefireing                         #
##########################################################################################
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

if "2016" in options.DataProcessing:
    process.prefiringweight.DataEra = cms.string("2016BtoH")

print 'L1 prefireing weight is saved for runs  ' + str(process.prefiringweight.DataEra)


##########################################################################################
#                                   2018 MET corrections                                 #
##########################################################################################
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
        process, 
        isData = "data" in options.DataProcessing, # false for MC
        fixEE2017 = "2017" in options.DataProcessing,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = ""
)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336, 872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339, 872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185, 872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694, 872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609, 872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413,872421414,872421031,872423083,872421439]
)
process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal       = baddetEcallist,
    taggingMode      = cms.bool(True),
    debug            = cms.bool(False)
)


##########################################################################################
#                            MY analysis input!                                       ####
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag  = cms.string(globalTag)
process.IIHEAnalysis.isData     = cms.untracked.bool("data" in options.DataProcessing)
process.IIHEAnalysis.isMC       = cms.untracked.bool("mc"   in options.DataProcessing)
#effectivearea file
process.IIHEAnalysis.EAFile = cms.FileInPath("UserCode/IIHETree/test/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt")
# Collections added before the analysis
process.IIHEAnalysis.particleFlowEGammaGSFixedCollection     = cms.InputTag("particleFlowEGammaGSFixed"     ,"dupECALClusters"                )
process.IIHEAnalysis.ecalMultiAndGSGlobalRecHitEBCollection  = cms.InputTag("ecalMultiAndGSGlobalRecHitEB"  ,"dupESClusters"   ,"PAT"         )
process.IIHEAnalysis.METsMuEGCleanCollection                 = cms.InputTag("slimmedMETsMuEGClean"                                            )
process.IIHEAnalysis.discardedMuonCollection                 = cms.InputTag("packedPFCandidatesDiscarded"                                     )
# AK4 CHS Jet Collections
process.IIHEAnalysis.JetCollectionPrecor                     = cms.InputTag("slimmedJets"                                                     )
process.IIHEAnalysis.JetCollection                           = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"                                        )
process.IIHEAnalysis.JetCollectionSmeared                    = cms.InputTag("mySmearedJets"                 ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.JetCollectionSmearedJetResUp            = cms.InputTag("mySmearedJetsUp"               ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.JetCollectionSmearedJetResDown          = cms.InputTag("mySmearedJetsDown"             ,""                ,"IIHEAnalysis")
# AK4 Puppi Jet Collections
process.IIHEAnalysis.PuppiJetCollectionPrecor                = cms.InputTag("slimmedJetsPuppi"                                                )
process.IIHEAnalysis.PuppiJetCollection                      = cms.InputTag("updatedPatJetsPuppiJEC"                                          )
process.IIHEAnalysis.PuppiJetCollectionSmeared               = cms.InputTag("mySmearedPuppiJets"            ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.PuppiJetCollectionSmearedJetResUp       = cms.InputTag("mySmearedPuppiJetsUp"          ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.PuppiJetCollectionSmearedJetResDown     = cms.InputTag("mySmearedPuppiJetsDown"        ,""                ,"IIHEAnalysis")
# FatJet Collections
process.IIHEAnalysis.DeepAK8JetCollection                    = cms.InputTag("selectedUpdatedPatJetsDeepAK8WithUserData",   "",               "IIHEAnalysis")
process.IIHEAnalysis.DeepAK8JetCollectionSmeared             = cms.InputTag("mySmearedDeepAK8Jets"          ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.DeepAK8JetCollectionSmearedJetResUp     = cms.InputTag("mySmearedDeepAK8JetsUp"        ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.DeepAK8JetCollectionSmearedJetResDown   = cms.InputTag("mySmearedDeepAK8JetsDown"      ,""                ,"IIHEAnalysis")
# MET Collections
process.IIHEAnalysis.patPFMetCollection                      = cms.InputTag("patPFMet"                      ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1Collection                    = cms.InputTag("patPFMetT1"                    ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1JetEnDownCollection           = cms.InputTag("patPFMetT1JetEnDown"           ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1JetEnUpCollection             = cms.InputTag("patPFMetT1JetEnUp"             ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1SmearCollection               = cms.InputTag("patPFMetT1Smear"               ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1SmearJetEnDownCollection      = cms.InputTag("patPFMetT1SmearJetEnDown"      ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1SmearJetEnUpCollection        = cms.InputTag("patPFMetT1SmearJetEnUp"        ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1SmearJetResDownCollection     = cms.InputTag("patPFMetT1SmearJetResDown"     ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1SmearJetResUpCollection       = cms.InputTag("patPFMetT1SmearJetResUp"       ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetT1TxyCollection                 = cms.InputTag("patPFMetT1Txy"                 ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.patPFMetFinalCollection                 = cms.InputTag("slimmedMETs"                   ,""                ,"IIHEAnalysis")
# Particle Level Fiducial Collections
process.IIHEAnalysis.particleLevelJetsCollection             = cms.InputTag("pseudoTop"                     ,"jets"            ,"IIHEAnalysis")
process.IIHEAnalysis.particleLevelak1DressedLeptonCollection = cms.InputTag("pseudoTop"                     ,"leptons"         ,"IIHEAnalysis")
process.IIHEAnalysis.particleLevelMETCollection              = cms.InputTag("pseudoTop"                     ,"mets"            ,"IIHEAnalysis")
# Corrected Electron Collection
process.IIHEAnalysis.electronCollection                      = cms.InputTag("slimmedElectrons"              ,""                ,"IIHEAnalysis")
process.IIHEAnalysis.photonCollection                        = cms.InputTag("slimmedPhotons"                ,""                ,"IIHEAnalysis")
#
process.IIHEAnalysis.includeSkimEventsModule      = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeEventModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule             = cms.untracked.bool(True)
process.IIHEAnalysis.includeJetModule             = cms.untracked.bool(True)
#process.IIHEAnalysis.includePuppiJetModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeFatJetModule          = cms.untracked.bool(True)
# process.IIHEAnalysis.includeTauModule             = cms.untracked.bool(True)
process.IIHEAnalysis.includePhotonModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule         = cms.untracked.bool("mc" in options.DataProcessing)
process.IIHEAnalysis.includeLHEWeightModule       = cms.untracked.bool("mc" in options.DataProcessing)
process.IIHEAnalysis.includeAutoAcceptEventModule = cms.untracked.bool(False)
##########################################################################################
#                            Woohoo!  We"re ready to start!                              #
##########################################################################################
if "mc" in options.DataProcessing:
    process.IIHE = cms.Sequence(
    process.egammaPostRecoSeq   *
    process.prefiringweight     *
    process.rerunMvaIsolationSequence *
    process.NewTauIDsEmbedded   *
    process.AK4CHSJetSequence  * 
    process.fullPatMetSequence  *
    process.ecalBadCalibReducedMINIAODFilter * 
    process.SmearedJetsSequence     * 
    process.SmearedDeepAK8Sequence   *
    process.IIHEAnalysis
    )
else:
    process.IIHE = cms.Sequence(
    process.egammaPostRecoSeq   *
    process.rerunMvaIsolationSequence *
    process.NewTauIDsEmbedded   *
    process.AK4CHSJetSequence   *
    process.fullPatMetSequence  *
    process.ecalBadCalibReducedMINIAODFilter *
    process.IIHEAnalysis
    )

process.p1 = cms.Path(process.IIHE)
process.p1.associate(process.patAlgosToolsTask)

#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("EDMiii.root")
#    )

#process.outpath = cms.EndPath(process.out)
##
