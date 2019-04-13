import FWCore.ParameterSet.Config as cms
import getpass, os
#see https://github.com/cms-sw/cmssw/blob/2e3e3823af73ee521584fc9888ee77ec5adbc8d1/PhysicsTools/NanoAOD/python/nano_cff.py#L60
btagWeightTable = cms.EDProducer("BTagSFProducer",
    src = cms.InputTag("slimmedJets"),
    cut = cms.string("pt > 30. && abs(eta) < 2.5"),
    discNames = cms.vstring(
        "pfCombinedInclusiveSecondaryVertexV2BJetTags",
        "pfDeepCSVJetTags:probb+pfDeepCSVJetTags:probbb",       #if multiple MiniAOD branches need to be summed up (e.g., DeepCSV b+bb), separate them using '+' delimiter
    ),
    discShortNames = cms.vstring(
        "CSVV2",
        "DeepCSVB",
    ),
    weightFiles = cms.vstring(                                  #default settings are for 2017 94X. toModify function is called later for other eras.
        "unavailable",                    
        "unavailable"                                           #if SFs for an algorithm in an era is unavailable, the corresponding branch will not be stored
    ),
    operatingPoints = cms.vstring("3","3"),                 #loose = 0, medium = 1, tight = 2, reshaping = 3
    measurementTypesB = cms.vstring("iterativefit","iterativefit"),     #e.g. "comb", "incl", "ttbar", "iterativefit"
    measurementTypesC = cms.vstring("iterativefit","iterativefit"),
    measurementTypesUDSG = cms.vstring("iterativefit","iterativefit"),
    sysTypes = cms.vstring("central","central")  #e.g. central, plus, minus, plus_JEC, plus_JER, ...
)
