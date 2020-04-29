import FWCore.ParameterSet.Config as cms
from RecoBTag.MXNet.pfDeepBoostedJet_cff   import _pfDeepBoostedJetTagsAll

##########################################################################################
#                                   Deep AK8 Jet Tagging                                 #
##########################################################################################
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging

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

# Tag Task for AK8 Jets
# DeepAK8TagTask = cms.Task(
#     pfImpactParameterTagInfosDeepAK8,
#     pfImpactParameterAK8TagInfosDeepAK8,
#     pfInclusiveSecondaryVertexFinderTagInfosDeepAK8,
#     pfInclusiveSecondaryVertexFinderAK8TagInfosDeepAK8,
#     pfBoostedDoubleSVAK8TagInfosDeepAK8,
#     pfBoostedDoubleSecondaryVertexAK8BJetTagsDeepAK8,
#     pfDeepCSVTagInfosDeepAK8,
#     pfDeepCSVJetTagsDeepAK8,
#     pfDeepFlavourTagInfosDeepAK8,
#     pfDeepFlavourJetTagsDeepAK8,
#     pfDeepDoubleXTagInfosDeepAK8,
#     pfDeepDoubleBvLJetTagsDeepAK8,
#     pfDeepDoubleCvLJetTagsDeepAK8,
#     pfDeepDoubleCvBJetTagsDeepAK8,
#     pfMassIndependentDeepDoubleBvLJetTagsDeepAK8,
#     pfMassIndependentDeepDoubleCvLJetTagsDeepAK8,
#     pfMassIndependentDeepDoubleCvBJetTagsDeepAK8,
#     pfDeepBoostedJetTagInfosDeepAK8,
#     pfDeepBoostedJetTagsDeepAK8,
#     pfDeepBoostedDiscriminatorsJetTagsDeepAK8,
#     pfMassDecorrelatedDeepBoostedJetTagsDeepAK8,
#     pfMassDecorrelatedDeepBoostedDiscriminatorsJetTagsDeepAK8
# )


##########################################################################################
#                                       Jet Smearing                                     #
##########################################################################################
# see https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/PhysicsTools/PatUtils/python/patPFMETCorrections_cff.py#L106
from RecoMET.METProducers.METSigParams_cfi import *
JetSmearing = cms.EDProducer("SmearedPATJetProducer",
    src             = cms.InputTag("slimmedJets"),
    enabled         = cms.bool(True),   # If False, no smearing is performed
    rho             = cms.InputTag("fixedGridRhoFastjetAll"),
    skipGenMatching = cms.bool(False),  # If True, always skip gen jet matching and smear jet with a random gaussian
    algopt          = cms.string('AK4PFchs_pt'),
    algo            = cms.string('AK4PFchs'),
    genJets         = cms.InputTag("slimmedGenJets"),
    dRMax           = cms.double(0.2),  # = cone size (0.4) / 2
    dPtMaxFactor    = cms.double(3),    # dPt < 3 * resolution
    variation       = cms.int32(0),     # If not specified, default to 0
    seed            = cms.uint32(37428479), # If not specified, default to 37428479
    debug           = cms.untracked.bool(False)
)

##############################
# AK4 CHS Jets
mySmearedJets = JetSmearing.clone(
    src       = cms.InputTag("updatedPatJetsUpdatedJEC"),
    variation = cms.int32(0)
)
mySmearedJetsUp = mySmearedJets.clone(
    variation = cms.int32(+1)
)
mySmearedJetsDown = mySmearedJets.clone(
    variation = cms.int32(-1)
)

##############################
# AK4 Puppi Jets
mySmearedPuppiJets = JetSmearing.clone(
    src       = cms.InputTag("updatedPatJetsPuppiJEC"),
    algopt    = cms.string('AK4PFPuppi_pt'),
    algo      = cms.string('AK4PFPuppi'),
    genJets   = cms.InputTag("slimmedGenJets"),
    variation = cms.int32(0)
)
mySmearedPuppiJetsUp = mySmearedPuppiJets.clone(
    variation = cms.int32(+1)
)
mySmearedPuppiJetsDown = mySmearedPuppiJets.clone(
    variation = cms.int32(-1)
)

##############################
# DeepAK8 Jets
mySmearedDeepAK8Jets = JetSmearing.clone(
    src       = cms.InputTag("selectedUpdatedPatJetsDeepAK8"),
    algopt    = cms.string('AK8PFPuppi_pt'),
    algo      = cms.string('AK8PFPuppi'),
    genJets   = cms.InputTag("slimmedGenJetsAK8"),
    dRMax     = cms.double(0.4),  # = cone size (0.8) / 2
    variation = cms.int32(0)
)
mySmearedDeepAK8JetsUp = mySmearedDeepAK8Jets.clone(
    variation = cms.int32(+1)
)
mySmearedDeepAK8JetsDown = mySmearedDeepAK8Jets.clone(
    variation = cms.int32(-1)
)