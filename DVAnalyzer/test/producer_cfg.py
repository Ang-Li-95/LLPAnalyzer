import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/GluinoGluinoToNeutralinoNeutralinoTo2T2B2S_M-800_CTau-1mm_TuneCP2_13TeV_2018-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/E7677EE5-376B-6F4B-B504-BC820873679A.root'
    )
)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)

#from ProdTutorial.TrackAndPointsProducer.trackandpointsproducer_cfi import *
process.prod = cms.EDProducer('trackRefProd',
                              packed_candidates_src = cms.InputTag('packedPFCandidates'),
                              pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              skip_weirdos = cms.bool(True)
)

process.rescale = cms.EDProducer('RescaleTracks',
                                 trackTag = cms.InputTag('prod'),
                                 doRescale = cms.bool(True)
)

#process.ana = cms.EDAnalyzer('testAnlzr',
#                             tracks = cms.untracked.InputTag('test'),
#                             cands = cms.untracked.InputTag('packedPFCandidates'),
#                             tracks_ref = cms.untracked.InputTag('test')
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
    #,outputCommands = cms.untracked.vstring('drop *',
    #  "keep *_generalTracks_*_*",
    #  "keep *_globalMuons_*_*",
    #  "keep *_MuonTrackPoints_*_*",
    #  "keep *_TrackTrackPoints_*_*")

)

#process.TFileService = cms.Service("TFileService", fileName = cms.string("Demo_sig_test.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v19'

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

process.p = cms.Path(process.jecSequence*process.prod*process.rescale)

process.e = cms.EndPath(process.out)
