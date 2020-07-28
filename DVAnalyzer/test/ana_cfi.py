import FWCore.ParameterSet.Config as cms

process = cms.Process("TES1")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                              'file:myOutputFile.root'
                            )
                           )

process.ana = cms.EDAnalyzer('DVAnalyzer',
                      jetTag = cms.untracked.InputTag("slimmedJets"),
                      trackRefTag = cms.untracked.InputTag("test"),
                      #tracksTag = cms.untracked.InputTag("packedPFCandidates"),
                      muonTag = cms.untracked.InputTag("slimmedMuons"),
                      vertexTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                      triggerTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
                      processName = cms.untracked.string("HLT"),
                      triggerNameTag = cms.untracked.string("HLT_PFHT1050_v18"),
)


process.TFileService = cms.Service("TFileService", fileName = cms.string("Demo_sig_test.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v15'

DemoAnalyzer = process.ana
process.runseq = cms.Sequence()
process.runseq += DemoAnalyzer
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)
