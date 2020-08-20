import FWCore.ParameterSet.Config as cms

DVAOD = cms.EDAnalyzer('DVAnalyzerAOD',
                      beamspotTag = cms.untracked.InputTag("offlineBeamSpot"),
                      #genJetTag = cms.untracked.InputTag("ak4GenJets"),
                      jetTag = cms.untracked.InputTag("ak4PFJets"),
                      trackTag = cms.untracked.InputTag("generalTracks"),
                      muonTag = cms.untracked.InputTag("muons"),
                      vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                      triggerTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
                      processName = cms.untracked.string("HLT"),
                      triggerNameTag = cms.untracked.string("HLT_PFHT1050_v14"),
                      debug = cms.untracked.bool(False),
                    )
