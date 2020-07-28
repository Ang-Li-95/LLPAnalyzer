import FWCore.ParameterSet.Config as cms

Demo = cms.EDAnalyzer('DemoAnalyzer',
                      #genJetTag = cms.untracked.InputTag("ak4GenJets"),
                      pfJetTag = cms.untracked.InputTag("ak4PFJets"),
                      tracksTag = cms.untracked.InputTag("generalTracks"),
                      muonTag = cms.untracked.InputTag("muons"),
                      vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                      triggerTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
                      processName = cms.untracked.string("HLT"),
                      triggerNameTag = cms.untracked.string("HLT_PFHT1050_v18"),
                    )
