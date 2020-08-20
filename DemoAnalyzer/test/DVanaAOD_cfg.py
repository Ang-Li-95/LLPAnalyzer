import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                  '/store/mc/RunIIFall17DRPremix/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/70001/FE846B27-DA84-E811-A32E-FA163EC6168E.root'
            #'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18DRPremix/GluinoGluinoToNeutralinoNeutralinoTo2T2B2S_M-1600_CTau-1mm_TuneCP2_13TeV_2018-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/210000/C8385554-2D21-EE4E-8AF5-9E829F45B6E8.root'
            #'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18DRPremix/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/80000/A68DAF33-F18B-094A-906E-91B56567D3F3.root'
                                  #'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18DRPremix/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/120000/F5EECD04-62D2-1545-AD98-5197674877FD.root'
                )
                            )

process.load('LLPAnalyzer.DemoAnalyzer.DVanalyzerAOD_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v17'
#process.demo = cms.EDAnalyzer('DemoAnalyzer'
#                              )

#process.p = cms.Path(process.demo)
DVAOD = process.DVAOD

process.TFileService = cms.Service("TFileService", fileName = cms.string("DVAOD.root") )

process.runseq = cms.Sequence()
process.runseq += DVAOD
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)
