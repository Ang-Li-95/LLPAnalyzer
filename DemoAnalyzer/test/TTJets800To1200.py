from CRABClient.UserUtilities import config, getUsernameFromCRIC
config = config()

config.General.requestName = 'TTJets800To1200'
config.General.workArea = 'crab_0726'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'DVanaAOD_cfg.py'
#config.JobType.maxMemoryMB = 8000
#config.JobType.numCores = 4
#config.JobType.maxJobRuntimeMin = 2700
#make crab work for CMSSW_9_X on slc7 nodes
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'
config.Data.inputDataset = '/TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
#config.Data.ignoreLocality = True
#config.Data.unitsPerJob = 1
#NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.totalUnits = 100
config.Data.outLFNDirBase = '/store/user/%s' % (getUsernameFromCRIC())
config.Data.publication = False
config.Data.outputDatasetTag = 'TTJets800To1200_AOD'

#config.Site.whitelist = ['T2_US_Purdue','T1_US_FNAL','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'
