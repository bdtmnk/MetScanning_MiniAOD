from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#submit with 'python crab.py'
#Don't write to my directory (schoef), though

config = config()

config.General.requestName = "SingleMuon_Run2018D-PromptReco-v2_MINIAOD"
config.General.workArea = "SingleMuon_Run2018D-PromptReco-v2/MINIAOD"

config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/skimMINIAOD.py'
config.JobType.outputFiles = ['tuple.root']
config.Data.inputDataset = "/SingleMuon/Run2018D-PromptReco-v2/MINIAOD"


#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 18
#config.Data.runRange = '272818'
config.Data.publication = False


#config.Data.ignoreLocality = True                                                                                                                                                                     
#config.Site.whitelist = ['T1_US_FNAL_Disk'] 
#config.Data.outLFNDirBase = '' 
#config.Data.publishDataName = ''

config.Data.outLFNDirBase = '/store/user/%s' % (getUsernameFromSiteDB())
config.Site.storageSite = "T2_DE_DESY"
