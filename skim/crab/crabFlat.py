from CRABClient.UserUtilities import config
#submit with 'python crab.py'
#Don't write to my directory (schoef), though

## define only these variables here

production = "/afs/cern.ch/work/c/cheidegg/crab3/2015-09-29_filters_runc_sarahskim_flat"
json       = "json/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON.txt" # RunC 4T
#json       = "json/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt" # RunD 4T
site       = "T3_CH_PSI"
outdir     = "/store/user/cheidegg/crab3/2015-09-29_filters_runc_sarahskim_flat"


## do not touch beyond this point

config = config()
config.General.requestName   = 'ZeroBias1_Run2015A-PromptReco-v1_RECO'
config.General.workArea      = production

config.JobType.outputFiles   = ['tuple.root']
config.JobType.pluginName    = 'Analysis'
config.JobType.psetName      = '../python/skimFlat.py'

config.Data.inputDataset     = '/ZeroBias1/Run2015A-PromptReco-v1/RECO'
config.Data.inputDBS         = 'global'
config.Data.lumiMask         = json
config.Data.splitting        = 'LumiBased'
config.Data.unitsPerJob      = 20

config.Data.publication      = False
#config.Data.outLFNDirBase   = '' 
#config.Data.publishDataName = ''

config.Data.outLFNDirBase    = outdir
config.Site.storageSite      = site

datasets=[
##'/BTagCSV/Run2015D-PromptReco-v1/RECO',
##'/BTagMu/Run2015D-PromptReco-v1/RECO',
##'/Charmonium/Run2015D-PromptReco-v1/RECO',
##'/DoubleEG/Run2015D-PromptReco-v1/RECO',
##'/DoubleMuon/Run2015D-PromptReco-v1/RECO',
##'/EGamma/Run2015D-PromptReco-v1/RECO',
#'/ExpressPhysics/Run2015D-Express-v3/FEVT',
##'/Jet/Run2015D-PromptReco-v1/RECO',
#'/JetHT/Run2015D-PromptReco-v3/RECO',
#'/HighMultiplicity/Run2015D-PromptReco-v1/RECO',
'/HTMHT/Run2015D-PromptReco-v3/RECO',
#'/MET/Run2015D-PromptReco-v3/RECO',
##'/MinimumBias/Run2015D-PromptReco-v1/RECO',
##'/MuonEG/Run2015D-PromptReco-v1/RECO',
#'/SingleElectron/Run2015D-PromptReco-v3/RECO',
#'/SingleMuon/Run2015D-PromptReco-v3/RECO',
#'/SinglePhoton/Run2015D-PromptReco-v3/RECO',
##'/StreamExpress/Tier0_Test_SUPERBUNNIES_vocms015-Hotline-Express-v28/ALCARECO', 
##'/Tau/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias1/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias2/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias3/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias4/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias5/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias6/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias7/Run2015D-PromptReco-v1/RECO',
##'/ZeroBias8/Run2015D-PromptReco-v1/RECO',
]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.rstrip('/').lstrip('/').replace('/','_') + "_flat"
#        print config.General.requestName
        crabCommand('submit', config = config)



