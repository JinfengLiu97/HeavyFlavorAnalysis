from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'Data_Ntuple_XXXXXXX'

config.section_('JobType')
config.JobType.psetName = './BPH_SKIM_Ntuple_data.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['mymultilep.root']

config.section_('Data')
config.Data.inputDataset = '/Charmonium/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD' #2016B
#config.Data.inputDataset = '/Charmonium/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD' #2016C
#config.Data.inputDataset = '/Charmonium/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD' #2016D
#config.Data.inputDataset = '/Charmonium/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD' #2016E
#config.Data.inputDataset = '/Charmonium/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD' #2016FH
#config.Data.inputDataset = '/Charmonium/Run2016F-21Feb2020_UL2016-v1/AOD' #2016F
#config.Data.inputDataset = '/Charmonium/Run2016G-21Feb2020_UL2016-v1/AOD' #2016G
#config.Data.inputDataset = '/Charmonium/Run2016H-21Feb2020_UL2016-v1/AOD' #2016H
#config.Data.inputDataset = '/Charmonium/Run2017C-09Aug2019_UL2017-v1/AOD' #2017C
#config.Data.inputDataset = '/Charmonium/Run2017D-09Aug2019_UL2017-v1/AOD' #2017D
#config.Data.inputDataset = '/Charmonium/Run2017E-09Aug2019_UL2017-v1/AOD' #2017E
#config.Data.inputDataset = '/Charmonium/Run2017F-09Aug2019_UL2017-v1/AOD' #2017F
#config.Data.inputDataset = '/Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/AOD' #2018A
#config.Data.inputDataset = '/Charmonium/Run2018B-12Nov2019_UL2018-v1/AOD' #2018B
#config.Data.inputDataset = '/Charmonium/Run2018C-12Nov2019_UL2018 rsb v2-v2/AOD' #2018C
#config.Data.inputDataset = '/Charmonium/Run2018D-12Nov2019_UL2018-v1/AOD' #2018D

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30

#NJOBS = 5000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/XXXXXX/' #lxplus username here
config.Data.lumiMask = 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt' #2016
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_UL2017_Collisions17_JSON_MuonJSON.txt' #2017
#config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON MuonPhys.txt' #2018

config.section_('User')
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX' 
