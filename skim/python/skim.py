import FWCore.ParameterSet.Config as cms

##____________________________________________________________________________||
process = cms.Process("SKIM")

##____________________________________________________________________________||
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


##___________________________Global_Tag_______________________________________||
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_P_V56::All'
from RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi import *
process.GlobalTag.globaltag = "100X_dataRun2_Prompt_v1"    #'80X_dataRun2_Prompt_v8'
#80X_dataRun2_Express_v5'


##___________________________Input_Files______________________________________||

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/Run2016B_int/filter_comp/crab_SingleMuon_2016B_RECO_GlobalAndPFcheck/results/crab_pickevents_TL/crab_pickEvents_TL/results/pickevents_1.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/Run2016B_int/filter_comp/crab_SingleMuon_2016B_RECO_GlobalAndPFcheck/results/crab_pickevents_TL/crab_pickEvents_TL/results/pickevents_2.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/Run2016B_int/filter_comp/crab_SingleMuon_2016B_RECO_GlobalAndPFcheck/results/crab_pickevents_TL/crab_pickEvents_TL/results/pickevents_3.root'
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/252/00000/A4345B88-404B-E811-894C-FA163EECC2C3.roor"
	"root://cms-xrd-global.cern.ch///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/264/00000/12AB3330-A44B-E811-92A8-02163E012E66.root"
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/257/00000/1C3A2CCE-604B-E811-99DA-FA163E1FBD5A.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/257/00000/50AE2707-524B-E811-8B45-FA163E39DB61.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/257/00000/9618A647-6F4B-E811-8C23-FA163E355512.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/258/00000/8041ACE6-344B-E811-92B0-FA163E6E32F8.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/259/00000/78C5BDC7-614B-E811-882D-FA163E47F3B0.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/259/00000/8250BBFF-694B-E811-839F-FA163E10001E.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/259/00000/CE9909F1-5D4B-E811-BB1A-FA163E761665.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/259/00000/D46C9006-664B-E811-AA92-FA163EC9355E.root",
	#"root://xrootd-cms.infn.it///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/259/00000/D46C9006-664B-E811-AA92-FA163EC9355E.root"
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/python/crab_pickevents_FRB/crab_pickEvents/results/skim.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/Run2016B_int/filter_comp/pickevents/crab_pickevents_filterbyNew/crab_pickEvents/results/skim_FNew.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/RECO_RAW2DIGI_L1Reco_RECO_ALCA_SKIM_EI_PAT_DQM_inRECO.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/crab_pickevents_overlapHM/crab_pickEvents/results/pickSkim.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/MetScanning/skim/crab/crab_pickevents_filterOverlap/crab_pickEvents/results/pickSkim.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/remainlastRB/RECO_RAW2DIGI_L1Reco_RECO_ALCA_SKIM_EI_PAT_DQM_inRECO.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/lastRB/RECO_RAW2DIGI_L1Reco_RECO_ALCA_SKIM_EI_PAT_DQM_inRECO.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/edgeMuons/RECO_RAW2DIGI_L1Reco_RECO_ALCA_SKIM_EI_PAT_DQM_inRECO.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/rereco/CMSSW_8_0_13_patch1/src/crab_pickevents_EdgeRAW/crab_pickEvents/results/pickevents_EdgeMuonRAW_1.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/edmpick/crab_pickevents_AOD/crab_pickEvents/results/pickevents_1.root' #to check error message for AOD
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/edmpick/crab_pickevents_FEVT/crab_pickEvents/results/pickevents_1.root' #GlobalReport
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/edmpick/crab_pickevents_muonRB/crab_pickEvents/results/pickevents_1.root' #GlobalReport
        #'file:/nfs/dust/cms/user/singha/SUSY_2016C/CMSSW_8_0_11/src/CMGTools/TTHAnalysis/python/plotter/pick_DMD.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/crab/crab_pickevents_22jul-jun_allpass/crab_pickEvents_26jul/results/skim_22Jul-Jun.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/crab/crab_pickevents_22jul-jun_allpass/crab_pickEvents_26jul/results/skim_22Jul-Jun.root'
        #'file:/nfs/dust/cms/group/susy-desy/Run2/Data/2016/ICHEP_12p88fb_2016BCD_HT350Skim_fromDESY/Friends/teste/pickevents.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/zeynep/pickevents_blackholes.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/ece/pickevents_4Muonevent_rereco.root'
       #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/python/ece/ece_single_Muon_rereco.root'
        #'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW808_bugfix/CMSSW_8_0_8/src/MetScanning/skim/crab/crab_pickevents_allPass/crab_pickEvents_all/results/Pickskim_allPass.root' 
#        'file:/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_8/src/MetScanning/skim/crab/raman_events/pickeventsAllFilters.root'
        #'file:/hdfs/store/user/khurana/MET/crab_MET_Run2016B-PromptReco-v2_RECO/160524_134615/0000/skim_406.root'
        #'root://xrootd.unl.edu//store/express/Run2016B/ExpressPhysics/FEVT/Express-v1/000/272/007/00000/84DC1CFA-830D-E611-9BFC-02163E012058.root'
       # 'root://xrootd.unl.edu//store/express/Run2016B/ExpressPhysics/FEVT/Express-v1/000/272/008/00000/00100D61-870D-E611-9B62-02163E01470B.root'
)
    )


##___________________________EDM_Output_File__________________________________||
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('skim.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring(
        'keep *'
#        'keep *_pfClusterMet_*_*', 'keep *_CSCTightHaloFilter_*_*', 'keep *_HBHENoiseFilterResultProducer_*_*'
        )
    )


##____________________________________________________________________________||

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


##___________________________CSC_Halo_Filter__________________________________||
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.CSCTightHalo2015Filter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.CSCTightHaloTrkMuUnvetoFilter_cfi')
process.CSCTightHaloTrkMuUnvetoFilter.taggingMode = cms.bool(True)

##___________________________Global_Halo_Filter__________________________________||
process.load('RecoMET.METFilters.globalTightHalo2016Filter_cfi')
process.globalTightHalo2016Filter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)

#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#    reverseDecision = cms.bool(False)
#)
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('RecoMET.METFilters.HcalStripHaloFilter_cfi')
process.HcalStripHaloFilter.taggingMode = cms.bool(True)


##___________________________PFClusterMet_____________________________________||
process.load('RecoMET.METProducers.PFClusterMET_cfi')
process.pfClusterRefsForJetsHCAL = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHCAL'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsECAL = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterECAL'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsHF = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHF'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsHO = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHO'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJets = cms.EDProducer("PFClusterRefCandidateMerger",
  src = cms.VInputTag("pfClusterRefsForJetsHCAL", "pfClusterRefsForJetsECAL", "pfClusterRefsForJetsHF", "pfClusterRefsForJetsHO")
)
process.pfClusterMetSequence = cms.Sequence(
 process.particleFlowRecHitECAL*
 process.particleFlowRecHitHBHE*
 process.particleFlowRecHitHF*
 process.particleFlowRecHitHO*
 process.particleFlowClusterECALUncorrected*
 process.particleFlowClusterECAL*
 process.particleFlowClusterHBHE*
 process.particleFlowClusterHCAL*
 process.particleFlowClusterHF*
 process.particleFlowClusterHO*
 process.pfClusterRefsForJetsHCAL*
 process.pfClusterRefsForJetsECAL*
 process.pfClusterRefsForJetsHF*
 process.pfClusterRefsForJetsHO*
 process.pfClusterRefsForJets*
#   process.ak4PFClusterJets
 process.pfClusterMet
)


##___________________________PFCaloMet_____________________________________||
process.hltParticleFlowBlock = cms.EDProducer("PFBlockProducer",
  debug = cms.untracked.bool(False),
  verbose = cms.untracked.bool(False),
  elementImporters = cms.VPSet(
      cms.PSet(
          source = cms.InputTag("particleFlowClusterECAL"),
          #source = cms.InputTag("particleFlowClusterECALUncorrected"), #we use uncorrected
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHCAL"),
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHO"),
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHF"),
          importerName = cms.string('GenericClusterImporter')
      )
  ),
  linkDefinitions = cms.VPSet(
      cms.PSet(
          linkType = cms.string('ECAL:HCAL'),
          useKDTree = cms.bool(False),
          #linkerName = cms.string('ECALAndHCALLinker')
          linkerName = cms.string('ECALAndHCALCaloJetLinker') #new ECal and HCal Linker for PFCaloJets
      ),
      cms.PSet(
          linkType = cms.string('HCAL:HO'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('HCALAndHOLinker')
      ),
      cms.PSet(
          linkType = cms.string('HFEM:HFHAD'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('HFEMAndHFHADLinker')
      ),
      cms.PSet(
          linkType = cms.string('ECAL:ECAL'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('ECALAndECALLinker')
      )
   )
)
from RecoParticleFlow.PFProducer.particleFlow_cfi import particleFlowTmp
process.hltParticleFlow = particleFlowTmp.clone(
    GedPhotonValueMap = cms.InputTag(""),
    useEGammaFilters = cms.bool(False),
    useEGammaElectrons = cms.bool(False), 
    useEGammaSupercluster = cms.bool(False),
    rejectTracks_Step45 = cms.bool(False),
    usePFNuclearInteractions = cms.bool(False),  
    blocks = cms.InputTag("hltParticleFlowBlock"), 
    egammaElectrons = cms.InputTag(""),
    useVerticesForNeutral = cms.bool(False),
    PFEGammaCandidates = cms.InputTag(""),
    useProtectionsForJetMET = cms.bool(False),
    usePFConversions = cms.bool(False),
    rejectTracks_Bad = cms.bool(False),
    muons = cms.InputTag(""),
    postMuonCleaning = cms.bool(False),
    usePFSCEleCalib = cms.bool(False)
    )
process.load("RecoMET.METProducers.PFMET_cfi")
process.pfCaloMet = process.pfMet.clone(
  src = cms.InputTag("hltParticleFlow"),
  alias = cms.string('pfCaloMet')
)
process.pfCaloMetSequence = cms.Sequence(
   process.hltParticleFlowBlock *
   process.hltParticleFlow *
   process.pfCaloMet
)

process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.taggingMode = cms.bool(True)
#process.BadPFMuonFilter.debug = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateSummer16Filter_cfi')
process.BadChargedCandidateSummer16Filter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')
process.BadPFMuonSummer16Filter.taggingMode = cms.bool(True)
#process.BadPFMuonFilter.debug = cms.bool(True)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB=cms.vint32(12, 13, 14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE=cms.vint32(12, 13, 14)

process.condMETSelector = cms.EDProducer(
   "CandViewShallowCloneCombiner",
   decay = cms.string("caloMet pfMet"),
   cut = cms.string("(daughter(0).pt > 200) || (daughter(1).pt > 200)" ) 
   )

process.metCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("condMETSelector"),
    minNumber = cms.uint32(1),
    )


##___________________________Flat_Tuple________________________________________||
process.metScanNtupleMaker = cms.EDAnalyzer("METScanningNtupleMaker",
                                            isReco = cms.bool(True),
                                            rootOutputFile=cms.string("tuple.root"),
                                            muonCandidates = cms.InputTag("muons"),
                                            pfCandidates=cms.InputTag("particleFlow"),
                                            pfJets=cms.InputTag("ak4PFJets"),
                                            caloMET=cms.InputTag("caloMet"),
                                            pfCaloMET=cms.InputTag("pfCaloMet"),
                                            pfClusterMET=cms.InputTag("pfClusterMet"),
                                            pfMET=cms.InputTag("pfMet"),
                                            EcalPFClusterCollection=cms.InputTag("particleFlowClusterECAL"),
                                            HcalPFClusterCollection=cms.InputTag("particleFlowClusterHCAL"),
                                            HBHEPFClusterCollection=cms.InputTag("particleFlowClusterHBHE"),
                                            HOPFClusterCollection=cms.InputTag("particleFlowClusterHO"),
                                            HFPFClusterCollection=cms.InputTag("particleFlowClusterHF"),
                                            tracksCollection=cms.InputTag("generalTracks"),
                                            TRKfilterLETMC=cms.InputTag("logErrorTooManyClusters"),
                                            TRKfilterLETMS=cms.InputTag("logErrorTooManySeeds"),
                                            TRKfilterMSC=cms.InputTag("manystripclus53X"),
                                            TRKfilterTMSC=cms.InputTag("toomanystripclus53X"),
                                            CSC2015filter=cms.InputTag("CSCTightHalo2015Filter"),
                                            GlobalHalofilterTight=cms.InputTag("globalTightHalo2016Filter"),
                                            GlobalHalofilterSuperTight=cms.InputTag("globalSuperTightHalo2016Filter"),
                                            HcalStripHaloFilter=cms.InputTag("HcalStripHaloFilter"),
                                            HBHEfilterR1=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1"),
                                            HBHEfilterR2L=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
                                            HBHEfilterR2T=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
                                            HBHEfilterISO=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1"),
                                            ECALTPfilter=cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"),
                                            ECALSCfilter=cms.InputTag("eeBadScFilter"),
                                            EBRecHits=cms.InputTag("reducedEcalRecHitsEB"),
                                            EERecHits=cms.InputTag("reducedEcalRecHitsEE"),
                                            ESRecHits=cms.InputTag("reducedEcalRecHitsES"),
                                            HcalNoise=cms.InputTag("hcalnoise"),
                                            BadChCandFilter=cms.InputTag("BadChargedCandidateFilter"),
                                            EcalBadCalibSummer17Filter = cms.InputTag("EcalBadCalibSummer17Filter"),
					    BadPFMuon=cms.InputTag("BadPFMuonFilter"),
                                            BadChCandSummer16Filter=cms.InputTag("BadChargedCandidateSummer16Filter"),
                                            BadPFMuonOld=cms.InputTag("BadPFMuonSummer16Filter"),
                                            OfflinePrimaryVertices = cms.InputTag("offlinePrimaryVertices")
                                            

)

# This part is needed if you want to update the BeamHaloSummary information
from RecoMET.METProducers.CSCHaloData_cfi import *
from RecoMET.METProducers.EcalHaloData_cfi import *
from RecoMET.METProducers.HcalHaloData_cfi import *
from RecoMET.METProducers.GlobalHaloData_cfi import *
from RecoMET.METProducers.BeamHaloSummary_cfi import *

process.BeamHaloId = cms.Sequence(CSCHaloData*EcalHaloData*HcalHaloData*GlobalHaloData*BeamHaloSummary)




##___________________________PATH______________________________________________||
process.p = cms.Path(
#    process.BeamHaloId* #Uncomment this if you want to rerun the BeamHaloSummary. By default this line should remain commented
    process.primaryVertexFilter*
    process.bunchSpacingProducer *
    process.condMETSelector *
    process.metCounter* #uncomment this line to apply a met cut
    process.CSCTightHaloFilter*
    process.HBHENoiseFilterResultProducer* #produces bools    
#    process.ApplyBaselineHBHENoiseFilter* 
    process.EcalDeadCellTriggerPrimitiveFilter*
    process.pfClusterMetSequence*
    process.pfCaloMetSequence*
    process.eeBadScFilter*
    process.goodVertices*
    process.trackingFailureFilter*
    process.EcalDeadCellBoundaryEnergyFilter*
    process.CSCTightHalo2015Filter*
    process.CSCTightHaloTrkMuUnvetoFilter*
    process.globalTightHalo2016Filter * 
    process.globalSuperTightHalo2016Filter * 
    process.HcalStripHaloFilter*
    process.BadChargedCandidateFilter*
    process.BadPFMuonFilter*
    process.BadChargedCandidateSummer16Filter*
    process.BadPFMuonSummer16Filter*
    process.metScanNtupleMaker ##CH: writes a flat tree
    )

#process.e1 = cms.EndPath(
 #   process.out ##CH: write the skimmed edm file 
  #  )

##____________________________________________________________________________||
