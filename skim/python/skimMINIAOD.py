import FWCore.ParameterSet.Config as cms


##____________________________________________________________________________||
process = cms.Process("SKIM")


##____________________________________________________________________________||
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")



##___________________________Global_Tag_______________________________________||
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")
process.GlobalTag.globaltag = "100X_dataRun2_Prompt_v1"


##___________________________Input_Files______________________________________||
HEADER = "root://cms-xrd-global.cern.ch///"
HEADER = "root://srm-cms.cern.ch///" 
HEADER = "root://gridftp.echo.stfc.ac.uk///"

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
	"/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/316/199/00000/AC548C74-6658-E811-A139-02163E019FD9.root"
))

##___________________________EDM_Output_File__________________________________||
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('skim.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring(
        'keep *'
#	"root://cms-xrd-global.cern.ch///store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/264/00000/12AB3330-A44B-E811-92A8-02163E012E66.root"
#        'keep *_pfClusterMet_*_*', 'keep *_CSCTightHaloFilter_*_*', 'keep *_HBHENoiseFilterResultProducer_*_*'
        )
    )


##____________________________________________________________________________||
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )


process.load("RecoMET.METFilters.ecalBadCalibFilter_cfi")
process.ecalBadCalibFilter.taggingMode = cms.bool(True)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )



process.EcalDeadCellBoundaryEnergyFilter = cms.EDFilter('EcalDeadCellBoundaryEnergyFilter',
	recHitsEB = cms.InputTag("reducedEgamma:reducedEBRecHits"),
	recHitsEE = cms.InputTag("reducedEgamma:reducedEERecHits"),
	FilterAlgo= cms.untracked.string("FilterMode"),
	#### the following parameters skimGap, skimDead are only used in TuningMode
	#### switch bool to True to turn on filter: only Events with chosen signature pass, otherwise all events pass
	skimGap = cms.untracked.bool(False),
	skimDead  = cms.untracked.bool(False),
	#### cuts for finding energy deposit near Gaps
	## min. boundary energy (RecHit next to Gap) (abs value)
	cutBoundEnergyGapEE=cms.untracked.double(100),
	cutBoundEnergyGapEB=cms.untracked.double(100),
	#### cuts for finding energy deposit near dead region
	## min. boundary energy (RecHits next to Dead Region) (abs value)
	cutBoundEnergyDeadCellsEB=cms.untracked.double(10),
	cutBoundEnergyDeadCellsEE=cms.untracked.double(10),
	#### Limit complete filter processing to EE or EB, if both are 'True' nothing will happen in the filter at all...
	limitFilterToEB=cms.untracked.bool(False),
	limitFilterToEE=cms.untracked.bool(False),
	#### Limit dead cells to channel status, only rec hits around channel with channel status given are
	#### considered. E.g to sum only energy around dead cells with stati 12 & 14 in EB, but all dead cells
	#### in EB, do:
	#### limitDeadCellToChannelStatusEB=cms.vint32(12,14)
	#### limitDeadCellToChannelStatusEE=cms.vint32()
	#### for negative values all status>=abs(given value) are used (e.g. limitDeadCellToChannelStatusEE=cms.vint32(-13)--->limitDeadCellToChannelStatusEE=cms.vint32(13,14,15,16,17,...))
	limitDeadCellToChannelStatusEB=cms.vint32(12, 13, 14),
	limitDeadCellToChannelStatusEE=cms.vint32(12, 13, 14),
	#### enable calculation of energy deposits next to cracks/gaps
	enableGap=cms.untracked.bool(False),
        taggingMode   = cms.bool(False),
        debug = cms.bool(False),
)






process.condMETSelector = cms.EDProducer(
   "CandViewShallowCloneCombiner",
   isReco = cms.bool(False),
   decay = cms.string("slimmedMETs slimmedMETsPuppi"),
   cut = cms.string("(daughter(0).pt > 200) || (daughter(1).pt > 200)" ) 
   )

process.metCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("condMETSelector"),
    minNumber = cms.uint32(1),
    )
    




process.metScanNtupleMaker = cms.EDAnalyzer("METScanningNtupleMaker",
                                            isReco = cms.bool(False),
                                            RecPrimVertex = cms.untracked.bool(True),
                                            RecPFMet = cms.untracked.bool(True),
                                            RecPuppiMet = cms.untracked.bool(True),
                                            IsData = cms.untracked.bool(True),
                                           rootOutputFile=cms.string("tuple.root"),
                                           #collections:
                                           MuonCollectionTag = cms.InputTag("slimmedMuons"),
                                           MetCollectionTag = cms.InputTag("slimmedMETs"),
                                           PuppiMetCollectionTag = cms.InputTag("slimmedMETsPuppi"),
                                           PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
					   TriggerProcess = cms.untracked.string("HLT"),
					   HLTriggerPaths = cms.untracked.vstring(
								#SingleMuon
							   'HLT_IsoMu20_v',
							   'HLT_IsoMu24_v',
							   'HLT_IsoMu24_eta2p1_v',
							   'HLT_IsoMu27_v',
							   'HLT_IsoMu30_v',
							   'HLT_Mu50_v'),
	                                   FlagsProcesses = cms.untracked.vstring("RECO","PAT"), 
					   Flags = cms.untracked.vstring(
  							'Flag_HBHENoiseFilter',
							'Flag_HBHENoiseIsoFilter',
							'Flag_CSCTightHalo2015Filter',
						        'Flag_EcalDeadCellTriggerPrimitiveFilter',
							'Flag_goodVertices',
					  		'Flag_eeBadScFilter',
							'Flag_chargedHadronTrackResolutionFilter',
							'Flag_muonBadTrackFilter',
						 	'Flag_globalTightHalo2016Filter',
							'Flag_METFilters',
							'Flag_BadChargedCandidateFilter',
							'Flag_BadPFMuonFilter',
							'Flag_ecalBadCalibFilter',
							#'allMetFilterPaths'
					   	), 

#					    BadChargedCandidateFilter =  cms.InputTag("BadChargedCandidateFilter"),
#					    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
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
                                            GlobalHalofilterTight=cms.InputTag("globalTightHalo2016Filter"),
                                            GlobalHalofilterSuperTight=cms.InputTag("globalSuperTightHalo2016Filter"),
                                            HcalStripHaloFilter=cms.InputTag("HcalStripHaloFilter"),
                                            HBHEfilterR1=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1"),
                                            HBHEfilterR2L=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
                                            HBHEfilterR2T=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
                                            HBHEfilterISO=cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1"),
                                            ECALTPfilter=cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"),
                                            ECALSCfilter=cms.InputTag("eeBadScFilter"),
                                            EBRecHits=cms.InputTag("educedEgamma:reducedEERecHits"),
                                            EERecHits=cms.InputTag("reducedEcalRecHitsEE"),
                                            ESRecHits=cms.InputTag("reducedEcalRecHitsES"),
                                            BadChCandFilter=cms.InputTag("BadChargedCandidateFilter"),
                                            EcalBadCalibFilter=cms.InputTag("ecalBadCalibFilter"),
					    OfflinePrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            HcalNoise=cms.InputTag("hcalnoise")
                                            
)



##___________________________PATH______________________________________________||
process.p = cms.Path(
    process.metScanNtupleMaker ##CH: writes a flat tree
    )

#process.e1 = cms.EndPath(
#    process.out ##CH: write the skimmed edm file 
#    )
