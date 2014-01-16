import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")
runOnMC = False # to run on DATA
#runOnMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
       $inputFileNames
#'file:/hdfs/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/002612AC-463D-E111-802A-E0CB4E19F9BC.root'
#'file:/hdfs/store/data/Run2011A/DoubleElectron/AOD/08Nov2011-v1/0001/00288607-B51B-E111-A61D-003048FFCB9E.root'
)
)
if(runOnMC == False):
	import PhysicsTools.PythonAnalysis.LumiList as LumiList
	import FWCore.ParameterSet.Types as CfgTypes
	myLumis = LumiList.LumiList(filename = '/afs/hep.wisc.edu/home/anil79/exercise06/CMSSW_4_4_2/src/PhysicsTools/PatAlgos/test/excbadCert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt').getCMSSWString().split(',')
	process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
	process.source.lumisToProcess.extend(myLumis)

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1))

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
if(runOnMC == False):
	process.GlobalTag.globaltag = cms.string('GR_R_44_V13::All')# autoCond[ 'startup' ] )
else :
        process.GlobalTag.globaltag = cms.string('START44_V12::All')# autoCond[ 'startup' ] )
print process.GlobalTag.globaltag
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
                               )

process.outpath = cms.EndPath(process.out)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.out.fileName = cms.untracked.string('$outputFileName')
#process.out.fileName = cms.untracked.string('/scratch/anil79/patTuple_PF2PAT_442.root')


# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. It is currently 
# not possible to run PF2PAT+PAT and standart PAT at the same time
from PhysicsTools.PatAlgos.tools.pfTools import *
#runOnMC = False

postfix = "PFlow"
jetAlgo="AK5"
if(runOnMC==False):
	usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=runOnMC, postfix = postfix, jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']) )
else:
	usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix,jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']))

process.patElectronsPFlow.pfElectronSource = "pfElectronsPFlow"

from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT

# the following is advocated by JetMET, but leads to include very far tracks in the no pile up collection
enablePileUpCorrectionInPF2PAT( process, postfix=postfix)


if runOnMC == False:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatchingPF2PAT( process, '' )


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
#                                debugOn = cms.untracked.bool(True),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

########################################
# HBHENoiseFilter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

from RecoJets.JetProducers.kt4PFJets_cfi import *
# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True, voronoiRfact = 0.9 )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

if(runOnMC==False):
	process.patseq = cms.Sequence(
		process.noscraping*
		process.primaryVertexFilter*
		process.HBHENoiseFilter *
		getattr(process,"patPF2PATSequence"+postfix)*
                process.kt6PFJetsForIsolation 
		)
	
else:
	process.patseq = cms.Sequence(    
		process.noscraping*
		process.primaryVertexFilter*
		process.HBHENoiseFilter *
		getattr(process,"patPF2PATSequence"+postfix)*
                process.kt6PFJetsForIsolation 
		)


## let it run
process.p = cms.Path(
    process.patseq
    )

# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning

#process.out.outputCommands =  cms.untracked.vstring('drop *')
process.out.outputCommands = cms.untracked.vstring('drop *')
#                                                   'keep recoPFCandidates_particleFlow_*_*',
 #                                                  *patEventContentNoCleaning )

process.out.outputCommands.extend(cms.untracked.vstring(
#    'keep *_goodOfflinePrimaryVertices*_*_*',    
    'keep double_*PFlow*_*_PAT',
    'keep *_offlinePrimaryVertices_*_*',
   # 'keep *_met*_*_*',
    'keep *_selectedPatElectronsPFlow_*_*',
    'keep *_selectedPatJetsPFlow_*_*',
#    'keep *_gsfElectron*_*_*',
    #'keep *_scalersRawToDigi_*_*',
    'keep *_addPileupInfo_*_*',
    'keep *_*kt6PFJets*_*_PAT',
        # 'keep *_genParticles_*_*',
#    'keep *_*generalTracks*_*_*',
 #        'keep *_generator_*_*',
    #     'keep *_particleFlow_*_*',
         # Trigger
         'keep *_TriggerResults_*_*',
         'keep *_hltTriggerSummaryAOD_*_*',
     #    'keep *_pfElectronTranslator_*_*',
 #        'keep *_generalTracks_*_*',
       #  'keep *_ak5PFJets*_*_*',
         'keep *_offlineBeamSpot_*_*',
  #       'keep *_dimu_*_*',
   #      'keep *_zMuMuCands*_*_*',
    #     'keep *_zMuMuCandsMuEta*_*_*',
#         'keep *_hcalnoise_*_*',
 #        'keep *_gtDigis_*_*',
  #       'keep *_muon*METValueMapProducer_*_*',
         #'keep *_offlinePrimaryVertices*_*_*',
    #     'keep *_*pfNoPileUp*_*_*',
         #'keep *_*scalersRawToDigi*_*_*'
                'keep *_superClusters_*_*', 
                'keep *_patElecPtEta*_*_*',
                'keep *_patJetPtEta*_*_*',
                'keep *_genElectrons*_*_*',
                'keep *_ak5CaloJets*_*_*',
                'keep *_*Rcd*_*_*',
                'keep *_*hltL1GtObjectMap*_*_*',
                'keep edmTriggerResults_TriggerResults*_*_*',
                'keep *_hltTriggerSummaryAOD_*_*',
                'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
))

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

# enable delta beta correction for muon selection in PF2PAT? 
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
