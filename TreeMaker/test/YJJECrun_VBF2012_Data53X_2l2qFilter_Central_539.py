import FWCore.ParameterSet.Config as cms

process = cms.Process("COMBINE")
runOnMC = False # to run on DATA
#runOnMC = True
#DoJEC=0,1,-1
DoJEC = 0

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tag for data
#######################################################3
if (runOnMC):  ##########used in 53X
        process.GlobalTag.globaltag = 'START53_V27::All'
else:
        process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB2013")
#process.load("RecoBTag.PerformanceDB.BTagPerformanceDB2013")

# needed for TransientTrackBuilder
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')



## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
#  'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/VBF_HToZZTo2L2Q_M-200_8TeV-powheg-pythia6_h2l2qSkimData_1_1_ClF.root'
#   'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/GluGluToHToZZTo2L2Q_M-200_8TeV-powheg-pythia6_SkimPAT_H2l2q_523_v3_l_h2l2qSkimData_1_1_JVm.root'
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_51_1_2MJ.root',
#  'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_205_1_stI.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/ZZ_TuneZ2star_8TeV_pythia6_tauolah2l2qSkimData_382_1_foR.root' 
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_104_1_sZg.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50filter_8TeV-madgraph_h2l2qSkimData_3_1_v1i.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_15_1_Jz1.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_219_1_RZs.root'
#'root://eoscms//eos/cms/store/mc/Summer12/GluGluToHToZZTo2L2Q_M-200_8TeV-powheg-pythia6/AODSIM/PU_S7_START1B1A-3AB8-E111-8D00-003048FFCB8C.root'
#'file:/data4/yunju/VBF2012/skim2l2q/SkimPAT_H2l2q_523_v3_l_GluGluToHToZZTo2L2Q_M-200h2l2qSkimData_22_1_2KD.root'
#'file:/home/yunju/cdxfe_Higgs/Test_525skim/skimsample/YJTest_h2l2qSkimData.root'
#'file:/data4/yunju/VBF2012/SkimPAT_H2l2q_533_v1_GluGluToHToZZTo2L2Q_M-200_h2l2qSkimData_2_1_5SE_2nd.root'
#'file:/data4/yunju/VBF2012/DoubleElectron_Run2012C-24Aug2012-v1h2l2qSkimData_26_1_DhY.root'
#'file:/data4/yunju/VBF2012/DoubleMuSkimPAT_H2l2q_533_v1-Run2012C-PromptReco-v2h2l2qSkimData_101_1_EhS.root'
#'file:/afs/cern.ch/work/y/yunju/private/Samples2l2q533/SkimPAT_H2l2q_533_v1_GluGluToHToZZTo2L2Q_M-200_h2l2qSkimData_2_1_5SE_2nd.root'
#'file:/scratch/yunju/VBF2012/JECUN/CMGTools/CMSSW_5_3_3_patch3/src/h2l2qSkimData.root'   
#'file:/afs/cern.ch/user/m/mmozer/public/h2l2qSkimData.root'
#'file:/afs/cern.ch/work/y/yunju/private/Collection_Sample/ZZ2l2q_539/SkimPAT_H2l2q_539_v4_GluGluToHToZZTo2L2Q_M-300h2l2qSkimData_11_1_kCT.root'
#'file:/afs/cern.ch/work/y/yunju/private/Collection_Sample/ZZ2l2q_539/SkimPAT_H2l2q_539_v4_GluGluToHToZZTo2L2Q_M-600h2l2qSkimData_10_1_pJs.root'
#'file:/afs/cern.ch/work/y/yunju/private/Collection_Sample/ZZ2l2q_539/DiEleSkimPAT_H2l2q_539_v4-Run2012B-22Jan2013-v1h2l2qSkimData_100_1_AEO.root'

'file:/afs/cern.ch/work/y/yunju/private/Collection_Sample/ZZ2l2q_539/DiMuSkimPAT_H2l2q_539_v4-Run2012B-22Jan2013-v1h2l2qSkimData_100_1_Mzq.root'


   )
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1))



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#baseJetSel = cms.PSet(
#  Jets=cms.InputTag("cleanPatJetsNoPUIsoLept")
#   JetsPY=cms.InputTag("customPFJets")
#)

#baseJetSel = cms.PSet(
#  Jets=cms.InputTag("pfJetsAK5")
#)


from DelPanj.TreeMaker.eSelYJVBF_cff import *
from DelPanj.TreeMaker.muSelYJVBF_cff import *


process.tree1 = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(False),
	fillEventInfo_ = cms.bool(True),
	fillGenInfo_   = cms.bool(False),
	fillMuonInfo_  = cms.bool(False),
	fillElecInfo_  = cms.bool(False),
	fillElecIsoInfo_ = cms.bool(False),
	fillJetInfo_   = cms.bool(False),
	fillMetInfo_   = cms.bool(False),
	fillTrigInfo_  = cms.bool(False),
	fillPhotInfo_  = cms.bool(False),
	fillZJetPlant_ = cms.bool(False),
	fillZZInfo_    = cms.bool(False),
        fillYJHiggInfo_= cms.bool(False),       
        fillYJTestSkim_= cms.bool(False), 
        fillJecYJ_     = cms.bool(True),
        eleRhoIso = cms.InputTag("kt6PFJets","rho"),# for rho in eSelector
        muoRhoIso = cms.InputTag("kt6PFJetsCentralNeutral", "rho"),#not using rho in muob any more
        
        e2012IDSet  = eSelYJVBF, # eSelector
        mu2012IDSet = muSelYJVBF,# muSelector
        
        hzzeejjTag = cms.InputTag("hzzeejj:h"),
        hzzmmjjTag = cms.InputTag("hzzmmjj:h"),
	
        genPartLabel=cms.InputTag("genParticles"),
	
        patMuonsPY=cms.InputTag("userDataSelectedMuons"),
	patElectronsPY = cms.InputTag("userDataSelectedElectrons"),
        JetsPY=cms.InputTag("customPFJetsNoPUSub"), 
        

        EleRhoPY =cms.InputTag("kt6PFJets","rho"),#for rho in patEletree
        leadElecPset_ = eSelYJVBF,#in pateleisotree not using now 
	patMet=cms.InputTag("patMETs"),
	beamSpotLabel=cms.InputTag("offlineBeamSpot"),
       	outFileName=cms.string('outputFileName.root'),
        SetJEC_PY=cms.int32(DoJEC)  	
)



process.counter_original = cms.EDAnalyzer('YJEventCounter',
   instance = cms.int32(1) 
)


process.TFileService = cms.Service("TFileService",
      #fileName = cms.string("DYJ200YJTest_v2.root"),
    #  fileName = cms.string("ZZ200YJTest_v2.root"),
      #fileName = cms.string("GGH200YJTest_v2.root"),
      #fileName = cms.string("Mjj_VBF200YJTest_v2.root"),
      fileName = cms.string("JECVBFTree_Data_539_YJAna2l2qfil.root"),




       closeFileFast = cms.untracked.bool(True)
  )

#process.metInfoProducer = cms.EDProducer(
#     "MetVariablesProducer",
#     metTag = cms.InputTag("patMETsAK5"),
#     t1CorrMetTag = cms.InputTag("patType1CorrectedPFMetAK5")
#)

ZZ2l2qFilter = cms.EDFilter('ZZ2l2qFilter')
process.ZZ2l2qfil = ZZ2l2qFilter

process.p = cms.Path(
	process.counter_original*
 #       process.metInfoProducer*
        process.ZZ2l2qfil*
	process.tree1##Trigger Applied.
	)
  
 



