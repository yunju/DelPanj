import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.muSelBase_cff import *

###########################################################
#identification stuff                                     # 
###########################################################
muIdBase2012 = muIdBase.clone()
muIdBase2012.normalizedChi2= cms.double(10.)
muIdBase2012.muonHits=cms.double(0.)
muIdBase2012.nMatches=cms.double(1)
muIdBase2012.dxy = cms.double(0.2)
muIdBase2012.dz = cms.double(0.5)
muIdBase2012.pixelHits=cms.double(0.)
muIdBase2012.trackerHits=cms.double(5.)



muSel2012NoIso = muSelBase.clone()
muSel2012NoIso.pt = cms.double(20)
muSel2012NoIso.eta= cms.double(2.4)
muSel2012NoIso.requireTrigMatch= cms.bool(False)
muSel2012NoIso.idPar  = muIdBase2012
muSel2012NoIso.isoPar = muIsoBase
muSel2012NoIso.usePfIso  = cms.bool(True)
muSel2012NoIso.useRelIso = cms.bool(True)
muSel2012NoIso.useCombIso= cms.bool(True)
muSel2012NoIso.coneRad   = cms.double(0.4)#this parameter is not used anymore.
muSel2012NoIso.delBCorr  = cms.bool(True) #True for deltaBeta correction
                                        #False for rho correction




