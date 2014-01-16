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


##########################################################################
##usePfIso=True) iso1=chHadIso;iso2=neHadIso;iso3=gamIso;iso4=combIso    # 
##else iso1=tkIso;iso2=ecalIso;iso3=hcalIso;iso4=combIso                 #
##########################################################################
muIsoBase2012 = muIsoBase.clone()
muIsoBase2012.iso1 = cms.double(999.99)##tkIso
muIsoBase2012.iso2 = cms.double(999.99)##ecIso
muIsoBase2012.iso3 = cms.double(999.99)##hcIso
muIsoBase2012.iso4 = cms.double(0.12)


muSel2012HZZ = muSelBase.clone()
muSel2012HZZ.pt = cms.double(20)
muSel2012HZZ.eta= cms.double(2.4)
muSel2012HZZ.requireTrigMatch= cms.bool(False)
muSel2012HZZ.idPar  = muIdBase2012
muSel2012HZZ.isoPar = muIsoBase2012
muSel2012HZZ.usePfIso  = cms.bool(True)
muSel2012HZZ.useRelIso = cms.bool(False)
muSel2012HZZ.useCombIso= cms.bool(False)
muSel2012HZZ.coneRad   = cms.double(0.4)#this parameter is not used anymore.
#muSel2012HZZ.delBCorr  = cms.bool(True) #True for deltaBeta correction
muSel2012HZZ.delBCorr  = cms.bool(False) #True for deltaBeta correction
                                        #False for rho correction




