import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.muSelBase_cff import *

##The islation: non-combined, relative, detector based.
muIsoLvdp = muIsoBase.clone()
muIsoLvdp.iso1 = cms.double(999.09)##tkIso
muIsoLvdp.iso2 = cms.double(999.07)##ecIso
muIsoLvdp.iso3 = cms.double(999.10)##hcIso
muIsoLvdp.iso4 = cms.double(0.2)


muIdLvdp = muIdBase.clone()
muIdLvdp.dxy = cms.double(0.2)
muIdLvdp.normalizedChi2= cms.double(10.)
muIdLvdp.trackerHits=cms.double(8.)
muIdLvdp.pixelHits=cms.double(0.)
muIdLvdp.muonHits=cms.double(0.)
muIdLvdp.nMatches=cms.double(1)

muSelLvdp2011 = muSelBase.clone()
muSelLvdp2011.pt = cms.double(20)
muSelLvdp2011.coneRad   = cms.double(0.3)#this parameter is not used anymore.
muSelLvdp2011.usePfIso  = cms.bool(True)
muSelLvdp2011.useRelIso = cms.bool(True)
muSelLvdp2011.isoPar = muIsoLvdp
muSelLvdp2011.idPar  = muIdLvdp


muSelLvdp2012 = muSelLvdp2011.clone()#i dont see no difference between common and above wp.


