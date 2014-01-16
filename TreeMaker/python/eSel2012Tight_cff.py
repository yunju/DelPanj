
import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.eSelBase_cff import *

##The isolation: non-combined, relative, detector based.

#Barrel
eIdBase2012Brl = eIdBase.clone()
eIdBase2012Brl.detain = cms.double(0.004)
eIdBase2012Brl.delphi= cms.double(0.030)
eIdBase2012Brl.sieie = cms.double(0.010)
eIdBase2012Brl.hoe =  cms.double(0.120)
eIdBase2012Brl.d0vtx = cms.double(0.020)
eIdBase2012Brl.dzvtx = cms.double(0.100)
eIdBase2012Brl.ooemoop = cms.double(0.050)
eIdBase2012Brl.dist = cms.double(0.020)
eIdBase2012Brl.dcot = cms.double(0.020)
eIdBase2012Brl.hasConv = cms.double(1e-6)
eIdBase2012Brl.nmisHit=cms.double(0.)


eIsoBase2012Brl = eIsoBase.clone()
eIsoBase2012Brl.iso1 = cms.double(999.99)     
eIsoBase2012Brl.iso2 = cms.double(999.99)     
eIsoBase2012Brl.iso3 = cms.double(999.99)      
eIsoBase2012Brl.iso4 = cms.double(0.10)##hcIso

#Endcap
eIdBase2012Ecp = eIdBase.clone()
eIdBase2012Ecp.detain = cms.double(0.005)
eIdBase2012Ecp.delphi= cms.double(0.020)
eIdBase2012Ecp.sieie = cms.double(0.030)
eIdBase2012Ecp.hoe =  cms.double(0.100)
eIdBase2012Ecp.d0vtx = cms.double(0.020)
eIdBase2012Ecp.dzvtx = cms.double(0.100)
eIdBase2012Ecp.ooemoop = cms.double(0.050)
eIdBase2012Ecp.dist = cms.double(0.020)
eIdBase2012Ecp.dcot = cms.double(0.020)
eIdBase2012Ecp.hasConv = cms.double(1e-6)
eIdBase2012Ecp.nmisHit=cms.double(0.)


eIsoBase2012Ecp = eIsoBase.clone()
eIsoBase2012Ecp.iso1 = cms.double(999.99)     
eIsoBase2012Ecp.iso2 = cms.double(999.99)     
eIsoBase2012Ecp.iso3 = cms.double(999.99)      
eIsoBase2012Ecp.iso4 = cms.double(0.07)##hcIso


#Total 2012 Tight ID
eSel2012Tight  = eSelBase.clone()
eSel2012Tight.ptx = 20.
eSel2012Tight.etax = 2.4 
eSel2012Tight.idBrl = eIdBase2012Brl
eSel2012Tight.idEcp = eIdBase2012Ecp
eSel2012Tight.isoBrl = eIsoBase2012Brl
eSel2012Tight.isoEcp = eIsoBase2012Ecp
eSel2012Tight.usePfIso  = cms.bool(True)
eSel2012Tight.useRelIso = cms.bool(True)
eSel2012Tight.useCombIso= cms.bool(True)
eSel2012Tight.coneRad   = cms.double(0.3)#this parameter is not used anymore.
eSel2012Tight.rhoCorr   = cms.bool(True)


