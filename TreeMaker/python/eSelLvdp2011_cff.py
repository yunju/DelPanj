
import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.eSelBase_cff import *

##The isolation: non-combined, relative, detector based.
eIsoLvdpBrl = eIsoBase.clone()
##eIsoLvdpBrl.iso1 = cms.double(0.09)##tkIso
##eIsoLvdpBrl.iso2 = cms.double(0.07)##ecIso
##eIsoLvdpBrl.iso3 = cms.double(0.10)##hcIso
eIsoLvdpBrl.iso4 = cms.double(0.20)##hcIso



eIsoLvdpEcp = eIsoLvdpBrl.clone()
#eIsoLvdpEcp.iso1 = cms.double(0.04)##tkIso
#eIsoLvdpEcp.iso2 = cms.double(0.05)##ecIso
#eIsoLvdpEcp.iso3 = cms.double(0.025)##hcIso
eIsoLvdpEcp.iso4 = cms.double(0.20)##hcIso



##These are the gsf like id cuts, 
##ANIL: MVA cut not included yet.
eIdLvdpBrl = eIdBase.clone()
eIdLvdpBrl.sieie = cms.double(0.01)
eIdLvdpBrl.delphi= cms.double(0.06)
eIdLvdpBrl.detain = cms.double(0.004)
eIdLvdpBrl.dist = cms.double(0.02)
eIdLvdpBrl.dcot = cms.double(0.02)
eIdLvdpBrl.hoe =  cms.double(0.04)
eIdLvdpBrl.nmisHit=cms.double(0.)

eIdLvdpEcp = eIdLvdpBrl.clone()
eIdLvdpEcp.sieie = cms.double(0.03)
eIdLvdpEcp.delphi= cms.double(0.03)
eIdLvdpEcp.detain = cms.double(0.007)
eIdLvdpEcp.hoe =  cms.double(0.15)


#Make the complete selection PSet.
eSelLvdp2011 = eSelBaseWithNonCombRelDetIso.clone()
#rint eSelLvdp2011
eSelLvdp2011.ptx = cms.double(20)
eSelLvdp2011.etax = cms.double(2.4) 
eSelLvdp2011.idBrl = eIdLvdpBrl
eSelLvdp2011.idEcp  = eIdLvdpEcp#OAp
eSelLvdp2011.isoBrl = eIsoLvdpBrl
eSelLvdp2011.isoEcp = eIsoLvdpEcp

#LvdpSel2011 without Iso, Id
eSelLvdp2011NoIsoNoId = eSelLvdp2011.clone()
eSelLvdp2011NoIsoNoId.ptx = cms.double(20.)
eSelLvdp2011NoIsoNoId.etax = cms.double(2.4)
eSelLvdp2011NoIsoNoId.isoBrl = eIsoBase
eSelLvdp2011NoIsoNoId.isoEcp = eIsoBase
eSelLvdp2011NoIsoNoId.idBrl = eIdBase
eSelLvdp2011NoIsoNoId.idEcp = eIdBase

##LvdpSel2011 without Isolation.
eSelLvdp2011NoIso = eSelLvdp2011.clone()
eSelLvdp2011NoIso.isoBrl = eIsoBase
eSelLvdp2011NoIso.isoEcp = eIsoBase
#print eSelLvdp2011NoIso
##LvdpSel2011 without Identification.
eSelLvdp2011NoId = eSelLvdp2011.clone()
eSelLvdp2011NoId.idBrl = eIdBase
eSelLvdp2011NoId.idEcp = eIdBase


#LvdpSel2011 without Iso, Id
eSelLvdp2011NoPtEtaNoIsoNoId = eSelLvdp2011.clone()
eSelLvdp2011NoPtEtaNoIsoNoId.ptx = 0.
eSelLvdp2011NoPtEtaNoIsoNoId.etax = 200.
eSelLvdp2011NoPtEtaNoIsoNoId.isoBrl = eIsoBase
eSelLvdp2011NoPtEtaNoIsoNoId.isoEcp = eIsoBase
eSelLvdp2011NoPtEtaNoIsoNoId.idBrl = eIdBase
eSelLvdp2011NoPtEtaNoIsoNoId.idEcp = eIdBase


 

