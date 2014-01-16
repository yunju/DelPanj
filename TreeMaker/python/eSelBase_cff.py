import FWCore.ParameterSet.Config as cms


###########################################################
#identification stuff                                     # 
###########################################################
eIdBase = cms.PSet(
    #Need to put in the variable values.
    detain=cms.double(99.),
    delphi= cms.double(99.),
    sieie= cms.double(99.),
    hoe=cms.double(99.),
    d0vtx=cms.double(-1.),
    dzvtx=cms.double(-1.),
    ooemoop=cms.double(99.),    
    dist=cms.double(0.0),
    dcot=cms.double(0.0),
    hasConv=cms.double(-1.),
    nmisHit=cms.double(99.)
    )


##########################################################################
##usePfIso=True) iso1=chHadIso;iso2=neHadIso;iso3=gamIso;iso4=combIso    # 
##else iso1=tkIso;iso2=ecalIso;iso3=hcalIso;iso4=combIso                 #
##########################################################################
eIsoBase = cms.PSet(
    iso1 = cms.double(9999.),
    iso2 = cms.double(9999.),
    iso3 = cms.double(9999.),
    iso4 = cms.double(9999.)
    )


#######################################
#Base Electron Selector Configuration #
#######################################
eSelBase = cms.PSet(
    ptx= cms.double(20),
    etax= cms.double(2.4),
    ##Identification
    idBrl = eIdBase,
    idEcp = eIdBase,
    ##Isolation: Detector Based Rel Iso.
    coneRad= cms.double(0.3),
    usePfIso= cms.bool(True),
    useCombIso = cms.bool(True),
    useRelIso= cms.bool(True),
    rhoCorr  = cms.bool(True),
    isoBrl = eIsoBase,
    isoEcp = eIsoBase
    )

#####################################
#Some Over Fastidiousness           # 
#####################################


#non-combined relative detector based iso
eSelBaseWithNonCombRelDetIso = eSelBase.clone()

#non-combined absolute detector based iso
eSelBaseWithNonCombAbsDetIso = eSelBase.clone()
eSelBaseWithNonCombAbsDetIso.useRelIso= cms.bool(False)


#combined relative detector based iso
eSelBaseWithCombRelDetIso = eSelBaseWithNonCombRelDetIso.clone()
eSelBaseWithCombRelDetIso.useCombIso = cms.bool(True)

#combined absoulte detector based iso
eSelBaseWithCombAbsRelDetIso = eSelBaseWithCombRelDetIso.clone()
eSelBaseWithCombAbsRelDetIso.useRelIso = cms.bool(True)


















