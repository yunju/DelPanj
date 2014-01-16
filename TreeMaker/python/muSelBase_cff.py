import FWCore.ParameterSet.Config as cms


###########################################################
#identification stuff                                     # 
###########################################################
muIdBase = cms.PSet(
    #Need to put in the variable values.
    dxy= cms.double(99.),
    dz = cms.double(-1),
    normalizedChi2= cms.double(99.),
    trackerHits=cms.double(0.),
    pixelHits=cms.double(0.),
    muonHits=cms.double(0.),
    nMatches=cms.double(0.),
    applyGlbTrk=cms.double(1)
    )


##########################################################################
##usePfIso=True) iso1=chHadIso;iso2=neHadIso;iso3=gamIso;iso4=combIso    # 
##else iso1=tkIso;iso2=ecalIso;iso3=hcalIso;iso4=combIso                 #
##########################################################################
muIsoBase = cms.PSet(
    iso1 = cms.double(9999.),
    iso2 = cms.double(9999.),
    iso3 = cms.double(9999.),
    iso4 = cms.double(9999.)
    )


#######################################
#Base Electron Selector Configuration #
#######################################
muSelBase = cms.PSet(
    pt= cms.double(20),
    eta= cms.double(2.4),
    requireTrigMatch= cms.bool(False),
    ##Isolation: Detector Based Rel Iso.
    coneRad= cms.double(0.4),
    usePfIso= cms.bool(True),
    useCombIso = cms.bool(False),
    useRelIso= cms.bool(False),
    delBCorr = cms.bool(True),
    idPar = muIdBase,
    isoPar = muIsoBase    
    )
