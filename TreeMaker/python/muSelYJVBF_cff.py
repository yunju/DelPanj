import FWCore.ParameterSet.Config as cms


muIdBaseYJVBF = cms.PSet(
    #Need to put in the variable values.
    dxy= cms.double(0.2),
    dz = cms.double(0.5),
    normalizedChi2= cms.double(10.),
    trackerHits=cms.double(5.),
    pixelHits=cms.double(0.),
    muonHits=cms.double(0.),
    nMatches=cms.double(1),
    applyGlbTrk=cms.double(1)# not using now
    )
##########################################################################
##usePfIso=True) iso1=chHadIso;iso2=neHadIso;iso3=gamIso;iso4=combIso    # 
##else iso1=tkIso;iso2=ecalIso;iso3=hcalIso;iso4=combIso                 #
##########################################################################
muIsoBaseYJVBF = cms.PSet(
    iso1 = cms.double(9999.),
    iso2 = cms.double(9999.),
    iso3 = cms.double(9999.),
    iso4 = cms.double(0.12)
    )

muSelYJVBF = cms.PSet(
    pt= cms.double(20),
    eta= cms.double(2.4),
    requireTrigMatch= cms.bool(False),
    ##Isolation: Detector Based Rel Iso.
    coneRad= cms.double(0.4), #this parameter is not used anymore.
    usePfIso= cms.bool(True),
    useCombIso = cms.bool(True),
    useRelIso= cms.bool(True),
    delBCorr = cms.bool(True),
    idPar = muIdBaseYJVBF,
    isoPar = muIsoBaseYJVBF    
    )


