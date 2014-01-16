#include "DelPanj/TreeMaker/interface/patMuonTree.h"

patMuonTree::patMuonTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig): mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")){
 tree_=tree; 
  patMuonLabel_ = iConfig.getParameter<edm::InputTag>("patMuonsPY");
  beamSpotLabel_ = iConfig.getParameter<edm::InputTag> ("beamSpotLabel");
  
SetBranches();
}


patMuonTree::~patMuonTree(){
  delete tree_;

} 

double
patMuonTree::IsoDeposit(const pat::Muon Mu,std::string type,double radius){
  double SumPt = 0;
  if(std::string::npos != type.find("char") && std::string::npos != type.find("Had"))
    SumPt=Mu.isoDeposit(pat::PfChargedHadronIso)
      ->depositAndCountWithin(radius).first;

  if(std::string::npos != type.find("neu") && std::string::npos != type.find("Had"))
    SumPt=Mu.isoDeposit(pat::PfNeutralHadronIso)
      ->depositAndCountWithin(radius).first;


  if(std::string::npos != type.find("phot") || std::string::npos != type.find("gam"))
    SumPt=Mu.isoDeposit(pat::PfGammaIso)
      ->depositAndCountWithin(radius).first;

  return SumPt;
}

void
patMuonTree::Fill(const edm::Event& iEvent){
Clear();
//  cout<<"========Event :"<<iEvent.id().event()<<endl;
  edm::Handle<pat::MuonCollection> patMuonHandle;
  if(not iEvent.getByLabel(patMuonLabel_,patMuonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patMuonLabel_<<std::endl; exit(0);}

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
  reco::BeamSpot beamSpot = (*beamSpotHandle);
  
  pat::MuonCollection muColl(*(patMuonHandle.product()));
  std::sort(muColl.begin(),muColl.end(),PtGreater());

  nmu=muColl.size();
  pat::MuonCollection::const_iterator mu;

/*
  pat::MuonCollection::const_iterator mu1;
  pat::MuonCollection::const_iterator mu2;
   
 
  int nmu1=0;
  for(mu1=muColl.begin(); mu1!=muColl.end(); mu1++){
   nmu1++;
  // cout<<nmu1<<" ---1st  "<<endl;
   TLorentzVector Vmu1;
   Vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),mu1->mass()); 
   int nmu2=0; 
   for(mu2=muColl.begin(); mu2!=muColl.end(); mu2++){
     
      //cin.get();  
      nmu2++;
    //   cout<<nmu2<<" ---2nd  "<<endl;
      if(nmu1>=nmu2) continue;  
    //  cout<<nmu1<<" ---  "<<nmu2<<endl;
      
     TLorentzVector Vmu2;     
      Vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),mu2->mass());
      TLorentzVector Vmumu;
      Vmumu=Vmu1+Vmu2; 
      cout<<nmu1<<" "<<nmu2<<" mll:"<<Vmumu.M()<<endl;
   
   }
 }

*/
 
int nMu=0;
  for(mu=muColl.begin(); mu!=muColl.end(); mu++){

   nMu++; 
  // const pat::Muon* myMuo
 // = mu;
  std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*mu);
  int passOrNot = PassAll(Pass);
   //cout <<"passOrNot: "<<passOrNot<<endl;
    patMuonID_.push_back(passOrNot);
   // cout<<"UserDataMu Pt :"<<mu->pt()<<endl;
    patMuonPt_.push_back(mu->pt());
    patMuonEta_.push_back(mu->eta());
    patMuonPhi_.push_back(mu->phi());
    patMuonM_.push_back(mu->mass());
    patMuonTrkIso_.push_back(mu->trackIso());
    patMuonHcalIso_.push_back(mu->hcalIso());
    patMuonEcalIso_.push_back(mu->ecalIso());
    patMuonCharge_.push_back(mu->charge());
    patMuonNumChambers_.push_back(mu->numberOfChambers());
    patMuonNumMatches_.push_back(mu->numberOfMatches());
    patMuonStationMask_.push_back(mu->stationMask());
    
    bool isTrackMuon =mu->isTrackerMuon();
    bool isGlobalMuon =mu->isGlobalMuon();
    bool isStandAloneMuon =mu->isStandAloneMuon();
  
    isGloMu_.push_back(isGlobalMuon);
    isSTDMu_.push_back(isStandAloneMuon);
    isTraMu_.push_back(isTrackMuon); 
/*
    reco::TrackRef muTrk;
    if(isStandAloneMuon)
    {    
    muTrk = mu->standAloneMuon();
    MuType_.push_back(1);
    }   
  else if(isTrackMuon)
   { 
     muTrk = mu->track();
     MuType_.push_back(2);
   }
    else if(isGlobalMuon)
    {
      muTrk = mu->globalTrack();
       MuType_.push_back(3);
    }  
*/
    //patMuonSumTrkPt_.push_back(muTrk->);
    //patMuonNumSegments_.push_back(muTrkTrk->);
   
    //cout<<"which muon : "<<isGlobalMuon<<" "<<isStandAloneMuon<<" "<<isTrackMuon <<endl;
    

    if(isGlobalMuon)
    {
  //  if(passOrNot&&fabs(mu->eta())<2.4) cout<<"UserDataMu Pt :"<<iEvent.id().event()<<" "<<nMu<<" ;pass:"<<passOrNot <<" ;pt:"<<mu->pt()<<" "<<mu->eta()<<" Charge:"<<mu->charge()<<endl;
    
    patMuonChi2Ndoff_.push_back(mu->globalTrack()->normalizedChi2());
    patMuonNhits_.push_back(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
    patMuonNMatStat_.push_back(mu->numberOfMatchedStations()); 
    patMuondB_.push_back(mu->dB()); 
    //patMuonDxy_.push_back(mu->innerTrack()->dxy(vertex->position()));
    patMuonDz_.push_back(mu->userFloat("dzVtx"));
    patMuonNumPix_.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
    patMuonNumTlayer_.push_back(mu->track()->hitPattern().trackerLayersWithMeasurement());
    }
    else
    { 
      patMuonChi2Ndoff_.push_back(-999);
      patMuonNhits_.push_back(-999);
      patMuonNMatStat_.push_back(-999); 
      patMuondB_.push_back(-999); 
      //patMuonDxy_.push_back(mu->innerTrack()->dxy(vertex->position()));
      patMuonDz_.push_back(-999);
      patMuonNumPix_.push_back(-999);
      patMuonNumTlayer_.push_back(-999);

    }

    patMuonChHadSumPt04_.push_back(mu->pfIsolationR04().sumChargedHadronPt);
    patMuonNeHadSumPt04_.push_back(mu->pfIsolationR04().sumNeutralHadronEt);
    patMuonGamSumPt04_.push_back(mu->pfIsolationR04().sumPhotonEt);
    patMuonPUSumPt04_.push_back(mu->pfIsolationR04().sumPUPt);
/*
    patMuonQoverP_.push_back(muTrk->qoverp());
    isTraMu_patMuonTheta_.push_back(muTrk->theta());
    patMuonLambda_.push_back(muTrk->lambda());
    patMuonDxy_.push_back(muTrk->dxy());
    patMuonD0_.push_back(muTrk->d0());
    patMuonDsz_.push_back(muTrk->dsz());
    patMuonDxyBS_.push_back(muTrk->dxy(beamSpot.position()));
    patMuonDzBS_.push_back(muTrk->dsz(beamSpot.position()));
    patMuonDszBS_.push_back(muTrk->dz(beamSpot.position()));
    patMuonVx_.push_back(muTrk->vx());
    patMuonVy_.push_back(muTrk->vy());
    patMuonVz_.push_back(muTrk->vz());
*/
//    bool useGenInfo = 1;
  //  if(useGenInfo){    
      reco::GenParticleRef genMu = mu->genParticleRef();
      if(genMu.isNonnull()){
	patMuonGenMatchPt_.push_back(genMu->pt());
	patMuonGenMatchEta_.push_back(genMu->eta());
	patMuonGenMatchPhi_.push_back(genMu->phi());
	//  patMuonGenMotherId_.push_back(genMu->);
	patMuonGenMatchCharge_.push_back(genMu->charge());
    }

//==================
      else{
        patMuonGenMatchPt_.push_back(0);
        patMuonGenMatchEta_.push_back(0);
        patMuonGenMatchPhi_.push_back(0);
        //  patMuonGenMotherId_.push_back(genMu->);
        patMuonGenMatchCharge_.push_back(0);
     }


   // }
//   //particle based isolation.

    patMuonChHadSumPt03_.push_back(IsoDeposit(*mu,"charHad",0.3));
    patMuonNeHadSumPt03_.push_back(IsoDeposit(*mu,"neuHad",0.3));
    patMuonGamSumPt03_.push_back(IsoDeposit(*mu,"gam",0.3));
    
    
    patMuonChHadSumPt05_.push_back(IsoDeposit(*mu,"charHad",0.5));
    patMuonNeHadSumPt05_.push_back(IsoDeposit(*mu,"neuHad",0.5));
    patMuonGamSumPt05_.push_back(IsoDeposit(*mu,"gam",0.5));
  }//muon Loop
}
void
patMuonTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}



void
patMuonTree::AddBranch(std::vector<double>* vec, std::string name){
  //std::string brName="patMuon"+name;
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void 
patMuonTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(name.c_str(),x,(brName+"/D").c_str());
}

void 
patMuonTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName="patMuon"+name;
  tree_->Branch(name.c_str(),vec);
}





void 
patMuonTree::SetBranches(){
AddBranch(&nmu,"NumMu");
AddBranch(&patMuonPt_, "patMuonPt_");
AddBranch(&patMuonEta_, "patMuonEta_");
AddBranch(&patMuonPhi_, "patMuonPhi_");
AddBranch(&patMuonM_, "patMuonM_");
AddBranch(&patMuonSumTrkPt_, "patMuonSumTrkPt_");
AddBranch(&patMuonTrkIso_, "patMuonTrkIso_");
AddBranch(&patMuonHcalIso_, "patMuonHcalIso_");
AddBranch(&patMuonEcalIso_, "patMuonEcalIso_");
AddBranch(&patMuonCharge_, "patMuonCharge_");
AddBranch(&patMuonNumChambers_, "patMuonNumChambers_");
AddBranch(&patMuonNumMatches_, "patMuonNumMatches_");
AddBranch(&patMuonStationMask_, "patMuonStationMask_");
AddBranch(&patMuonNumSegments_, "patMuonNumSegments_");
AddBranch(&patMuonChi2Ndoff_, "patMuonChi2Ndoff_");
AddBranch(&patMuonNhits_, "patMuonNhits_");
AddBranch(&patMuonQoverP_, "patMuonQoverP_");
AddBranch(&patMuonTheta_, "patMuonTheta_");
AddBranch(&patMuonLambda_, "patMuonLambda_");
AddBranch(&patMuonD0_, "patMuonD0_");
AddBranch(&patMuonDsz_, "patMuonDsz_");
AddBranch(&patMuonDxyBS_, "patMuonDxyBS_");
AddBranch(&patMuonDzBS_, "patMuonDzBS_");
AddBranch(&patMuonDszBS_, "patMuonDszBS_");
AddBranch(&patMuonVx_, "patMuonVx_");
AddBranch(&patMuonVy_, "patMuonVy_");
AddBranch(&patMuonVz_, "patMuonVz_");
AddBranch(&patMuonGenMatchPt_, "patMuonGenMatchPt_");
AddBranch(&patMuonGenMatchEta_, "patMuonGenMatchEta_");
AddBranch(&patMuonGenMatchPhi_, "patMuonGenMatchPhi_");
AddBranch(&patMuonGenMotherId_, "patMuonGenMotherId_");
AddBranch(&patMuonGenMatchCharge_, "patMuonGenMatchCharge_");
AddBranch(&patMuonChHadSumPt03_, "patMuonChHadSumPt03_");
AddBranch(&patMuonNeHadSumPt03_, "patMuonNeHadSumPt03_");
AddBranch(&patMuonGamSumPt03_, "patMuonGamSumPt03_");

AddBranch(&patMuonChHadSumPt05_, "patMuonChHadSumPt05_");
AddBranch(&patMuonNeHadSumPt05_, "patMuonNeHadSumPt05_");
AddBranch(&patMuonGamSumPt05_, "patMuonGamSumPt05_");

AddBranch(&patMuonID_,"patMuonID_");
AddBranch(&isGloMu_,"isGloMu_");
AddBranch(&isSTDMu_,"isSTDMu_");
AddBranch(&isTraMu_,"isTraMu_");
AddBranch(&patMuonNMatStat_,"patMuonNMatStat_");
AddBranch(&patMuondB_,"patMuondB_");
//AddBranch(&patMuonDxy_,"patMuonDxy_");
AddBranch(&patMuonDz_, "patMuonDz_");
AddBranch(&patMuonNumPix_, "patMuonNumPix_");
AddBranch(&patMuonNumTlayer_, "patMuonNumTlayer_");

AddBranch(&patMuonChHadSumPt04_, "patMuonChHadSumPt04_");
AddBranch(&patMuonNeHadSumPt04_, "patMuonNeHadSumPt04_");
AddBranch(&patMuonGamSumPt04_, "patMuonGamSumPt04_");
AddBranch(&patMuonPUSumPt04_,"patMuonPUSumPt04_");
}


void
patMuonTree::Clear(){
patMuonPt_.clear();
patMuonEta_.clear();
patMuonPhi_.clear();
patMuonM_.clear();
patMuonSumTrkPt_.clear();
patMuonTrkIso_.clear();
patMuonHcalIso_.clear();
patMuonEcalIso_.clear();
patMuonCharge_.clear();
patMuonNumChambers_.clear();
patMuonNumMatches_.clear();
patMuonStationMask_.clear();
patMuonNumSegments_.clear();
patMuonChi2Ndoff_.clear();
patMuonNhits_.clear();
patMuonQoverP_.clear();
patMuonTheta_.clear();
patMuonLambda_.clear();
patMuonD0_.clear();
patMuonDsz_.clear();
patMuonDxyBS_.clear();
patMuonDzBS_.clear();
patMuonDszBS_.clear();
patMuonVx_.clear();
patMuonVy_.clear();
patMuonVz_.clear();
patMuonGenMatchPt_.clear();
patMuonGenMatchEta_.clear();
patMuonGenMatchPhi_.clear();
patMuonGenMotherId_.clear();
patMuonGenMatchCharge_.clear();
patMuonChHadSumPt03_.clear();
patMuonNeHadSumPt03_.clear();
patMuonGamSumPt03_.clear();
patMuonChHadSumPt05_.clear();
patMuonNeHadSumPt05_.clear();
patMuonGamSumPt05_.clear();

patMuonID_.clear();
isGloMu_.clear();
isSTDMu_.clear();
isTraMu_.clear();
patMuonNMatStat_.clear();
patMuondB_.clear();
//patMuonDxy_.clear();
patMuonDz_.clear();
patMuonNumPix_.clear();
patMuonNumTlayer_.clear();

patMuonChHadSumPt04_.clear();
patMuonNeHadSumPt04_.clear();
patMuonGamSumPt04_.clear();
patMuonPUSumPt04_.clear();


nmu=0;
}
