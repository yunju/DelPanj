#include "DelPanj/TreeMaker/interface/jetTree.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TMath.h"
typedef math::XYZTLorentzVector LorentzVector;

jetTree::jetTree(std::string desc, TTree* tree, const edm::ParameterSet& iConfig)
{
  tree_=tree;  
  JetLabel_ = iConfig.getParameter<edm::InputTag>("JetsPY");
  genPartLabel_ = iConfig.getParameter<edm::InputTag>("genPartLabel");
  SetBranches();
}


jetTree::~jetTree(){
  delete tree_;
}

bool jetTree::passLooseJetID(const pat::Jet* recjet)
{
  double eta = recjet->eta();
  if(recjet->getPFConstituents().size() <= 1)return false;                                                           
  if(recjet->neutralHadronEnergyFraction() >= 0.99)return false;
  if(recjet->neutralEmEnergyFraction() >= 0.99)return false;
  //   // for the tracker region
  if(fabs(eta)<2.4 && recjet->chargedHadronEnergyFraction()<= 0.0)return false;
  if(fabs(eta)<2.4 && recjet->chargedEmEnergyFraction() >= 0.99)return false;
  if(fabs(eta)<2.4 && recjet->chargedMultiplicity() <= 0)return false;
  return true;

}

void
jetTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup){
  Clear();
   bool isData = iEvent.isRealData();  
      double dummy = -99999.0;
  
edm::Handle<std::vector<pat::Jet> > JetHandle;
  if(not iEvent.getByLabel(JetLabel_,JetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<JetLabel_<<std::endl; exit(0);}


  const std::vector<pat::Jet>* jets = JetHandle.product();
  std::vector<pat::Jet>::const_iterator jet =jets->begin();

  for(;jet!=jets->end();jet++){

    //Stuff common for all jets.
    JetPt_.push_back(jet->pt());
    JetEta_.push_back(jet->eta());
    JetPhi_.push_back(jet->phi());
    JetM_.push_back(jet->mass());
    JetEn_.push_back(jet->energy());

     
    if ( !passLooseJetID(&*jet) )
    { 
         JetID_.push_back(0); 
         //cout<<"Fail: "<<jet->eta()<<"jet->chargedMultiplicity() "<<jet->chargedMultiplicity()<<" "<<jet->neutralEmEnergyFraction ()<<" "<<jet->chargedHadronEnergyFraction()<<endl; 
    //cin.get();   
    }
    else
    {
         JetID_.push_back(1); 
        // cout<<jet->eta()<<"jet->chargedMultiplicity() "<<jet->chargedMultiplicity()<<" "<<jet->neutralEmEnergyFraction ()<<" "<<jet->chargedHadronEnergyFraction()<<endl;
          //cin.get();  
    }    
    //jetProb btag
  JetJetProb_.push_back(jet->bDiscriminator("jetProbabilityBJetTags"));

    if(jet->isPFJet()){      
      JetCharMulti_.push_back(jet->chargedMultiplicity());
      JetNeutEmEFr_.push_back(jet->neutralEmEnergyFraction ());
      JetCharHadEFr_.push_back(jet->chargedHadronEnergyFraction ());
      JetNeutHadEFr_.push_back(jet->neutralHadronEnergyFraction ());
      JetCharEmEFr_.push_back(jet->chargedEmEnergyFraction ());
    }
    else{
      JetCharMulti_.push_back(dummy);
      JetNeutEmEFr_.push_back(dummy);
      JetCharHadEFr_.push_back(dummy);
      JetNeutHadEFr_.push_back(dummy);
      JetCharEmEFr_.push_back(dummy);
    }
     
     if(!isData)
      {
        edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle; 
        iEvent.getByLabel(genPartLabel_,genParticlesHandle);
   
    //    cout<<"begin"<<endl; 
        if ( genParticlesHandle.isValid()&&jet->genParton()) 
        {
             JetGenPartonID_.push_back(jet->genParton()->pdgId());
             JetGenPartonEn_.push_back(jet->genParton()->energy());
             JetGenPartonPt_.push_back(jet->genParton()->pt());
             JetGenPartonEta_.push_back(jet->genParton()->eta());
             JetGenPartonPhi_.push_back(jet->genParton()->phi());
             JetGenPartonStatus_.push_back(jet->genParton()->status());
             JethasGenParton_.push_back(1);
             if(jet->genParton()->mother())
             {
               JetGenPartonU1ID_.push_back(jet->genParton()->mother()->pdgId());
               if(jet->genParton()->mother()->mother())
               {
                 JetGenPartonU2ID_.push_back(jet->genParton()->mother()->mother()->pdgId());
                 if(jet->genParton()->mother()->mother()->mother())     
                 {
                   JetGenPartonU3ID_.push_back(jet->genParton()->mother()->mother()->mother()->pdgId()); 
                 }//has u3
                 else
                 {
                   JetGenPartonU3ID_.push_back(dummy);       
                 }//no u3       

               }//has u2
               else
               {
                 JetGenPartonU2ID_.push_back(dummy);
                 JetGenPartonU3ID_.push_back(dummy);
               }//no u2
   

             }//has u1
             else
             {
             JetGenPartonU1ID_.push_back(dummy);
             JetGenPartonU2ID_.push_back(dummy);
             JetGenPartonU3ID_.push_back(dummy);
             }//no u1 
   


         }//has genparton
         else
         {
             JetGenPartonID_.push_back(dummy);
             JetGenPartonEn_.push_back(dummy);
             JetGenPartonPt_.push_back(dummy);
             JetGenPartonEta_.push_back(dummy);
             JetGenPartonPhi_.push_back(dummy);
             JetGenPartonStatus_.push_back(dummy);
             
             JetGenPartonU1ID_.push_back(dummy);
             JetGenPartonU2ID_.push_back(dummy);
             JetGenPartonU3ID_.push_back(dummy);
             JethasGenParton_.push_back(0);  
         }//no genparton

  //    cout<<"end"<<endl;
         
    }//is MC

 

  }//jet loop
}

void
jetTree::SetBranches(){
  
  AddBranch(&JetPt_, "JetPt_");
  AddBranch(&JetEta_, "JetEta_");
  AddBranch(&JetPhi_, "JetPhi_");
  AddBranch(&JetM_, "JetM_");
  AddBranch(&JetEn_, "JetEn_");
  
  AddBranch(&JetCharMulti_, "JetCharMulti_");
  AddBranch(&JetNeutEmEFr_, "JetNeutEmEFr_");
  AddBranch(&JetCharHadEFr_, "JetCharHadEFr_");
  AddBranch(&JetNeutHadEFr_, "JetNeutHadEFr_");
  AddBranch(&JetCharEmEFr_, "JetCharEmEFr_");
  AddBranch(&JetID_,"JetID_");  
  AddBranch(&JetJetProb_,"JetJetProb_");

   AddBranch(&JetNConstituents_,"JetNConstituents_");
    AddBranch(&JetGenPartonID_,"JetGenPartonID_");
    AddBranch(&JetGenPartonEn_,"JetGenPartonEn_");
    AddBranch(&JetGenPartonPt_,"JetGenPartonPt_");
    AddBranch(&JetGenPartonEta_,"JetGenPartonEta_") ;
    AddBranch(&JetGenPartonPhi_,"JetGenPartonPhi_") ;
    AddBranch(&JetGenPartonStatus_,"JetGenPartonStatus_");
    AddBranch(&JetGenPartonU3ID_,"JetGenPartonU3ID_");
    AddBranch(&JetGenPartonU2ID_,"JetGenPartonU2ID_");
    AddBranch(&JetGenPartonU1ID_,"JetGenPartonU1ID_");
    AddBranch(&JethasGenParton_,"JethasGenParton_") ;









}


void
jetTree::Clear(){

  JetPt_.clear();
  JetEta_.clear();
  JetPhi_.clear();
  JetM_.clear();
  JetEn_.clear();

  JetCharMulti_.clear();
  JetNeutEmEFr_.clear();
  JetCharHadEFr_.clear();
  JetNeutHadEFr_.clear();
  JetCharEmEFr_.clear();
 JetID_.clear();
 JetJetProb_.clear();
  
JetNConstituents_.clear();
 JetGenPartonID_.clear();
 JetGenPartonEn_.clear();
 JetGenPartonPt_.clear();
 JetGenPartonEta_.clear();
 JetGenPartonPhi_.clear();
 JetGenPartonStatus_.clear();
 JetGenPartonU3ID_.clear();
 JetGenPartonU2ID_.clear();
 JetGenPartonU1ID_.clear();
 JethasGenParton_.clear();













}
