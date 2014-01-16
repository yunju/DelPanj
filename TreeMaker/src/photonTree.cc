
#include "DelPanj/TreeMaker/interface/photonTree.h"

photonTree::photonTree(std::string name, TTree* tree, const pset& iConfig):baseTree(name,tree){
  photonLabel_  = iConfig.getParameter<edm::InputTag> ("photonLabel");
  usePFObjects_ = iConfig.getParameter<bool> ("usePFlow");
  SetBranches();
}

photonTree::~photonTree(){
  delete tree_;
}

void photonTree::Fill(const edm::Event& iEvent){
  Clear();
  //fetch the input collection
  edm::Handle<std::vector<pat::Photon> > photonHandle;
  if(not iEvent.getByLabel(photonLabel_,photonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<photonLabel_<<std::endl; 
    exit(0);
  }  
  pat::PhotonCollection phColl(*(photonHandle.product()));

 //sort the objects by transverse momentum
  std::sort(phColl.begin(),phColl.end(),PtGreater());

  nPhoton_ = phColl.size();
  pat::PhotonCollection::const_iterator ph;
  for(ph=phColl.begin(); ph!=phColl.end(); ph++){
    photonPt_.push_back(ph->pt());
    photonEta_.push_back(ph->eta());
    photonPhi_.push_back(ph->eta());
  }  

}

void photonTree::SetBranches(){
  AddBranch(&nPhoton_  ,"NumPh_");
  AddBranch(&photonPt_ ,"Pt");
  AddBranch(&photonEta_,"Eta");
  AddBranch(&photonPhi_,"Phi");
}

void photonTree::Clear(){
  photonPt_.clear();
  photonEta_.clear();
  photonPhi_.clear();
  }
