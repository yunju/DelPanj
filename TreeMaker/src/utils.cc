#include "DelPanj/TreeMaker/interface/utils.h"
//#include "TLorentzVector.h"
//#include <vector>
//#include <string>
//#include <map>

/*
class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.pt() > j.pt());
  }
};
*/
TLorentzVector Part2LorVec(reco::Candidate& cand){
  TLorentzVector* l = new TLorentzVector();
  l->SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass());
  return (*l);
}

//When Selectors are fully developed,this must move to baseSelector class.
int PassAll(std::map<std::string, bool> cutrecd){
  std::map<std::string, bool>::iterator iter= cutrecd.begin();
  bool decision =1 ;
  for(;iter!=cutrecd.end();iter++){
//     std::cout<<"-->"<<iter->first<<"\t"<<iter->second<<std::endl;	   
    decision = decision&&iter->second;     
  }
  return (int)decision;
}


int PassAllBut(std::string tag, std::map<std::string, bool> cutrecd){
  std::map<std::string, bool>::iterator iter= cutrecd.begin();
  bool decision =1 ;
  for(;iter!=cutrecd.end();iter++){
    if(iter->first==tag)continue;
    decision = decision&&iter->second;
  }
 return (int)decision;  
}



