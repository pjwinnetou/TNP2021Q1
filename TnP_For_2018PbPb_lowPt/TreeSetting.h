#ifndef TreeSetting_h
#define TreeSetting_h

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"

using namespace std;

const long int maxBranchSize = 100000;

// import the tree to the RooDataSet
UInt_t          runNb;
UInt_t          eventNb, LS;
float           zVtx;
Int_t           Centrality;
Float_t         SumET_HF;
ULong64_t       HLTriggers;
Int_t           Reco_QQ_size;
Int_t           Reco_mu_size;
TClonesArray    *Reco_QQ_4mom;
TClonesArray    *Reco_mu_4mom;
ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_runNb;   //!
TBranch        *b_eventNb;   //!
TBranch        *b_LS;
TBranch        *b_zVtx;   //!
TBranch        *b_Centrality;   //!
TBranch        *b_SumET_HF;
TBranch        *b_HLTriggers;   //!
TBranch        *b_Reco_QQ_size;   //!
TBranch        *b_Reco_QQ_4mom;   //!
TBranch        *b_Reco_mu_size;   //!
TBranch        *b_Reco_mu_4mom;   //!
TBranch        *b_Reco_QQ_trig;   //!
TBranch        *b_Reco_mu_trig;   //!
TBranch        *b_Reco_QQ_VtxProb;   //!

Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_mu_highPurity;   //!

//  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nTrkHits;   //!
Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_normChi2_global;   //!
Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nMuValHits;   //!
Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_StationsMatched;   //!
Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_dxy;   //!
TBranch        *b_Reco_mu_dxyErr;   //!
Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_dz;   //!
TBranch        *b_Reco_mu_dzErr;   //!
Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nTrkWMea;   //!
Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_TMOneStaTight;   //!
Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nPixWMea;   //!
Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_sign;   //!

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;

Float_t         Reco_QQ_ctau[maxBranchSize];
Float_t         Reco_QQ_ctau3D[maxBranchSize];
TBranch        *b_Reco_QQ_ctau;
TBranch        *b_Reco_QQ_ctau3D;

Float_t         Reco_QQ_MassErr[maxBranchSize];
Float_t         Reco_mu_normChi2_inner[maxBranchSize];
Float_t         Reco_mu_ptErr_inner[maxBranchSize];
TBranch        *b_Reco_QQ_MassErr;
TBranch        *b_Reco_mu_normChi2_inner;
TBranch        *b_Reco_mu_ptErr_inner;


/////////////////////////////////////////
////// Gen QQ 
/////////////////////////////////////////
Int_t           Gen_QQ_size;
Int_t           Gen_QQ_type[maxBranchSize];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_4mom;
TClonesArray    *Gen_mu_4mom;
TBranch        *b_Gen_QQ_size;   //!
TBranch        *b_Gen_QQ_type;   //!
TBranch        *b_Gen_QQ_4mom;   //!
TBranch        *b_Gen_mu_4mom;   //!
Int_t Reco_mu_whichGen[maxBranchSize];
TBranch *b_Reco_mu_whichGen;
Float_t Gen_weight;
TBranch *b_Gen_weight;

Int_t Gen_QQ_mupl_idx[maxBranchSize];
Int_t Gen_QQ_mumi_idx[maxBranchSize];
TBranch *b_Gen_QQ_mupl_idx;
TBranch *b_Gen_QQ_mumi_idx;

Int_t Gen_mu_charge[maxBranchSize];
TBranch *b_Gen_mu_charge;



//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*Conversion~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//

const int MaxTrackBranchSize = 5000;
Int_t          Reco_trk_size;
TClonesArray   *Reco_trk_4mom;
TClonesArray   *Reco_trk_vtx;
Int_t          Reco_trk_charge[MaxTrackBranchSize];
Float_t        Reco_trk_dxyError[MaxTrackBranchSize];
Float_t        Reco_trk_dzError[MaxTrackBranchSize];

TBranch   *b_Reco_trk_size;
TBranch   *b_Reco_trk_4mom;
TBranch   *b_Reco_trk_vtx;
TBranch   *b_Reco_trk_charge;
TBranch   *b_Reco_trk_dxyError;
TBranch   *b_Reco_trk_dzError;


class SetTree
{
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree, bool isMC, bool isConversion);
    Bool_t SoftMuIdCut(int irqq);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* tree, bool isMC, bool isConversion)
{
  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;

  tree->SetBranchAddress("runNb", &runNb, &b_runNb);
  tree->SetBranchAddress("LS", &LS, &b_LS);
  tree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  tree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  tree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  tree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  tree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  tree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  tree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  tree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  tree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);
  tree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  tree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  tree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  tree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  tree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  tree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  tree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  tree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  tree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  tree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  tree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  tree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  tree->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau,&b_Reco_QQ_ctau);
  tree->SetBranchAddress("Reco_QQ_ctau3D",Reco_QQ_ctau3D,&b_Reco_QQ_ctau3D);
  tree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  tree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
  tree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
  tree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  tree->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
  tree->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  tree->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);



  if (isMC) { 
    Gen_QQ_4mom = 0;
    Gen_mu_4mom = 0;
    tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    tree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
    tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    tree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
    tree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    tree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
    tree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
    tree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
    tree->SetBranchAddress("Gen_mu_charge",Gen_mu_charge, &b_Gen_mu_charge);
  }

  if(isConversion) {
    Reco_trk_4mom=0;
    Reco_trk_vtx=0;
    tree->SetBranchAddress("Reco_trk_size", &Reco_trk_size, &b_Reco_trk_size);
    tree->SetBranchAddress("Reco_trk_4mom", &Reco_trk_4mom, &b_Reco_trk_4mom);
    tree->SetBranchAddress("Reco_trk_vtx", &Reco_trk_vtx, &b_Reco_trk_vtx);
    tree->SetBranchAddress("Reco_trk_charge", Reco_trk_charge, &b_Reco_trk_charge);
    tree->SetBranchAddress("Reco_trk_dxyError", Reco_trk_dxyError, &b_Reco_trk_dxyError);
    tree->SetBranchAddress("Reco_trk_dzError", Reco_trk_dzError, &b_Reco_trk_dzError);
  }

};

Bool_t SetTree::SoftMuIdCut(int irqq)
{
  bool passMuonTypePl = true;
  passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
  passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

  bool passMuonTypeMi = true;
  passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
  passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

  bool muplSoft = (  //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
      (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
      (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
      (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
      (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
      passMuonTypePl        //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
      ) ; 

  bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
      (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
      (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
      (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
      (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
      passMuonTypeMi       //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
      ) ; 

  return (muplSoft && mumiSoft);
}  



#endif 
