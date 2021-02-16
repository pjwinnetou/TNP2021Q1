#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <dirent.h>
#include <memory>
#include "TreeSetting.h"
#include "cutsAndBinUpsilonV2.h"

using namespace std;

void getGenEff_L2L3_ratio(int idx = 1){

  TString filename = Form("/eos/cms/store/group/phys_heavyions/dileptons/TNPTagAndProbe2018/MC2018/PbPb502TeV/Onia/20191027/OniaTree_HINPbPbAutumn18DR-mva98_JPsi_MC_%d.root",idx);

  string triglist = "triggerList.txt";
  ifstream triglistfile;
  vector<string> TRIGLIST;
  string line;

  triglistfile.open(Form("%s",triglist.c_str()));
  while (getline(triglistfile, line))
  {
    TRIGLIST.push_back(line);
  }
    
  TChain* myTree = new TChain("hionia/myTree");
  myTree->Add(filename.Data());

  SetTree settree_;
  settree_.TreeSetting(myTree,1,0);

  const int sizelist = TRIGLIST.size();
  TChain* OnlineTR[sizelist];
  map<TString,map<TString,vector<double>*>> valVect;
  map<TString,TH1D*> hDen;
  map<TString,TH1D*> hNum;

  const int nEta=4;
  const int nPt[nEta] = {10, 10, 11, 9};
  double etaRange[nEta+1] = {0,1.2,1.8,2.1,2.4};
  double dRL3Cut[nEta] = {0.1,0.2,0.3,0.4};
  double dRL2Cut[nEta] = {0.3,0.4,0.5,0.6};
  double dPtL3Cut[nEta] = {0.3,0.4,0.5,0.6};
  double dPtL2Cut[nEta] = {0.5,0.6,0.7,0.8};

  double dimuFWCut = 3.0;
  double dimuMidCut = 6.5;

  vector<vector<double>> ptBins = {
    {3.5, 4, 4.5, 5, 5.5, 6.5, 8., 10.5, 14, 18, 30.},
    {2.07, 3.0, 3.5, 4, 4.5, 5., 6., 7.5, 10, 15, 30},
    {1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 6.5, 8, 9.5, 13, 20},
    {1.5, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 20},
  };


  int itree=0;
  for(auto& t : TRIGLIST){
    const string tt = t.substr(0, t.size());
    OnlineTR[itree] = new TChain(Form("hltobject/%s",tt.c_str()));
    OnlineTR[itree]->Add(filename.Data());

    int ivar=0;
    for( const auto& var : {"pt", "eta", "phi", "mass"}){     
      valVect[itree][ivar] = 0;
      OnlineTR[itree] -> SetBranchAddress(var,&valVect[itree][ivar]);
      ivar++;
    }
    itree++;
  } 
  int itree_L3Mu5 = itree-1;
  
  for(int ieta=0; ieta<nEta; ieta++){
    const int nSize = ptBins[ieta].size();
    double ptBins_[nSize];
    for(int ib=0;ib<nSize;ib++){ptBins_[ib] = ptBins[ieta][ib];}
    hDen[ieta] = new TH1D(Form("hDen_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum[ieta] = new TH1D(Form("hNum_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen[ieta]->Sumw2();
    hNum[ieta]->Sumw2();
  }

  
  double dRL3MuPlCut = 0.3;
  double dRL2MuPlCut = 0.5;
  double dRL3MuMiCut = 0.3;
  double dRL2MuMiCut = 0.5;
  double dPtL3MuPlCut = 0.3;
  double dPtL3MuMiCut = 0.3;
  double dPtL2MuPlCut = 0.5;
  double dPtL2MuMiCut = 0.5;

  const int nEntries = myTree->GetEntries();

  TLorentzVector* JP_Gen= new TLorentzVector;
  TLorentzVector* muPl_Gen= new TLorentzVector;
  TLorentzVector* muMi_Gen= new TLorentzVector;

  int cBin; double weight;

  for(int iev = 0; iev < nEntries; iev++){
    myTree->GetEntry(iev);
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << myTree->GetEntries() <<  " ("<<(int)(100.*iev/myTree->GetEntries()) << "%)" << endl;
    
    cBin = Centrality;
    weight = findNcoll(Centrality) * Gen_weight;

    if(cBin>=180) continue;
    
    for(int j=0; j<Gen_QQ_size; j++){
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(j);
      muPl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[j]);
      muMi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[j]);

      if(!(IsAcceptanceQQ(muPl_Gen->Pt(), muPl_Gen->Eta()) && IsAcceptanceQQ(muMi_Gen->Pt(), muMi_Gen->Eta())) ) continue;
      if(Gen_mu_charge[Gen_QQ_mupl_idx[j]]*Gen_mu_charge[Gen_QQ_mumi_idx[j]]>0) continue;
      if(abs(JP_Gen->Rapidity())>1.8 && JP_Gen->Pt()<dimuFWCut) continue;
      if(abs(JP_Gen->Rapidity())<1.8 && JP_Gen->Pt()<dimuMidCut) continue;
      
      double pt1 = muPl_Gen->Pt(); double pt2 = muMi_Gen->Pt(); 
      double eta1 = muPl_Gen->Eta(); double eta2 = muMi_Gen->Eta(); 

      for(int ieta=0;ieta<nEta;ieta++){
        if(abs(eta1)<etaRange[ieta+1] && abs(eta1)>=etaRange[ieta]){dRL3MuPlCut = dRL3Cut[ieta]; dRL2MuPlCut = dRL2Cut[ieta]; dPtL3MuPlCut = dPtL3Cut[ieta]; dPtL2MuPlCut = dPtL2Cut[ieta];}
        if(abs(eta2)<etaRange[ieta+1] && abs(eta2)>=etaRange[ieta]){dRL3MuMiCut = dRL3Cut[ieta]; dRL2MuMiCut = dRL2Cut[ieta]; dPtL3MuMiCut = dPtL3Cut[ieta]; dPtL2MuMiCut = dPtL2Cut[ieta];}
      }

      OnlineTR[itree_L3Mu5] ->GetEntry(iev);
      bool muPlTriggerMatchedL3Mu5 = isTriggerMatched(muPl_Gen, valVect, itree_L3Mu5, dRL3MuPlCut, dPtL3MuPlCut); 
      bool muMiTriggerMatchedL3Mu5 = isTriggerMatched(muMi_Gen, valVect, itree_L3Mu5, dRL3MuMiCut, dPtL3MuMiCut); 

      OnlineTR[0]->GetEntry(iev);
      bool muPlL2TriggerMatched = isTriggerMatched(muPl_Gen, valVect, 0, dRL2MuPlCut, dPtL2MuPlCut); 
      bool muMiL2TriggerMatched = isTriggerMatched(muMi_Gen, valVect, 0, dRL2MuMiCut, dPtL2MuMiCut); 
      OnlineTR[1]->GetEntry(iev);
      bool muPlL3TriggerMatched = isTriggerMatched(muPl_Gen, valVect, 1, dRL3MuPlCut, dPtL3MuPlCut); 
      bool muMiL3TriggerMatched = isTriggerMatched(muMi_Gen, valVect, 1, dRL3MuMiCut, dPtL3MuMiCut); 

      for(int ieta=0;ieta<nEta;ieta++){
        if(abs(eta1)<etaRange[ieta+1] && abs(eta1)>=etaRange[ieta] && muMiTriggerMatchedL3Mu5){
          if(muPlL2TriggerMatched) hDen[ieta] -> Fill(pt1,weight);
          if(muPlL2TriggerMatched && muPlL3TriggerMatched) hNum[ieta]->Fill(pt1,weight);
        }
        if(abs(eta2)<etaRange[ieta+1] && abs(eta2)>=etaRange[ieta] && muPlTriggerMatchedL3Mu5){
          if(muMiL2TriggerMatched) hDen[ieta] -> Fill(pt1,weight);
          if(muMiL2TriggerMatched && muMiL3TriggerMatched) hNum[ieta]->Fill(pt1,weight);
        }
      }
    }
  }


  TFile* wf = new TFile(Form("output_L2L3_ratio_%d.root",idx),"RECREATE");
  wf->cd();

  map<TString,TH1D*> hRatio;
  for(int ieta=0;ieta<nEta;ieta++){
    hRatio[ieta] = (TH1D*) hNum[ieta]->Clone(Form("hEff_eta%d",ieta));
    hRatio[ieta]->Divide(hDen[ieta]);
    hRatio[ieta]->Write();
    hDen[ieta]->Write();
    hNum[ieta]->Write();
  }
}
