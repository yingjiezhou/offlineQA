#include "syang/headers.h"
#include "syang/function.C"
//Draws event QA plots
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <string>

#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include <TLatex.h>
#include <THistPainter.h>
#include <TAttLine.h>
#include <TLegend.h>
#include <TStyle.h>

typedef map<Int_t, Int_t> IntMap;

using namespace std;
IntMap mRunId2RunIdx;
IntMap mRunIdx2RunId;

Bool_t Init(TString runlist="0911_production_7p7GeV_2020_runnumber_DD.dat");
Int_t  grabRunIdx(Int_t runId);
Int_t  grabRunId(Int_t runIdx);

void drawDay(Int_t runStart=21246012, Int_t runStop=21251016, TString runlist="0911_production_7p7GeV_2020_runnumber_DD.dat"){
  
  if(!Init(runlist)) cout<<"Failed to initialize code !";
  
  TString outfileName;
  outfileName = runlist;
  outfileName.ReplaceAll("_runnumber_DD.dat", "");
  
  
  Int_t binLow,binHi;
  IntMap::iterator iter = mRunId2RunIdx.find(runStart);
  if(iter != mRunId2RunIdx.end()){
    binLow = iter->second+1;
  }
  else{
    cout<<"Can not find the start runId"<<" in the run list !"<<endl;
    return;
  }
  
  iter = mRunId2RunIdx.find(runStop);
  if(iter != mRunId2RunIdx.end()){
    binHi = iter->second+1;
  }
  else{
    cout<<"Can not find the stop runId"<<" in the run list !"<<endl;
    return;
  }
  
  cout<<"start run bin:"<<binLow<<" \t stop run bin"<<binHi<<endl;
  //===============================================================
  
  vector<double> runNum;
  
  vector<double> runCount;
  vector<double> runCount_err;
  int run_count = 0;
  
  vector<double> newDay;
  vector<TString> dayN;
  
  
  TString dayStart;
  TString dayEnd;
  
  int dayI;
  int dayPrev = 0;
  bool alpha = false;
  
  double nCount = binHi;
  double graphMin = 0;
  double graphMax = nCount;
  //===============================================================
  // get run day info.
  for(int i=grabRunIdx(runStart); i<=grabRunIdx(runStop); i++){
    int runNumberI = grabRunId(i);
    //    //    if (runNumberI >= 20366000 && runNumberI < 21001001) runNumberI = 21001001;
    //    //    if (runNumberI >= 21246012 && runNumberI < 21251016) runNumberI = 21251016;
    TString runNumber(Form("%d",runNumberI));
    TSubString day = runNumber(2,3);
    TSubString run = runNumber(5,3);
    //    cout<<"day: "<<day<<"          run: "<<run<<endl;
    if(runNumberI<21000000) dayI = (runNumberI-20000000)/1000;
    else dayI = (runNumberI-21000000)/1000;
    if(runNumberI==runStart){
      dayStart = day;
    }
    dayEnd = day;
    
    
    if ((dayI != dayPrev) || (runNumberI == runStart)){
      newDay.push_back(run_count);
      dayN.push_back(day);
    }
    
    
    //    cout << "working on run " << runNumber << endl;
    //
    runCount.push_back(run_count);
    runNum.push_back(runNumberI);
    runNum.push_back(runNumberI);
    
    
    run_count++;
    dayPrev = dayI;
  }
  int nRuns = runCount.size();
  //===============================================================
  // draw day line
  vector<TLine*> linesGlobal;
  vector<TText*> textsGlobal;
  
  TText *dayText = new TText(0.01,0.93,"Day");
  dayText->SetNDC();
  dayText->SetX(0.03);
  dayText->SetY(0.97);
  dayText->SetTextSize(0.025);
  TText *globalTitle = new TText(0.93,0.4,"Global Tracks");
  globalTitle->SetNDC();
  globalTitle->SetX(0.96);
  globalTitle->SetY(0.35);
  globalTitle->SetTextAngle(90.);
  globalTitle->SetTextSize(0.05);
  
  // draw all day line
  for (int l=0; l<newDay.size(); l++){
    TLine *lineGlobal = new TLine(newDay.at(l),0,newDay.at(l),0.9); // one for the runs per day
    lineGlobal->SetNDC();
    TText *textGlobal = new TText(newDay.at(l)+2.5,0.35,dayN.at(l)); // one for the day number
    textGlobal->SetNDC();
    
    if (l==0){
      
      lineGlobal->SetX1(0.1);
      lineGlobal->SetX2(0.1);
      lineGlobal->SetY1(0.1);
      lineGlobal->SetY2(0.1);
      linesGlobal.push_back(lineGlobal);
      
      textGlobal->SetText(0.1,0.91,dayStart);
      textGlobal->SetX(0.09);
      textGlobal->SetY(0.93);
      textGlobal->SetTextAngle(90.);
      textGlobal->SetTextSize(0.02);
      textsGlobal.push_back(textGlobal);
      
    }
    else{
      
      lineGlobal->SetX1(0.09+0.81*newDay.at(l)/nCount);
      lineGlobal->SetX2(0.09+0.81*newDay.at(l)/nCount);
      lineGlobal->SetY1(0.1); //0.026
      lineGlobal->SetY2(0.9); //0.947
      lineGlobal->SetLineColor(kRed);
      linesGlobal.push_back(lineGlobal);
      
      textGlobal->SetX(0.09+0.81*newDay.at(l)/nCount);
      textGlobal->SetY(0.93);
      textGlobal->SetTextAngle(90.);
      textGlobal->SetTextSize(0.02);
      textsGlobal.push_back(textGlobal);
      
    }
  }
  
  //===============================================================
  // draw result
  TCanvas *globalCanvas = new TCanvas("globalCanvas","globalCanvas",1680,1200);
  globalCanvas->cd();
  //    gPad->SetTopMargin(0.35);
  //    gPad->SetBottomMargin(0.15);
  TH2D *rangeHist = new TH2D("rangeHist","",1,graphMin,graphMax,1,-2,2); // draw the x axis lable
  rangeHist->SetLineColor(kBlack);
  
  gStyle->SetOptStat(0);
  //  etaRunCountGlobal->SetAxisRange(graphMin,graphMax,"X");
  rangeHist->GetXaxis()->SetTitle("Run Count");
  rangeHist->GetYaxis()->SetTitle("#eta");
  rangeHist->GetYaxis()->CenterTitle();
  rangeHist->GetYaxis()->SetTitleOffset(0.3);
  rangeHist->GetYaxis()->SetTitleSize(0.1);
  rangeHist->GetYaxis()->SetLabelSize(0.08);
  rangeHist->Draw();
  //
  dayText->SetY(0.95);
  globalCanvas->cd(0);
  for (int l=0; l<newDay.size(); l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  
  globalTitle->Draw("same");
  dayText->Draw("same");
  globalCanvas->SaveAs(Form("%s.pdf", outfileName.Data()));
  //===============================================================
  
  return ;
}

//___________________________________________________________________
Int_t grabRunIdx(Int_t runId){
  IntMap::iterator iter = mRunId2RunIdx.find(runId);
  if(iter != mRunId2RunIdx.end())
    return iter->second;
  else
    return -1;
}
//___________________________________________________________________
Int_t grabRunId(Int_t runIdx){
  IntMap::iterator iter = mRunIdx2RunId.find(runIdx);
  if(iter != mRunIdx2RunId.end())
    return iter->second;
  else
    return -1;
}

Bool_t Init(TString runlist="0911_production_7p7GeV_2020_runnumber_DD.dat")
{
  ifstream indata;
  
  indata.open(Form("rootfile/%s", runlist.Data()));
  
  mRunId2RunIdx.clear();
  mRunIdx2RunId.clear();
  if(indata.is_open()){
    cout<<"read in total run number list and recode run number ...";
    Int_t oldId;
    Int_t newId=0;
    while(indata>>oldId){
      mRunId2RunIdx[oldId] = newId;
      mRunIdx2RunId[newId] = oldId;
      newId++;
    }
    cout<<" [OK]"<<endl;
  }else{
    cout<<"Failed to load the total run number list !!!"<<endl;
    return kFALSE;
  }
  indata.close();
  
  return kTRUE;
}
