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

using namespace std;

void primaryRuns_production(TString dir, int runStart, int runEnd, float energy, TString eventConfig, TString outputDir, Bool_t eToF_req){

  bool eToFTrig = false;

  //observables mean yields
  vector<double> meanEtaPrimary;
  vector<double> errorEtaPrimary;
  vector<double> meanPhiPrimary;
  vector<double> errorPhiPrimary;
  vector<double> meanPtPrimary;
  vector<double> errorPtPrimary;
  vector<double> meandEdxPrimary;
  vector<double> errordEdxPrimary;
  
  //  vector<double> meanEToFMatchGlobal;
  //  vector<double> errorEToFMatchGlobal;
  vector<double> meanEToFMatchPrimary;
  vector<double> errorEToFMatchPrimary;

  //  vector<double> meanEToFInvBetaGlobal;
  //  vector<double> errorEToFInvBetaGlobal;
  vector<double> meanEToFInvBetaPrimary;
  vector<double> errorEToFInvBetaPrimary;


  vector<double> meanEToFMassPrimary;
  vector<double> errorEToFMassPrimary;

  //  vector<double> meanToFInvBetaGlobal;
  //  vector<double> errorToFInvBetaGlobal;
  vector<double> meanToFInvBetaPrimary;
  vector<double> errorToFInvBetaPrimary;

  //number of events variable
  vector<int> nEventsVar;
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

  TH1::AddDirectory(kFALSE);
  TH2::AddDirectory(kFALSE);

  //  int dayStartI;
  // int dayEndI;

  TString energyStr(Form("%.1f",energy));
  if (energy == 5.75 ) energyStr = Form("%.2f",energy);
  Ssiz_t first = energyStr.First(".");
  energyStr.Replace(first,1,"p",1);

  TString FXTString="";
  if(eventConfig == "FixedTarget") FXTString = "_FXT";


  TH2D* etaRunCountPrimary = new TH2D("etaRunCount","#eta vs Run Count", 2201, -0.5, 2200.5, 100, -2.0, 2.0);
  

  for(int runNumberI=runStart; runNumberI<=runEnd; runNumberI++){
    if (runNumberI >= 20366000 && runNumberI < 21001001) runNumberI = 21001001;
    if (runNumberI >= 21080030 && runNumberI < 21169001) runNumberI = 21169001;
    TString runNumber(Form("%d",runNumberI));
    TSubString day = runNumber(2,3);
    TSubString run = runNumber(5,3);
    if(runNumberI<21000000) dayI = (runNumberI-20000000)/1000;
    else dayI = (runNumberI-21000000)/1000;
    if(runNumberI==runStart){
      dayStart = day;
    }
    dayEnd = day;

    if (runNumberI>=20352037) eToFTrig=true;


    TString QAFilePath;
    TString QAFileN;
    TString QAFileName;

    QAFilePath = dir;
    QAFileN = Form("AuAu_%sGeV%s_run_%s_runByRunQA_production.root",energyStr.Data(),FXTString.Data(),runNumber.Data());

    /*
    if (eventConfig == "ColliderCenter"){
      QAFilePath = "/star/data03/pwg/bkimel/run20_runbyrun_QA/lfsupc_qa/userfiles/AuAu_" + energyStr + "GeV/day" + day + "/run" + run + "/";
      QAFileN =  "run_" + runNumber + "_runByRunQA.root";
      if (runNumberI >= 21049000) QAFileN = "AuAu_" + energyStr + "GeV_run_" + runNumber + "_runByRunQA.root";
    }
    else if (eventConfig == "FixedTarget" && !eToF_req){
      QAFilePath = "/star/data03/pwg/bkimel/run20_runbyrun_QA/lfsupc_qa/userfiles/AuAu_" + energyStr + "GeV_FXT/day" + day + "/run" + run + "/";
      QAFileN = "run_" + runNumber + "_runByRunQA.root";
      if (runNumberI >= 21049000) QAFileN = "AuAu_" + energyStr + "GeV_FXT_run_" + runNumber + "_runByRunQA.root";
    }
    else if (eventConfig == "FixedTarget" && eToF_req){
      QAFilePath = "/star/data03/pwg/bkimel/run20_runbyrun_QA/lfsupc_qa/userfiles/FXT_eToF_req/AuAu_" + energyStr + "GeV_FXT/day" + day + "/run" + run + "/";
      QAFileN = "run_" + runNumber + "_runByRunQA.root";
      if (runNumberI >= 21049000) QAFileN = "AuAu_" + energyStr + "GeV_FXT_run_" + runNumber + "_runByRunQA.root";
    }*/

    QAFileName = QAFilePath + QAFileN;

    //    TString QAFileName = "/star/u/bkimel/lightflavorspectra/userfiles/AuAu_3p2_Dec6_sched_v1/run_20" + day + run + "_runByRunQA.root";
   //    TString trackQAFileCuts = "/star/data03/pwg/bkimel/run20/repoTest/picodstanalysis/userfiles/" + collision + "/day" + day + "/run" + run + "/trackQA_trackCuts.root";

    if (gSystem->AccessPathName(QAFileName)) continue;

    if ((dayI != dayPrev) || (runNumberI == runStart)){
      newDay.push_back(run_count);
      dayN.push_back(day);
    }

    /*
    if(runNumberI==runStart){
      cout << runNumberI << " == " << runStart << endl;
      newDay.push_back(run_count);                                                                                                                                                                         
      dayN.push_back(day);                                                                                                                                                                                 
    }
    */

    cout << "working on run " << runNumber << endl;

    runCount.push_back(run_count);
    runCount_err.push_back(0);
    runNum.push_back(runNumberI);
    
    TFile *QAFile     = TFile::Open(QAFileName,"READ");
    TH1D* nEventsHist = (TH1D*)QAFile->Get("nEventsHist_cuts");

    TH2D* etaPhiPrimary = (TH2D*)QAFile->Get("etaphiHist_primary");

    TH2D* PtPhiPrimary = (TH2D*)QAFile->Get("pTphiHist_primary");

    TH2D* dEdxPrimary = (TH2D*)QAFile->Get("dEdxHist_primary");

    //    TH1D* eToFMatchGlobal = (TH1D*)QAFile->Get("eTofMatch_global");
    TH1D* eToFMatchPrimary = (TH1D*)QAFile->Get("eTofMatch_primary");

    //    TH2D* eToFInvBetaGlobal = (TH2D*)QAFile->Get("eTofInvBeta_global");
    TH2D* eToFInvBetaPrimary = (TH2D*)QAFile->Get("eTofInvBeta_primary");

    TH2D* eToFMassPrimary = (TH2D*)QAFile->Get("eTofMassVsMom_primary");

    //    TH2D* ToFInvBetaGlobal = (TH2D*)QAFile->Get("tofBeta_global");
    TH2D* ToFInvBetaPrimary = (TH2D*)QAFile->Get("tofBeta_primary");

    TH1D* etaProjPrimary = etaPhiPrimary->ProjectionX();

    (nEventsHist->GetEntries() < nEventsHist->GetBinContent(2)) ? nEventsVar.push_back(nEventsHist->GetBinContent(2)) : nEventsVar.push_back(nEventsHist->GetEntries()); 

    double nEToFEvents;
    if (eToFTrig){
      TH1D* eToFEvents;
      if (runNumberI < 21028011) eToFEvents = (TH1D*)QAFile->Get("eToFEvents");
      else eToFEvents = (TH1D*)QAFile->Get("nEventsHist_eToF");
      nEToFEvents = eToFEvents->GetEntries();
      if (nEToFEvents < eToFEvents->GetBinContent(2)){
       nEToFEvents = eToFEvents->GetBinContent(2);
     }
      //      if (nEToFEvents == 0) nEToFEvents = nEventsVar.at(run_count);
   }
   else if ( !eToFTrig ){
    nEToFEvents = nEventsVar.at(run_count);
  }

  

  meanEtaPrimary.push_back(1.0*etaPhiPrimary->GetMean(1));
  errorEtaPrimary.push_back(1.0*etaPhiPrimary->GetMeanError(1));
  meanPhiPrimary.push_back(1.0*etaPhiPrimary->GetMean(2));
  errorPhiPrimary.push_back(1.0*etaPhiPrimary->GetMeanError(2));
  meanPtPrimary.push_back(1.0*PtPhiPrimary->GetMean(1));
  errorPtPrimary.push_back(1.0*PtPhiPrimary->GetMeanError(1));
  meandEdxPrimary.push_back(1.0*dEdxPrimary->GetMean(2));
  errordEdxPrimary.push_back(1.0*dEdxPrimary->GetMeanError(2));
  if (nEToFEvents != 0){
    meanEToFMatchPrimary.push_back(1.0*eToFMatchPrimary->GetMean(1)*nEventsVar.at(run_count)/nEToFEvents);
    errorEToFMatchPrimary.push_back(1.0*eToFMatchPrimary->GetMeanError(1)*nEventsVar.at(run_count)/nEToFEvents);
  }
  else{
    meanEToFMatchPrimary.push_back(0);
    errorEToFMatchPrimary.push_back(0);
  }
  meanEToFInvBetaPrimary.push_back(1.0*eToFInvBetaPrimary->GetMean(2));
  errorEToFInvBetaPrimary.push_back(1.0*eToFInvBetaPrimary->GetMeanError(2));
  meanEToFMassPrimary.push_back(1.0*eToFMassPrimary->GetMean(2));
  errorEToFMassPrimary.push_back(1.0*eToFMassPrimary->GetMeanError(2));
  meanToFInvBetaPrimary.push_back(1.0*ToFInvBetaPrimary->GetMean(2));
  errorToFInvBetaPrimary.push_back(1.0*ToFInvBetaPrimary->GetMeanError(2));
  
  for (int etaBin=1; etaBin<(etaProjPrimary->GetNbinsX()+1); etaBin++){
    etaRunCountPrimary->Fill(run_count,etaProjPrimary->GetBinCenter(etaBin),1.0*((double)etaProjPrimary->GetBinContent(etaBin))/((double)nEventsVar.at(run_count)));
  }


  run_count++;
  dayPrev = dayI;
  QAFile->Close();
  delete QAFile;
  delete nEventsHist;
  delete etaPhiPrimary;
  delete PtPhiPrimary;
  delete dEdxPrimary;
  delete eToFMatchPrimary;
  delete eToFInvBetaPrimary;
  delete eToFMassPrimary;
  delete ToFInvBetaPrimary;
  delete etaProjPrimary;
  
}
  //  QAFile->Close();

cout << "Finished with all data, making graphs" << endl;

int nRuns = runCount.size();

  //  TString base = "/star/data03/pwg/bkimel/run20_runbyrun_QA/lfsupc_qa/userfiles/AuAu_" + energyStr + "GeV";
TString base = outputDir;

  //  TString base = "/star/u/bkimel/lightflavorspectra/userfiles/AuAu_3p2_Dec6_sched_v1/";
  //  TString base = "/star/data03/pwg/bkimel/run20/repoTest/picodstanalysis/userfiles/" + collision + "/";

if( eventConfig == "ColliderCenter" ) base += "/AuAu_" + energyStr + "GeV_";
else if (eventConfig == "FixedTarget") base += "/AuAu_" + energyStr + "GeV_FXT_";
if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) base += "day" + dayStart;
else base += "days" + dayStart + "_" + dayEnd;
if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) alpha = true;


TFile *outFile = new TFile( base + "_QA_monitoring.root","UPDATE");


  //observable tracking
TGraphErrors *etaPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanEtaPrimary[0]),&(runCount_err[0]),&(errorEtaPrimary[0]));
etaPrimaryGr->SetMarkerStyle(20);
etaPrimaryGr->SetTitle("");
etaPrimaryGr->GetXaxis()->SetTitle("RunCount");
etaPrimaryGr->GetYaxis()->SetTitle("Mean #eta");
etaPrimaryGr->GetYaxis()->CenterTitle();
etaPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
etaPrimaryGr->GetYaxis()->SetTitleSize(0.1);
etaPrimaryGr->GetYaxis()->SetLabelSize(0.08);

TGraphErrors *phiPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanPhiPrimary[0]),&(runCount_err[0]),&(errorPhiPrimary[0]));
phiPrimaryGr->SetMarkerStyle(20);
phiPrimaryGr->SetTitle("");
phiPrimaryGr->GetXaxis()->SetTitle("RunCount");
phiPrimaryGr->GetYaxis()->SetTitle("Mean #phi");
phiPrimaryGr->GetYaxis()->CenterTitle();
phiPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
phiPrimaryGr->GetYaxis()->SetTitleSize(0.1);
phiPrimaryGr->GetYaxis()->SetLabelSize(0.08);

TGraphErrors *pTPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanPtPrimary[0]),&(runCount_err[0]),&(errorPtPrimary[0]));
pTPrimaryGr->SetMarkerStyle(20);
pTPrimaryGr->SetTitle("");
pTPrimaryGr->GetXaxis()->SetTitle("RunCount");
pTPrimaryGr->GetYaxis()->SetTitle("Mean p_{T}");
pTPrimaryGr->GetYaxis()->CenterTitle();
pTPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
pTPrimaryGr->GetYaxis()->SetTitleSize(0.1);
pTPrimaryGr->GetYaxis()->SetLabelSize(0.08);

TGraphErrors *dEdxPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meandEdxPrimary[0]),&(runCount_err[0]),&(errordEdxPrimary[0]));
dEdxPrimaryGr->SetMarkerStyle(20);
dEdxPrimaryGr->SetTitle("");
dEdxPrimaryGr->GetXaxis()->SetTitle("RunCount");
dEdxPrimaryGr->GetYaxis()->SetTitle("Mean dE/dx");
dEdxPrimaryGr->GetYaxis()->CenterTitle();
dEdxPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
dEdxPrimaryGr->GetYaxis()->SetTitleSize(0.1);
dEdxPrimaryGr->GetYaxis()->SetLabelSize(0.08);


TGraphErrors *eToFMatchPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanEToFMatchPrimary[0]),&(runCount_err[0]),&(errorEToFMatchPrimary[0]));
eToFMatchPrimaryGr->SetMarkerStyle(20);
eToFMatchPrimaryGr->SetTitle("");
eToFMatchPrimaryGr->GetXaxis()->SetTitle("RunCount");
eToFMatchPrimaryGr->GetYaxis()->SetTitle("Mean eToF N_{match}");
eToFMatchPrimaryGr->GetYaxis()->CenterTitle();
eToFMatchPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
eToFMatchPrimaryGr->GetYaxis()->SetTitleSize(0.1);
eToFMatchPrimaryGr->GetYaxis()->SetLabelSize(0.08);


TGraphErrors *eToFInvBetaPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanEToFInvBetaPrimary[0]),&(runCount_err[0]),&(errorEToFInvBetaPrimary[0]));
eToFInvBetaPrimaryGr->SetMarkerStyle(20);
eToFInvBetaPrimaryGr->SetTitle("");
eToFInvBetaPrimaryGr->GetXaxis()->SetTitle("RunCount");
eToFInvBetaPrimaryGr->GetYaxis()->SetTitle("Mean eToF 1/#beta");
eToFInvBetaPrimaryGr->GetYaxis()->CenterTitle();
eToFInvBetaPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
eToFInvBetaPrimaryGr->GetYaxis()->SetTitleSize(0.1);
eToFInvBetaPrimaryGr->GetYaxis()->SetLabelSize(0.08);


TGraphErrors *eToFMassPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanEToFMassPrimary[0]),&(runCount_err[0]),&(errorEToFMassPrimary[0]));
eToFMassPrimaryGr->SetMarkerStyle(20);
eToFMassPrimaryGr->SetTitle("");
eToFMassPrimaryGr->GetXaxis()->SetTitle("RunCount");
eToFMassPrimaryGr->GetYaxis()->SetTitle("Mean eToF Mass");
eToFMassPrimaryGr->GetYaxis()->CenterTitle();
eToFMassPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
eToFMassPrimaryGr->GetYaxis()->SetTitleSize(0.1);
eToFMassPrimaryGr->GetYaxis()->SetLabelSize(0.08);


TGraphErrors *ToFInvBetaPrimaryGr = new TGraphErrors(nRuns,&(runCount[0]),&(meanToFInvBetaPrimary[0]),&(runCount_err[0]),&(errorToFInvBetaPrimary[0]));
ToFInvBetaPrimaryGr->SetMarkerStyle(20);
ToFInvBetaPrimaryGr->SetTitle("");
ToFInvBetaPrimaryGr->GetXaxis()->SetTitle("RunCount");
ToFInvBetaPrimaryGr->GetYaxis()->SetTitle("Mean bToF 1/#beta");
ToFInvBetaPrimaryGr->GetYaxis()->CenterTitle();
ToFInvBetaPrimaryGr->GetYaxis()->SetTitleOffset(0.45);
ToFInvBetaPrimaryGr->GetYaxis()->SetTitleSize(0.1);
ToFInvBetaPrimaryGr->GetYaxis()->SetLabelSize(0.08);


etaPrimaryGr->Write("etaPrimary");
phiPrimaryGr->Write("phiPrimary");
pTPrimaryGr->Write("pTPrimary");
dEdxPrimaryGr->Write("dEdxPrimary");
eToFMatchPrimaryGr->Write("eToFMatchPrimary");
eToFInvBetaPrimaryGr->Write("eToFInvBetaPrimary");
ToFInvBetaPrimaryGr->Write("ToFInvBetaPrimary");

double nCount = etaPrimaryGr->GetXaxis()->GetXmax() - etaPrimaryGr->GetXaxis()->GetXmin();


gStyle->SetOptDate(kFALSE);

  //  TCanvas *eventCanvas = new TCanvas("eventCanvas","eventCanvas",2940,2100);
TCanvas *primaryCanvas = new TCanvas("primaryCanvas","primaryCanvas",1680,1200);
TCanvas *primaryCanvasPID = new TCanvas("primaryCanvasPID","primaryCanvasPID",1680,1200);

vector<TLine*> linesGlobal;
vector<TText*> textsGlobal;

vector<TLine*> linesGlobalPID;
vector<TText*> textsGlobalPID;

TText *dayText = new TText(0.01,0.95,"Day");
dayText->SetNDC();
dayText->SetX(0.03);
dayText->SetY(0.97);
dayText->SetTextSize(0.025);
TText *primaryTitle = new TText(0.93,0.4,"Primary Tracks");
primaryTitle->SetNDC();
primaryTitle->SetX(0.96);
primaryTitle->SetY(0.35);
primaryTitle->SetTextAngle(90.);
primaryTitle->SetTextSize(0.05);

  //  double nRuns = runCount.at(((int)(runCount.size())-1)-((int)(runCount.at(0))));

  /*  if(dayEnd != dayStart){
    TText *text0 = new TText(newDay.at(0)/4,0.35,dayStart); 
    text0->SetNDC();
    text0->SetX(0.1);
    text0->SetY(0.91);
    text0->SetTextAngle(60.);
    text0->SetTextSize(0.04);
    texts.push_back(text0);

    }*/

for (int l=0; l<newDay.size(); l++){


  TLine *lineGlobal = new TLine(newDay.at(l),0,newDay.at(l),0.9);
  lineGlobal->SetNDC();
  TText *textGlobal = new TText(newDay.at(l)+2.5,0.35,dayN.at(l));
  textGlobal->SetNDC();

  TLine *lineGlobalPID = new TLine(newDay.at(l),0,newDay.at(l),0.9);
  lineGlobalPID->SetNDC();
  TText *textGlobalPID = new TText(newDay.at(l)+2.5,0.35,dayN.at(l));
  textGlobalPID->SetNDC();


  if (l==0){

    lineGlobal->SetX1(0.1);
    lineGlobal->SetX2(0.1);
    lineGlobal->SetY1(0.1);
    lineGlobal->SetY2(0.1);
    linesGlobal.push_back(lineGlobal);

    textGlobal->SetText(0.1,0.91,dayStart);
    textGlobal->SetX(0.09);
    textGlobal->SetY(.96);
    textGlobal->SetTextAngle(90.);
      //      text->SetTextSize(0.02);
    textGlobal->SetTextSize(0.02);
    textsGlobal.push_back(textGlobal);

    lineGlobalPID->SetX1(0.1);
    lineGlobalPID->SetX2(0.1);
    lineGlobalPID->SetY1(0.1);
    lineGlobalPID->SetY2(0.1);
    linesGlobalPID.push_back(lineGlobalPID);

    textGlobalPID->SetText(0.1,0.91,dayStart);
    textGlobalPID->SetX(0.09);
    textGlobalPID->SetY(.94);
    textGlobalPID->SetTextAngle(90.);
      //      text->SetTextSize(0.02);
    textGlobalPID->SetTextSize(0.02);
    textsGlobalPID.push_back(textGlobalPID);
  }
  else{

    lineGlobal->SetX1(0.09+0.81*newDay.at(l)/nCount);
    lineGlobal->SetX2(0.09+0.81*newDay.at(l)/nCount);
    lineGlobal->SetY1(0.026);
    lineGlobal->SetY2(0.947);
    lineGlobal->SetLineColor(kRed);
    linesGlobal.push_back(lineGlobal);

    textGlobal->SetX(0.09+0.81*newDay.at(l)/nCount);
    textGlobal->SetY(0.96);
    textGlobal->SetTextAngle(90.);
    textGlobal->SetTextSize(0.02);
    textsGlobal.push_back(textGlobal);

    lineGlobalPID->SetX1(0.09+0.81*newDay.at(l)/nCount);
    lineGlobalPID->SetX2(0.09+0.81*newDay.at(l)/nCount);
    lineGlobalPID->SetY1(0.033);
    lineGlobalPID->SetY2(0.934);
    lineGlobalPID->SetLineColor(kRed);
    linesGlobalPID.push_back(lineGlobalPID);

    textGlobalPID->SetX(0.09+0.81*newDay.at(l)/nCount);
    textGlobalPID->SetY(0.94);
    textGlobalPID->SetTextAngle(90.);
    textGlobalPID->SetTextSize(0.02);
    textsGlobalPID.push_back(textGlobalPID);
  }
}


primaryCanvas->cd();
gPad->SetTopMargin(0.35);
gPad->SetBottomMargin(0.15);
primaryCanvas->Divide(1,4,0.01,0);

primaryCanvas->cd(1);
etaPrimaryGr->Draw("AP");

primaryCanvas->cd(2);
phiPrimaryGr->Draw("AP");

primaryCanvas->cd(3);
pTPrimaryGr->Draw("AP");

  //  double graphMax = etaGlobalGr->GetXaxis()->GetXmax()+1+etaGlobalGr->GetXaxis()->GetBinWidth(2);
double graphMax = etaPrimaryGr->GetXaxis()->GetXmax();
  //  double graphMin = etaGlobalGr->GetXaxis()->GetXmin()+1+etaGlobalGr->GetXaxis()->GetBinWidth(2);
double graphMin = etaPrimaryGr->GetXaxis()->GetXmin();

TH2D *rangeHist = new TH2D("rangeHist","",1,graphMin,graphMax,1,-2,2);
rangeHist->SetLineColor(kBlack);

primaryCanvas->cd(4);
gStyle->SetOptStat(0);
  //  etaRunCountGlobal->SetAxisRange(graphMin,graphMax,"X");
rangeHist->GetXaxis()->SetTitle("Run Count");
rangeHist->GetYaxis()->SetTitle("#eta");
rangeHist->GetYaxis()->CenterTitle();
rangeHist->GetYaxis()->SetTitleOffset(0.3);
rangeHist->GetYaxis()->SetTitleSize(0.1);
rangeHist->GetYaxis()->SetLabelSize(0.08);
rangeHist->Draw();

primaryCanvas->cd(4);
gStyle->SetOptStat(0);
rangeHist->Draw();
etaRunCountPrimary->GetZaxis()->SetLabelSize(0.08);  
etaRunCountPrimary->Draw("SAME COLZ");


  //  dayText->SetY(0.93);
dayText->SetY(0.95);

primaryCanvas->cd(0);
for (int l=0; l<newDay.size(); l++){
  linesGlobalPID.at(l)->Draw("same");
    textsGlobalPID.at(l)->Draw("same"); //l+1?
  }
  primaryTitle->Draw("same");
  dayText->Draw("same");
  primaryCanvas->SaveAs( base + "_PrimaryTrack.png");


  dayText->SetY(0.95);

  primaryCanvasPID->cd();
  gPad->SetTopMargin(0.35);
  gPad->SetBottomMargin(0.15);
  primaryCanvasPID->Divide(1,5,0.01,0);


  primaryCanvasPID->cd(1);
  dEdxPrimaryGr->Draw("AP");

  primaryCanvasPID->cd(2);
  eToFMatchPrimaryGr->Draw("AP");

  primaryCanvasPID->cd(3);
  eToFInvBetaPrimaryGr->Draw("AP");

  primaryCanvasPID->cd(4);
  eToFMassPrimaryGr->Draw("AP");

  primaryCanvasPID->cd(5);
  ToFInvBetaPrimaryGr->Draw("AP");


  dayText->SetY(0.95);

  primaryCanvasPID->cd(0);
  for (int l=0; l<newDay.size(); l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  primaryTitle->Draw("same");
  dayText->Draw("same");
  primaryCanvasPID->SaveAs( base + "_PrimaryTrackPID.png");



  outFile->Close();

}

