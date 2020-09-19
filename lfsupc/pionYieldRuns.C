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

void pionYieldRuns(int runStart, int runEnd, float energy, TString eventConfig, TString outputDir, Bool_t eToF_req){

  bool eToFTrig = false;

  //observables mean yields
  //pion yield variables
  vector<double> runNum;
  vector<double> Amp_pi;
  vector<double> Amp_piY;
  vector<double> Amp_K;
  vector<double> Amp_pro;
  vector<double> piK_rat;

  vector<double> piW;
  vector<double> piW_err;

  vector<double> piM;
  vector<double> piM_err;

  vector<double> Amp_pi_err;
  vector<double> Amp_piY_err;
  vector<double> Amp_K_err;
  vector<double> Amp_pro_err;
  vector<double> piK_rat_err;


  vector<double> Amp_pi_tof;
  vector<double> Amp_pi_tof_err;

  vector<double> Amp_pi_etof;
  vector<double> Amp_pi_etof_err;


  //number of events variable
  vector<int> nEventsVar;

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

  TH1D *piYield_hist = new TH1D("piYield_hist","#pi Yield Histogram",100,0.,5.);
  TH1D *kYield_hist = new TH1D("kYield_hist","k Yield Histogram",100,0.,5.);
  TH1D *proYield_hist = new TH1D("proYield_hist","p Yield Histogram",100,0.,5.);

  TCanvas *cdEdx = new TCanvas("cdEdx","cdEdx",1200,600);
  cdEdx->Divide(2,1);

  TString energyStr(Form("%.1f",energy));
  if (energy == 5.75 ) energyStr = Form("%.2f",energy);
  Ssiz_t first = energyStr.First(".");
  energyStr.Replace(first,1,"p",1);


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
    }

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

    TH1D* dEdxSliceEta = (TH1D*)QAFile->Get("dEdxSliceEta_primary");
    TH1D* dEdxSliceY = (TH1D*)QAFile->Get("dEdxSliceY_primary");
    TH1D* etofSlice = (TH1D*)QAFile->Get("eTofSlice_primary");
    TH1D* tofSlice = (TH1D*)QAFile->Get("tofSlice_primary");

    (nEventsHist->GetEntries() < nEventsHist->GetBinContent(2)) ? nEventsVar.push_back(nEventsHist->GetBinContent(2)) : nEventsVar.push_back(nEventsHist->GetEntries()); 

    double nEToFEvents;
    TH1D* eToFEvents;
    if (eToFTrig){
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

    if (tofSlice->GetMaximum() == 0){
      Amp_pi_tof.push_back(0);
      Amp_pi_tof_err.push_back(0);
    }
    else{
      Amp_pi_tof.push_back(1.0*tofSlice->GetEntries()/nEventsVar.at(run_count));
      Amp_pi_tof_err.push_back(1.0*sqrt(tofSlice->GetEntries())/nEventsVar.at(run_count));     
    }


    if (etofSlice->GetMaximum() == 0){
      Amp_pi_etof.push_back(0);
      Amp_pi_etof_err.push_back(0);
    }
    else{
      if (nEToFEvents != 0){
	Amp_pi_etof.push_back(1.0*etofSlice->GetEntries()/nEToFEvents);
	Amp_pi_etof_err.push_back(1.0*sqrt(etofSlice->GetEntries())/nEToFEvents);
      }
      else{
	Amp_pi_etof.push_back(0);
	Amp_pi_etof_err.push_back(0);
      }
    }


    if (dEdxSliceEta->GetMaximum() == 0){
      
      Amp_pi.push_back(0);
      Amp_K.push_back(0);
      Amp_pro.push_back(0);
      piK_rat.push_back(0);
      
      piW.push_back(0);
      piW_err.push_back(0);
      
      piM.push_back(0);
      piM_err.push_back(0);
      
      Amp_pi_err.push_back(0);
      Amp_K_err.push_back(0);
      Amp_pro_err.push_back(0);
      piK_rat_err.push_back(0);
    }
    else{
      cdEdx->cd(1);
      dEdxSliceEta->Draw();
      TF1 *f1 = new TF1("f1","gaus(0)+gaus(3)+gaus(6)");
      f1->SetParameters(40,0.43,0.038,5,0.64,0.05,5,1,0.038);
      f1->SetParNames("Amp._{#pi}","#mu_{#pi}","#sigma_{#pi}","Amp._{K}","#mu_{K}","#sigma_{K}","Amp._{p}","#mu_{p}","#sigma_{p}");
      f1->SetParLimits(0,0,200000);
      f1->SetParLimits(1,0.4,0.46);
      f1->SetParLimits(2,0.02,0.1);
      f1->SetParLimits(3,0,20000);
      f1->SetParLimits(4,0.6,0.68);
      f1->SetParLimits(5,0.02,0.1);
      f1->SetParLimits(6,0,20000);
      f1->SetParLimits(7,0.97,1.3);
      f1->SetParLimits(8,0.02,0.1);
      dEdxSliceEta->Fit("f1");
      gStyle->SetOptFit(1);
      f1->Draw("same");
      
      double norm = 1.0/(nEventsVar.at(run_count)) * 1.0/(2*TMath::Pi()) * 1.0/(0.475) * 1.0/0.05 * 1.0/0.5 * 1.0/dEdxSliceEta->GetBinWidth(2);
      
      Amp_pi.push_back((f1->GetParameter(0))*(f1->GetParameter(2))*sqrt(2*TMath::Pi())*norm);
      Amp_K.push_back((f1->GetParameter(3))*(f1->GetParameter(5))*sqrt(2*TMath::Pi())*norm);
      Amp_pro.push_back((f1->GetParameter(6))**(f1->GetParameter(8))*sqrt(2*TMath::Pi())*norm);
      piK_rat.push_back(Amp_pi.at(run_count)/Amp_K.at(run_count));
      
      piW.push_back(f1->GetParameter(2));
      piW_err.push_back(f1->GetParError(2));
      piM.push_back(f1->GetParameter(1));
      piM_err.push_back(f1->GetParError(1));
      
      if(Amp_pi.at(run_count) != 0) piYield_hist->Fill(Amp_pi.at(run_count));
      if(Amp_K.at(run_count) != 0) kYield_hist->Fill(Amp_K.at(run_count));
      if(Amp_pro.at(run_count) != 0) proYield_hist->Fill(Amp_pro.at(run_count));
      
      Amp_pi_err.push_back(Amp_pi.at(run_count)*sqrt(pow(f1->GetParError(0)/f1->GetParameter(0),2)+pow(f1->GetParError(2)/f1->GetParameter(2),2)));
      Amp_K_err.push_back(Amp_K.at(run_count)*sqrt(pow(f1->GetParError(3)/f1->GetParameter(3),2)+pow(f1->GetParError(5)/f1->GetParameter(5),2)));
      Amp_pro_err.push_back(Amp_pro.at(run_count)*sqrt(pow(f1->GetParError(6)/f1->GetParameter(6),2)+pow(f1->GetParError(8)/f1->GetParameter(8),2)));
      if (Amp_pi.at(run_count) == 0 || Amp_pi_err.at(run_count) == 0 || Amp_K.at(run_count) == 0 || Amp_K_err.at(run_count) == 0){
	piK_rat_err.push_back(0);
      }
      else{
	piK_rat_err.push_back(piK_rat.at(run_count)*sqrt(pow(Amp_pi_err.at(run_count)/Amp_pi.at(run_count),2)+pow(Amp_K_err.at(run_count)/Amp_K.at(run_count),2)));
      }	
    }      
    
    if (dEdxSliceY->GetMaximum() == 0){
      Amp_piY.push_back(0);
      Amp_piY_err.push_back(0);
    }
    else{
      cdEdx->cd(2);
      dEdxSliceY->Draw();
      TF1 *f2 = new TF1("f2","gaus(0)+gaus(3)+gaus(6)");
      f2->SetParameters(40,0.43,0.038,5,0.64,0.05,5,1,0.038);
      f2->SetParNames("Amp._{#pi}","#mu_{#pi}","#sigma_{#pi}","Amp._{K}","#mu_{K}","#sigma_{K}","Amp._{p}","#mu_{p}","#sigma_{p}");
      f2->SetParLimits(0,0,200000);
      f2->SetParLimits(1,0.4,0.46);
      f2->SetParLimits(2,0.02,0.1);
      f2->SetParLimits(3,0,20000);
      f2->SetParLimits(4,0.6,0.68);
      f2->SetParLimits(5,0.02,0.1);
      f2->SetParLimits(6,0,20000);
      f2->SetParLimits(7,0.97,1.3);
      f2->SetParLimits(8,0.02,0.1);
      dEdxSliceY->Fit("f2");
      gStyle->SetOptFit(1);
      f2->Draw("same");
      double normY = 1.0/(nEventsVar.at(run_count)) * 1.0/(2*TMath::Pi()) * 1.0/(0.475) * 1.0/0.05 * 1.0/0.5 * 1.0/dEdxSliceY->GetBinWidth(2);
      
      
      Amp_piY.push_back((f2->GetParameter(0))*(f2->GetParameter(2))*sqrt(2*TMath::Pi())*normY);
      Amp_piY_err.push_back(Amp_piY.at(run_count)*sqrt(pow(f2->GetParError(0)/f2->GetParameter(0),2)+pow(f2->GetParError(2)/f2->GetParameter(2),2)));
    }
    

    run_count++;
    dayPrev = dayI;
    QAFile->Close();
    delete QAFile;
    delete nEventsHist;
    delete dEdxSliceEta;
    delete dEdxSliceY;
    delete etofSlice;
    delete tofSlice;
    delete eToFEvents;
   
  }
  //  QAFile->Close();

  cout << "Finished with all data, making graphs" << endl;
  
  int nRuns = runCount.size();

  //  TString base = "/star/data03/pwg/bkimel/run20_runbyrun_QA/lfsupc_qa/userfiles/AuAu_" + energyStr + "GeV";
  TString base = outputDir;

  //  TString base = "/star/u/bkimel/lightflavorspectra/userfiles/AuAu_3p2_Dec6_sched_v1/";
  //  TString base = "/star/data03/pwg/bkimel/run20/repoTest/picodstanalysis/userfiles/" + collision + "/";
  
  if (runEnd >= 21049000){
    if( eventConfig == "ColliderCenter" ) base += "/AuAu_" + energyStr + "GeV_";
    else if (eventConfig == "FixedTarget") base += "/AuAu_" + energyStr + "GeV_FXT_";
    if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) base += "day" + dayStart;
    else base += "days" + dayStart + "_" + dayEnd;
    if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) alpha = true;
  }
  else{
    if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) base += "/day" + dayStart;
    else base += "days" + dayStart + "_" + dayEnd;
    if ( (int)(runNum.at(0)/1000) == (int)(runNum.at(nRuns-1)/1000) ) alpha = true;
  }


  TFile *outFile = new TFile( base + "_QA_monitoring.root","UPDATE");


  //pion yields tracking
  TGraphErrors *yieldsPiY = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_piY[0]),&(runCount_err[0]),&(Amp_piY_err[0]));
  TGraphErrors *yieldsPi = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_pi[0]),&(runCount_err[0]),&(Amp_pi_err[0]));
  TGraphErrors *yieldsK = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_K[0]),&(runCount_err[0]),&(Amp_K_err[0]));
  TGraphErrors *yieldsPro = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_pro[0]),&(runCount_err[0]),&(Amp_pro_err[0]));
  TGraphErrors *piK_ratio = new TGraphErrors(nRuns,&(runCount[0]),&(piK_rat[0]),&(runCount_err[0]),&(piK_rat_err[0]));

  TGraphErrors *PiW = new TGraphErrors(nRuns,&(runCount[0]),&(piW[0]),&(runCount_err[0]),&(piW_err[0]));
  TGraphErrors *PiM = new TGraphErrors(nRuns,&(runCount[0]),&(piM[0]),&(runCount_err[0]),&(piM_err[0]));

  TGraphErrors *yieldsPi_run = new TGraphErrors(nRuns,&(runNum[0]),&(Amp_pi[0]),&(runCount_err[0]),&(Amp_pi_err[0]));
  TGraphErrors *yieldsK_run = new TGraphErrors(nRuns,&(runNum[0]),&(Amp_K[0]),&(runCount_err[0]),&(Amp_K_err[0]));
  TGraphErrors *yieldsPro_run = new TGraphErrors(nRuns,&(runNum[0]),&(Amp_pro[0]),&(runCount_err[0]),&(Amp_pro_err[0]));
  TGraphErrors *piK_ratio_run = new TGraphErrors(nRuns,&(runNum[0]),&(piK_rat[0]),&(runCount_err[0]),&(piK_rat_err[0]));


  TGraphErrors *yieldsPiToF = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_pi_tof[0]),&(runCount_err[0]),&(Amp_pi_tof_err[0]));
  yieldsPiToF->SetMarkerStyle(20);
  yieldsPiToF->SetTitle("");
  yieldsPiToF->GetXaxis()->SetTitle("RunCount");
  yieldsPiToF->GetYaxis()->SetTitle("bToF #pi Yield");
  yieldsPiToF->GetYaxis()->CenterTitle();
  yieldsPiToF->GetYaxis()->SetTitleOffset(0.45);
  yieldsPiToF->GetYaxis()->SetTitleSize(0.1);
  yieldsPiToF->GetYaxis()->SetLabelSize(0.08);

  TGraphErrors *yieldsPiEToF = new TGraphErrors(nRuns,&(runCount[0]),&(Amp_pi_etof[0]),&(runCount_err[0]),&(Amp_pi_etof_err[0]));
  yieldsPiEToF->SetMarkerStyle(20);
  yieldsPiEToF->SetTitle("");
  yieldsPiEToF->GetXaxis()->SetTitle("RunCount");
  yieldsPiEToF->GetYaxis()->SetTitle("eToF #pi Yield");
  yieldsPiEToF->GetYaxis()->CenterTitle();
  yieldsPiEToF->GetYaxis()->SetTitleOffset(0.45);
  yieldsPiEToF->GetYaxis()->SetTitleSize(0.1);
  yieldsPiEToF->GetYaxis()->SetLabelSize(0.08);

  

  //funCount
  yieldsPiY->SetMarkerStyle(20);
  //  yieldsPi->SetTitle("Uncorrected #pi Yield vs RunCount");
  yieldsPiY->SetTitle("");
  yieldsPiY->GetXaxis()->SetTitle("RunCount");
  //  yieldsPi->GetXaxis()->LabelsOption("v");
  yieldsPiY->GetYaxis()->SetTitle("Uncorrected #pi Yield y_{CM}");
  yieldsPiY->GetYaxis()->SetTitleOffset(1.0);
  
  
  //funCount
  yieldsPi->SetMarkerStyle(20);
  //  yieldsPi->SetTitle("Uncorrected #pi Yield vs RunCount");
  yieldsPi->SetTitle("");
  yieldsPi->GetXaxis()->SetTitle("RunCount");
  //  yieldsPi->GetXaxis()->LabelsOption("v");
  yieldsPi->GetYaxis()->SetTitle("Uncorrected #pi Yield");
  yieldsPi->GetYaxis()->SetTitleOffset(1.0);


  PiW->SetMarkerStyle(20);
  //  yieldsPi->SetTitle("Uncorrected #pi Yield vs RunCount");
  PiW->SetTitle("");
  PiW->GetXaxis()->SetTitle("RunCount");
  //  yieldsPi->GetXaxis()->LabelsOption("v");
  PiW->GetYaxis()->SetTitle("#pi Width");
  PiW->GetYaxis()->SetTitleOffset(1.3);

  PiM->SetMarkerStyle(20);
  //  yieldsPi->SetTitle("Uncorrected #pi Yield vs RunCount");
  PiM->SetTitle("");
  PiM->GetXaxis()->SetTitle("RunCount");
  //  yieldsPi->GetXaxis()->LabelsOption("v");
  PiM->GetYaxis()->SetTitle("#pi Centroid");
  PiM->GetYaxis()->SetTitleOffset(1.3);

  yieldsPro->SetMarkerStyle(20);
  // yieldsPro->SetTitle("Uncorrected p Yield vs RunCount");
  yieldsPro->SetTitle("");
  yieldsPro->GetXaxis()->SetTitle("RunCount");
  // yieldsPro->GetXaxis()->LabelsOption("v");
  yieldsPro->GetYaxis()->SetTitle("uncorrected p Yield");
  yieldsPro->GetYaxis()->SetTitleOffset(1.3);

  yieldsK->SetMarkerStyle(20);
  //  yieldsK->SetTitle("Uncorrected K Yield vs RunCount");
  yieldsK->SetTitle("");
  yieldsK->GetXaxis()->SetTitle("RunCount");
  // yieldsK->GetXaxis()->LabelsOption("v");
  yieldsK->GetYaxis()->SetTitle("Uncorrected K Yield");
  yieldsK->GetYaxis()->SetTitleOffset(1.3);

  piK_ratio->SetMarkerStyle(20);
  //  piK_ratio->SetTitle("#pi/K vs RunCount");
  piK_ratio->SetTitle("");
  piK_ratio->GetXaxis()->SetTitle("RunCount");
  // yieldsK->GetXaxis()->LabelsOption("v");
  piK_ratio->GetYaxis()->SetTitle("#pi/K");
  piK_ratio->GetYaxis()->SetTitleOffset(1.3);

  //runNum
  yieldsPi_run->SetMarkerStyle(20);
  //  yieldsPi_run->SetTitle("Uncorrected #pi Yield vs RunCount");
  yieldsPi_run->SetTitle("");
  yieldsPi_run->GetXaxis()->SetTitle("Runnumber");
  //  yieldsPi_run->GetXaxis()->LabelsOption("v");
  yieldsPi_run->GetYaxis()->SetTitle("Uncorrected #pi Yield");
  yieldsPi_run->GetYaxis()->SetTitleOffset(1.3);

  yieldsPro_run->SetMarkerStyle(20);
  // yieldsPro_run->SetTitle("Uncorrected p Yield vs RunCount");
  yieldsPro_run->SetTitle("");
  yieldsPro_run->GetXaxis()->SetTitle("Runnumber");
  // yieldsPro_run->GetXaxis()->LabelsOption("v");
  yieldsPro_run->GetYaxis()->SetTitle("uncorrected p Yield");
  yieldsPro_run->GetYaxis()->SetTitleOffset(1.3);

  yieldsK_run->SetMarkerStyle(20);
  //  yieldsK_run->SetTitle("Uncorrected K Yield vs RunCount");
  yieldsK_run->SetTitle("");
  yieldsK_run->GetXaxis()->SetTitle("Runnumber");
  // yieldsK_run->GetXaxis()->LabelsOption("v");
  yieldsK_run->GetYaxis()->SetTitle("Uncorrected K Yield");
  yieldsK_run->GetYaxis()->SetTitleOffset(1.3);

  piK_ratio_run->SetMarkerStyle(20);
  //  piK_ratio_run->SetTitle("#pi/K vs RunCount");
  piK_ratio_run->SetTitle("");
  piK_ratio_run->GetXaxis()->SetTitle("Runnumber");
  // yieldsK->GetXaxis()->LabelsOption("v");
  piK_ratio_run->GetYaxis()->SetTitle("#pi/K");
  piK_ratio_run->GetYaxis()->SetTitleOffset(1.3);


  yieldsPiEToF->Write("yieldsPiEToF");
  yieldsPiToF->Write("yieldsPiToF");
  yieldsPiY->Write("yieldsPi");
  yieldsPi->Write("yieldsPi");
  yieldsK->Write("yieldsK");
  yieldsPro->Write("yieldsPro");
  piK_ratio->Write("piK_ratio");
  yieldsPi_run->Write("yieldsPi_run");
  yieldsK_run->Write("yieldsK_run");
  yieldsPro_run->Write("yieldsPro_run");
  piK_ratio_run->Write("piK_ratio_run");
  PiW->Write("PiW");
  PiM->Write("PiM");

  double nCount = yieldsPi->GetXaxis()->GetXmax() - yieldsPi->GetXaxis()->GetXmin();


  gStyle->SetOptDate(kFALSE);
  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  //  TCanvas *eventCanvas = new TCanvas("eventCanvas","eventCanvas",2940,2100);
  TCanvas *primaryCanvasYields = new TCanvas("primaryCanvasYields","primaryCanvasYields",1680,1200);

  vector<TLine*> linesGlobalPID;
  vector<TText*> textsGlobalPID;

  vector<TLine*> linesSingle;
  vector<TText*> textsSingle;


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

    TLine *lineGlobalPID = new TLine(newDay.at(l),0,newDay.at(l),0.9);
    lineGlobalPID->SetNDC();
    TText *textGlobalPID = new TText(newDay.at(l)+2.5,0.35,dayN.at(l));
    textGlobalPID->SetNDC();


    TLine *lineSingle = new TLine(newDay.at(l),0,newDay.at(l),0.9);
    lineSingle->SetNDC();
    TText *textSingle = new TText(newDay.at(l)+2.5,0.35,dayN.at(l));
    textSingle->SetNDC();


    if (l==0){

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

      lineSingle->SetX1(0.1);
      lineSingle->SetX2(0.1);
      lineSingle->SetY1(0.1);
      lineSingle->SetY2(0.1);
      linesSingle.push_back(lineSingle);

      textSingle->SetText(0.1,0.91,dayStart);
      textSingle->SetX(0.1);
      textSingle->SetY(0.91);
      textSingle->SetTextAngle(90.);
      textSingle->SetTextSize(0.02);
      textsSingle.push_back(textSingle);

    }
    else{

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


      lineSingle->SetX1(0.1+0.8*newDay.at(l)/nCount);
      lineSingle->SetX2(0.1+0.8*newDay.at(l)/nCount);
      lineSingle->SetY1(0.1);
      lineSingle->SetY2(0.9);
      lineSingle->SetLineColor(kRed);
      linesSingle.push_back(lineSingle);

      textSingle->SetText(0.1,0.91,dayN.at(l));
      textSingle->SetX(0.1+0.8*newDay.at(l)/nCount);
      textSingle->SetY(0.91);
      textSingle->SetTextAngle(90.);
      textSingle->SetTextSize(0.02);
      textsSingle.push_back(textSingle);

    }
  }


  c1->cd();

  c1->Clear();
  yieldsPiY->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same"); //l+1?
  }
  c1->SaveAs( base + "_PionYields_yCM.png");


  c1->Clear();
  yieldsPi->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same"); //l+1?
  }
  c1->SaveAs( base + "_PionYields.png");

  c1->Clear();
  PiW->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same"); //l+1?
  }
  c1->SaveAs( base + "_PionWidths.png");

  c1->Clear();
  PiM->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same"); //l+1?
  }
  c1->SaveAs( base + "_PionCentroidss.png");

  c1->Clear();
  yieldsPro->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_ProtonYields.png");
  
  c1->Clear();
  yieldsK->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_KaonYields.png");

  c1->Clear();
  piK_ratio->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_pion_Kaon_ratio.png");



  c1->cd();
  c1->SetLogy();
  yieldsPi->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_PionYields_log.png");

  c1->Clear();
  c1->SetLogy();
  yieldsPro->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_ProtonYields_log.png");
  
  c1->Clear();
  c1->SetLogy();
  yieldsK->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_KaonYields_log.png");

  c1->Clear();
  c1->SetLogy();
  piK_ratio->Draw("AP");
  //  texts.at(0)->Draw("same");
  for (int l=0; l<newDay.size(); l++){
    linesSingle.at(l)->Draw("same");
    textsSingle.at(l)->Draw("same");
  }
  c1->SaveAs( base + "_pion_Kaon_ratio_log.png");




  //runnumber
  c1->SetLogy(0);
  c1->cd();
  yieldsPi_run->Draw("AP");
  c1->SaveAs( base + "_PionYields_run.png");

  c1->Clear();
  yieldsPro_run->Draw("AP");
  c1->SaveAs( base + "_ProtonYields_run.png");
  
  c1->Clear();
  yieldsK_run->Draw("AP");
  c1->SaveAs( base + "_KaonYields_run.png");

  c1->Clear();
  piK_ratio_run->Draw("AP");
  c1->SaveAs( base + "_pion_Kaon_ratio_run.png");

  
  TF1 *yieldGausPi = new TF1("yieldGausPi","gaus(0)+gaus(3)");
  yieldGausPi->SetLineColor(kRed);

  TF1 *yieldGausK = new TF1("yieldGausK","gaus(0)");
  yieldGausK->SetLineColor(kRed);

  TF1 *yieldGausPro = new TF1("yieldGausPro","gaus(0)");
  yieldGausPro->SetLineColor(kRed);

  c1->SetLogy(0);
  //  gStyle->SetOptStat(1);
  c1->Clear();
  piYield_hist->Draw();
  gPad->Update();
  TPaveStats *stPi = (TPaveStats*)piYield_hist->FindObject("stats");
  stPi->SetX1NDC(0.2);
  stPi->SetX2NDC(0.5);  
  stPi->SetY1NDC(0.5);
  stPi->SetY2NDC(0.8);  
  yieldGausPi->SetParameters(15,2.1,0.1,20,2.8,0.1);
  piYield_hist->Fit("yieldGausPi");
  yieldGausPi->Draw("same");
  c1->SaveAs( base + "_pionGaus.png");

  piYield_hist->Write();
  yieldGausPi->Write("PionYieldFit");

  c1->Clear();
  kYield_hist->Draw();
  yieldGausK->SetParameters(11,1,0.5);
  kYield_hist->Fit("yieldGausK");
  yieldGausK->Draw("same");
  c1->SaveAs( base + "_kaonGaus.png");

  kYield_hist->Write();
  yieldGausK->Write("KaonYieldFit");

  
  c1->Clear();
  proYield_hist->Draw();
  yieldGausPro->SetParameters(11,0.5,0.5);
  proYield_hist->Fit("yieldGausPro");
  yieldGausPro->Draw("same");
  c1->SaveAs( base + "_protonGaus.png");

  proYield_hist->Write();
  yieldGausPro->Write("ProtonYieldFit");


  primaryCanvasYields->cd();
  gPad->SetTopMargin(0.35);
  gPad->SetBottomMargin(0.15);
  primaryCanvasYields->Divide(1,4,0.01,0);

  primaryCanvasYields->cd(1);
  yieldsPiY->GetYaxis()->SetTitle("Raw #pi Yield y_{CM}");
  yieldsPiY->GetYaxis()->SetTitleOffset(0.45);
  yieldsPiY->GetYaxis()->CenterTitle();
  yieldsPiY->GetYaxis()->SetTitleSize(0.1);
  yieldsPiY->GetYaxis()->SetLabelSize(0.08);
  yieldsPiY->Draw("AP");


  primaryCanvasYields->cd(2);
  yieldsPi->GetYaxis()->SetTitle("Raw #pi Yield");
  yieldsPi->GetYaxis()->SetTitleOffset(0.45);
  yieldsPi->GetYaxis()->CenterTitle();
  yieldsPi->GetYaxis()->SetTitleSize(0.1);
  yieldsPi->GetYaxis()->SetLabelSize(0.08);
  yieldsPi->Draw("AP");

  primaryCanvasYields->cd(3);
  yieldsPiEToF->Draw("AP");

  primaryCanvasYields->cd(4);
  yieldsPiToF->Draw("AP");


  primaryCanvasYields->cd(0);
  for (int l=0; l<newDay.size(); l++){
    linesGlobalPID.at(l)->Draw("same");
    textsGlobalPID.at(l)->Draw("same"); //l+1?
  }
  primaryTitle->Draw("same");
  dayText->Draw("same");
  primaryCanvasYields->SaveAs( base + "_PrimaryTrackYields.png");


  outFile->Close();

}

