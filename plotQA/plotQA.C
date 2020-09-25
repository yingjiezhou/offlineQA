#include "syang/headers.h"
#include "syang/function.C"

//QA event level cuts
const Double_t mVrCut = 2.;
const Double_t mVzCut = 200.;
const Double_t mVzDiffCut= 3.;

const Double_t epsilon = 1.e-10;

Bool_t dumpEmptyHT2Runs = kFALSE;
Bool_t setProfileRange  = kTRUE;
Bool_t detailBadRunInfo = kTRUE;

typedef map<Int_t, Int_t> IntMap;

IntMap mRunId2RunIdx;
IntMap mRunIdx2RunId;
IntMap currentBadRuns;   // none-zero HT2 triggered events, and none-zero good variable events
IntMap currentEmptyRuns; // none-zero HT2 triggered events, but zero good variable event 
IntMap newBadRuns;
IntMap newEmptyRuns;
IntMap totalBadRuns;
IntMap totalEmptyRuns;
IntMap totalRejectRuns; // totalBadRuns + totalEmptyRuns (won't count overlap runs between these two categories)
vector<TLine*> linesGlobal; // draw day info.
vector<TText*> textsGlobal;

const Int_t rebRefMult = 10;

const Double_t nSigma = 3.;
const Int_t nTrgIds = 1;
Int_t runStart[nTrgIds] ={21226023}; // run20
Int_t runStop[nTrgIds] = {21232025};


const Int_t    MStyle = 20;
const Double_t MSize = 0.6;
const Int_t    MColor = 1;
const Int_t    badMStyle = MStyle;
const Double_t badMSize = MSize;
const Int_t    badMColor = 2;
const Int_t    emptyMStyle = MStyle;
const Double_t emptyMSize = MSize;
const Int_t    emptyMColor = 6;

Bool_t Init(TString runlist="0911_production_7p7GeV_2020_runnumber_DD.dat");
Int_t  grabRunIdx(Int_t runId);
Int_t  grabRunId(Int_t runIdx);
Bool_t produceShade(TProfile *h, Int_t binLow, Int_t binHi, Double_t &mean, Double_t &rms, TGraphErrors *shade);
void   addShade(TProfile *h, Int_t binLow, Int_t binHi, Int_t trgIdx);
Bool_t dumpBadRuns(TString dirName, TString varName, Int_t &varIdx);
Int_t drawDay(Int_t runStartT, Int_t runStopT, TString runlistT,Int_t binLow, Int_t binHi);
void pdfAction(TCanvas *c, TPDF *ps);

TH1D *hnEvtsvsRun;
void plotQA(Int_t runStartT=21226023, Int_t runStopT=21232025, TString runlist="0827_production_26p5_2020_runnumber_DD.dat"){
  
  runStart[0] = runStartT;
  runStop[0] = runStopT;
  TString outfileName;
  outfileName = runlist;
  outfileName.ReplaceAll("_runnumber_DD.dat", "");
  //==================================================================
  
  cout<<"Start to draw the QA plots ..."<<endl;
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(1);
  
  TFile *f = new TFile(Form("rootfile/%s.root", outfileName.Data()));
  if(!f->IsOpen()){
    cout<<"Fail to load the QA root file!"<<endl;
    return;
  }
  
  if(!Init(runlist)) cout<<"Failed to initialize code !";
  
  TString dir = outfileName;
  system(Form("mkdir -p %s", dir.Data()));
  system(Form("rm -rf %s/*", dir.Data()));
  
  TString badRunDir = outfileName;
  badRunDir.Append("/badruns");
  system(Form("mkdir -p %s", badRunDir.Data()));
  system(Form("rm -rf %s/*", badRunDir.Data()));
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  Int_t nRaws = 2;
  Int_t nColumns = 2;
  Int_t nPads = nRaws*nColumns;
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->Divide(nColumns, nRaws);
  
  TPDF *ps = new TPDF(Form("%s/%s.pdf", dir.Data(),dir.Data()),111);
  ps->Off();
  //*********   import histograms of inclusive QA - event level  *********//
  // default vertex
  TH1D *hEvent_DefVtx           = (TH1D *)f->Get("hEvent_DefVtx");
  TH2D *hVyvsVx_DefVtx          = (TH2D *)f->Get("hVyvsVx_DefVtx");
  TH2D *hVpdVzvsTpcVz_DefVtx    = (TH2D *)f->Get("hVpdVzvsTpcVz_DefVtx");
  TH2D *hVzDiffvsTpcVz_DefVtx   = (TH2D *)f->Get("hVzDiffvsTpcVz_DefVtx");
  TH2D *hVzDiffvsRefMult_DefVtx = (TH2D *)f->Get("hVzDiffvsRefMult_DefVtx");
  TH1D *hRefMult_DefVtx         = (TH1D *)hVzDiffvsRefMult_DefVtx->ProjectionX("hRefMult_DefVtx",-1,-1);
  TH1D *hRefMult_VzVrCut_DefVtx = (TH1D *)f->Get("hRefMult_VzVrCut_DefVtx");
  TH1D *hRefMult_EvtCut_DefVtx  = (TH1D *)f->Get("hRefMult_EvtCut_DefVtx");
  
  TH1D *hTpcVz_DefVtx  = (TH1D *)hVpdVzvsTpcVz_DefVtx->ProjectionX("hTpcVz_DefVtx",-1,-1);
  TH1D *hVzDiff_DefVtx = (TH1D *)hVzDiffvsTpcVz_DefVtx->ProjectionY("hVzDiff_DefVtx",-1,-1);
  
  // selected vertex
  TH1D *hEvent                     = (TH1D *)f->Get("hEvent");
  TH2D *hVyvsVx                    = (TH2D *)f->Get("hVyvsVx");
  TH2D *hVpdVzvsTpcVz              = (TH2D *)f->Get("hVpdVzvsTpcVz");
  TH2D *hTpcVzvsRefMult            = (TH2D *)f->Get("hTpcVzvsRefMult");
  TH2D *hRawVpdVzvsRefMult         = (TH2D *)f->Get("hRawVpdVzvsRefMult");
  TH2D *hVpdVzvsRefMult            = (TH2D *)f->Get("hVpdVzvsRefMult");
  TH2D *hVzDiffvsTpcVz             = (TH2D *)f->Get("hVzDiffvsTpcVz");
  TH2D *hVzDiffvsRefMult           = (TH2D *)f->Get("hVzDiffvsRefMult");
  TH2D *hRawVzDiffvsRefMult        = (TH2D *)f->Get("hRawVzDiffvsRefMult");
  TH1D *hRefMult                   = (TH1D *)hVzDiffvsRefMult->ProjectionX("hRefMult",-1,-1);
  TH1D *hRefMult_VzVrCut           = (TH1D *)f->Get("hRefMult_VzVrCut");
  TH1D *hRefMult_EvtCut            = (TH1D *)f->Get("hRefMult_EvtCut");
  TH2D *hVtxIdxvsRefMult           = (TH2D *)f->Get("hVtxIdxvsRefMult");
  TH2D *hVtxIdxvsRefMult_VzDiffCut = (TH2D *)f->Get("hVtxIdxvsRefMult_VzDiffCut");
  TH2D *hRefMultvsZdcX             = (TH2D *)f->Get("hRefMultvsZdcX");
  TH2D *hRefMultvsBbcX             = (TH2D *)f->Get("hRefMultvsBbcX");
  
  TH1D *hTpcVz  = (TH1D *)hVpdVzvsTpcVz->ProjectionX("hTpcVz",-1,-1);
  TH1D *hVzDiff = (TH1D *)hVzDiffvsTpcVz->ProjectionY("hVzDiff",-1,-1);
  
  TLegend *leg = new TLegend(0.14, 0.7, 0.55, 0.87);
  leg->SetTextSize(0.02);
  leg->SetFillStyle (0);
  leg->SetFillColor (0);
  leg->SetBorderSize(0);
  
  TF1 *g = new TF1("g", "gaus", -1000, 1000);
  setFun(g, 2, 2, 1);
  
  c1->cd();
  hEvent_DefVtx->GetXaxis()->SetBinLabel(1, "All events");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(3, "NPE_18");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(4, "NPE_25");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(5, "NPE_25_nozdc");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
  hEvent_DefVtx->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mVrCut));
  hEvent_DefVtx->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mVzCut));
  hEvent_DefVtx->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mVzDiffCut));
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(13, "NPE_18");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(14, "NPE_25");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(15, "NPE_25_nozdc");
  setHisto(hEvent_DefVtx, 20, 1., 1, 1, 2);
  setHisto(hEvent, 24, 1., 2, 2, 2);
  hEvent_DefVtx->Draw("hist");
  hEvent->Draw("histtextsame");
  //  leg->AddEntry(hEvent_DefVtx, "Default vertex", "l");
  //  leg->AddEntry(hEvent, "Selected vertex", "l"); // selected?
  //  leg->Draw("same");
  //  drawLatex(0.165, 0.45, "For selected vertex scenario:", 22, 0.05, 1);
  //  drawLatex(0.165, 0.4, Form("HT2 statistics: %1.1f M", hEvent->GetBinContent(3)/1.e6), 22, 0.045, 4);
  //  drawLatex(0.165, 0.35, Form("+ None-zero Vertex: %1.2f%%", hEvent->GetBinContent(8)*1./hEvent->GetBinContent(3)*100), 22, 0.045, 4);
  drawLatex(0.165, 0.30, Form("+ V_{r}<%1.1f cm: %1.2f%%", mVrCut, hEvent->GetBinContent(9)*1./hEvent->GetBinContent(1)*100), 22, 0.045, 4);
  drawLatex(0.165, 0.25, Form("+ |TPC V_{z}|<%1.1f cm: %1.2f%%", mVzCut, hEvent->GetBinContent(10)*1./hEvent->GetBinContent(1)*100), 22, 0.045, 4);
  drawLatex(0.165, 0.20, Form("+ |TPC V_{z} - VPD V_{z}|<%1.1f cm: %1.2f%%", mVzDiffCut, hEvent->GetBinContent(11)*1./hEvent->GetBinContent(1)*100), 22, 0.045, 4);
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/evtStats.png", dir.Data()));
  
  Int_t vzBinLow  = hTpcVz_DefVtx->GetXaxis()->FindBin(-mVzCut+1.e-6);
  Int_t vzBinHi   = hTpcVz_DefVtx->GetXaxis()->FindBin(mVzCut-1.e-6);
  Double_t vzFrac = hTpcVz_DefVtx->Integral(vzBinLow, vzBinHi)*1./hTpcVz_DefVtx->Integral(-1,-1);
  Int_t vzDiffBinLow  = hVzDiff_DefVtx->GetXaxis()->FindBin(-mVzDiffCut+1.e-6);
  Int_t vzDiffBinHi   = hVzDiff_DefVtx->GetXaxis()->FindBin(mVzDiffCut-1.e-6);
  Double_t vzDiffFrac = hVzDiff_DefVtx->Integral(vzDiffBinLow, vzDiffBinHi)*1./hVzDiff_DefVtx->Integral(-1,-1);
  
  c2->cd(1);
  gPad->SetLogz(1);
  hVyvsVx_DefVtx->Rebin2D(10, 10);
  hVyvsVx_DefVtx->Draw("colz");
  //  drawLatex(0.16, 0.8, "Default Vertex", 22, 0.12, 1);
  c2->cd(2);
  gPad->SetLogz(1);
  hVpdVzvsTpcVz_DefVtx->Rebin2D(5, 5);
  hVpdVzvsTpcVz_DefVtx->Draw("colz");
  c2->cd(3);
  gPad->SetLogy(1);
  hVzDiff_DefVtx->SetLineColor(1);
  hVzDiff_DefVtx->Fit(g, "", "", -2, 2);
  hVzDiff_DefVtx->Draw("same");
  drawLine(-mVzDiffCut, 0, -mVzDiffCut, hVzDiff_DefVtx->GetMaximum()*0.8, 2, 2, 2);
  drawLine( mVzDiffCut, 0,  mVzDiffCut, hVzDiff_DefVtx->GetMaximum()*0.8, 2, 2, 2);
  drawLatex(0.16, 0.82, Form("Vz Cut Eff.: %1.2f%%", vzFrac*100), 22, 0.05, 1);
  drawLatex(0.16, 0.76, Form("VzDiff Cut Eff.: %1.2f%%", vzDiffFrac*100), 22, 0.05, 1);
  pdfAction(c2, ps);
  c2->SaveAs(Form("%s/defVtx.png", dir.Data()));
  
  hVzDiffvsTpcVz->RebinX(10);
  hVzDiffvsTpcVz_DefVtx->RebinX(10);
  Int_t nVzPoints  = hVzDiffvsTpcVz->GetNbinsX();
  TGraphErrors *grVpdVzResvsTpcVz        = new TGraphErrors(nVzPoints);
  TGraphErrors *grVpdVzResvsTpcVz_DefVtx = new TGraphErrors(nVzPoints);
  //  for(Int_t i=0; i<nVzPoints; i++){
  //    Double_t nNumErr, nDenErr;
  //
  //    TH1D *hTempVzDiff = (TH1D *)hVzDiffvsTpcVz->ProjectionY("hTempVzDiff", i+1, i+1);
  //    hTempVzDiff->Fit(g, "Q", "", -2.5, 2.5);
  //    Double_t mean     = g->GetParameter(1);
  //    Double_t sigma    = g->GetParameter(2);
  //    Double_t sigmaErr = g->GetParError(2);
  //    if(mean-3*sigma>-2.5 && mean+3*sigma<2.5){
  //      hTempVzDiff->Fit(g, "Q", "", mean-3*sigma, mean+3*sigma);
  //    }
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    hTempVzDiff->Delete();
  //    grVpdVzResvsTpcVz->SetPoint(i, hVzDiffvsTpcVz->GetXaxis()->GetBinCenter(i+1), sigma);
  //    grVpdVzResvsTpcVz->SetPointError(i, hVzDiffvsTpcVz->GetXaxis()->GetBinWidth(i+1)/2., sigmaErr);
  //
  //    hTempVzDiff = (TH1D *)hVzDiffvsTpcVz_DefVtx->ProjectionY("hTempVzDiff", i+1, i+1);
  //    hTempVzDiff->Fit(g, "Q", "", -2.5, 2.5);
  //    mean     = g->GetParameter(1);
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    if(mean-3*sigma>-2.5 && mean+3*sigma<2.5){
  //      hTempVzDiff->Fit(g, "Q", "", mean-3*sigma, mean+3*sigma);
  //    }
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    hTempVzDiff->Delete();
  //    grVpdVzResvsTpcVz_DefVtx->SetPoint(i, hVzDiffvsTpcVz_DefVtx->GetXaxis()->GetBinCenter(i+1), sigma);
  //    grVpdVzResvsTpcVz_DefVtx->SetPointError(i, hVzDiffvsTpcVz_DefVtx->GetXaxis()->GetBinWidth(i+1)/2., sigmaErr);
  //  }
  
  hVzDiffvsRefMult->RebinX(rebRefMult);
  hVzDiffvsRefMult_DefVtx->RebinX(rebRefMult);
  Int_t nRefMultPoints = hVzDiffvsRefMult->GetXaxis()->FindBin(650-1.e-3) - hVzDiffvsRefMult->GetXaxis()->FindBin(0+1.e-3) + 1;
  TGraphErrors *grVpdVzResvsRefMult        = new TGraphErrors(nRefMultPoints);
  TGraphErrors *grVpdVzResvsRefMult_DefVtx = new TGraphErrors(nRefMultPoints);
  //  for(Int_t i=0; i<nRefMultPoints; i++){
  //    Double_t nNumErr, nDenErr;
  //
  //    TH1D *hTempVzDiff = (TH1D *)hVzDiffvsRefMult->ProjectionY("hTempVzDiff", i+1, i+1);
  //    hTempVzDiff->Fit(g, "Q", "", -2.5, 2.5);
  //    Double_t mean     = g->GetParameter(1);
  //    Double_t sigma    = g->GetParameter(2);
  //    Double_t sigmaErr = g->GetParError(2);
  //    if(mean-3*sigma>-2.5 && mean+3*sigma<2.5){
  //      hTempVzDiff->Fit(g, "Q", "", mean-3*sigma, mean+3*sigma);
  //    }
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    grVpdVzResvsRefMult->SetPoint(i, hVzDiffvsRefMult->GetXaxis()->GetBinCenter(i+1), sigma);
  //    grVpdVzResvsRefMult->SetPointError(i, hVzDiffvsRefMult->GetXaxis()->GetBinWidth(i+1)/2., sigmaErr);
  //
  //    hTempVzDiff = (TH1D *)hVzDiffvsRefMult_DefVtx->ProjectionY("hTempVzDiff", i+1, i+1);
  //    hTempVzDiff->Fit(g, "Q", "", -2.5, 2.5);
  //    mean     = g->GetParameter(1);
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    if(mean-3*sigma>-2.5 && mean+3*sigma<2.5){
  //      hTempVzDiff->Fit(g, "Q", "", mean-3*sigma, mean+3*sigma);
  //    }
  //    sigma    = g->GetParameter(2);
  //    sigmaErr = g->GetParError(2);
  //    hTempVzDiff->Delete();
  //    grVpdVzResvsRefMult_DefVtx->SetPoint(i, hVzDiffvsRefMult_DefVtx->GetXaxis()->GetBinCenter(i+1), sigma);
  //    grVpdVzResvsRefMult_DefVtx->SetPointError(i, hVzDiffvsRefMult_DefVtx->GetXaxis()->GetBinWidth(i+1)/2., sigmaErr);
  //  }
  
  TH2D *ddVpdVzResvsTpcVz = (TH2D *)histo("ddVpdVzResvsTpcVz", -200, 200, 0.3, 0.7,"TPC Vz (cm)","VpdVz Resolution (cm)");
  ddVpdVzResvsTpcVz->GetXaxis()->SetTitleSize(0.06);
  ddVpdVzResvsTpcVz->GetYaxis()->SetNdivisions(510);
  
  TH2D *ddVpdVzResvsRefMult = (TH2D *)histo("ddVpdVzResvsRefMult", 0, 800, 0, 1.6,"refMult","VpdVz Resolution (cm)");
  ddVpdVzResvsRefMult->GetXaxis()->SetTitleSize(0.06);
  
  c2->Clear();
  c2->Divide(nColumns, nRaws);
  
  c2->cd(1);
  gPad->SetLogy(0);
  gPad->SetLogz(1);
  hVzDiffvsTpcVz_DefVtx->Draw("colz");
  c2->cd(2);
  gPad->SetLogy(0);
  ddVpdVzResvsTpcVz->Draw("c");
  setGraph(grVpdVzResvsTpcVz_DefVtx, 20, 0.6, 1, 1, 2);
  grVpdVzResvsTpcVz_DefVtx->Draw("pzsame");
  c2->cd(3);
  gPad->SetLogy(0);
  gPad->SetLogz(1);
  hVzDiffvsRefMult_DefVtx->Draw("colz");
  c2->cd(4);
  gPad->SetLogy(0);
  ddVpdVzResvsRefMult->Draw("c");
  setGraph(grVpdVzResvsRefMult_DefVtx, 20, 0.6, 1, 1, 2);
  grVpdVzResvsRefMult_DefVtx->Draw("pzsame");
  pdfAction(c2, ps);
  c2->SaveAs(Form("%s/vpdVzRes_defVtx.png", dir.Data()));
  
  c1->cd();
  gPad->SetLogy(1);
  hRefMult_DefVtx->GetXaxis()->SetRangeUser(0, 800);
  hRefMult_DefVtx->SetLineColor(1);
  hRefMult_VzVrCut_DefVtx->SetLineColor(2);
  hRefMult_EvtCut_DefVtx->SetLineColor(4);
  hRefMult_DefVtx->Draw();
  hRefMult_VzVrCut_DefVtx->Draw("same");
  hRefMult_EvtCut_DefVtx->Draw("same");
  setLegend(leg, 0.36, 0.7, 0.55, 0.87, 0.05);
  leg->AddEntry(hRefMult_DefVtx, "w/o event cuts", "l");
  leg->AddEntry(hRefMult_VzVrCut_DefVtx,Form("|V_{z}|#leq%1.0f cm & |V_{r}|#leq%1.0f cm", mVzCut, mVrCut), "l");
  leg->AddEntry(hRefMult_EvtCut_DefVtx, Form("|V_{z}|#leq%1.0f cm & |V_{r}|#leq%1.0f cm & |V_{z}Diff|#leq%1.0f cm", mVzCut, mVrCut, mVzDiffCut), "l");
  leg->Draw("same");
  //  drawLatex(0.20, 0.48, "Default Vertex", 22, 0.06, 1);
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/refMult_defVtx.png", dir.Data()));
  //
  //  vzBinLow  = hTpcVz->GetXaxis()->FindBin(-mVzCut+1.e-6);
  //  vzBinHi   = hTpcVz->GetXaxis()->FindBin(mVzCut-1.e-6);
  //  vzFrac    = hTpcVz->Integral(vzBinLow, vzBinHi)*1./hTpcVz->Integral(-1,-1);
  //  vzDiffBinLow  = hVzDiff->GetXaxis()->FindBin(-mVzDiffCut+1.e-6);
  //  vzDiffBinHi   = hVzDiff->GetXaxis()->FindBin(mVzDiffCut-1.e-6);
  //  vzDiffFrac    = hVzDiff->Integral(vzDiffBinLow, vzDiffBinHi)*1./hVzDiff->Integral(-1,-1);
  //
  //  clearPad(c2, nPads);
  //  c2->cd(1);
  //  gPad->SetLogz(1);
  //  hVyvsVx->Rebin2D(10, 10);
  //  hVyvsVx->Draw("colz");
  //  drawLatex(0.16, 0.8, "Selected Vertex", 22, 0.12, 1);
  //  c2->cd(2);
  //  gPad->SetLogz(1);
  //  hVpdVzvsTpcVz->Rebin2D(5, 5);
  //  hVpdVzvsTpcVz->Draw("colz");
  //  c2->cd(3);
  //  gPad->SetLogy(1);
  //  hVzDiff->SetLineColor(1);
  //  hVzDiff->Fit(g, "", "", -2, 2);
  //  hVzDiff->Draw("same");
  //  drawLine(-mVzDiffCut, 0, -mVzDiffCut, hVzDiff->GetMaximum()*0.05, 2, 2, 2);
  //  drawLine( mVzDiffCut, 0,  mVzDiffCut, hVzDiff->GetMaximum()*0.05, 2, 2, 2);
  //  drawLatex(0.16, 0.82, Form("Vz Cut Eff.: %1.2f%%", vzFrac*100), 22, 0.05, 1);
  //  drawLatex(0.16, 0.76, Form("VzDiff Cut Eff.: %1.2f%%", vzDiffFrac*100), 22, 0.05, 1);
  //  pdfAction(c2, ps);
  //  c2->SaveAs(Form("%s/selVertex.png", dir.Data()));
  
  //  c1->cd();
  //  gPad->SetLogy(0);
  //  gPad->SetLogz(1);
  //  hTpcVzvsRefMult->Draw("colz");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/tpcVzvsRefMult_selVtx.png", dir.Data()));
  //
  c1->cd();
  gPad->SetLogy(0);
  gPad->SetLogz(1);
  hRawVpdVzvsRefMult->Draw("colz");
  drawLine(0, -500, 700, -500, 2, 2, 2);
  drawLine(0,  500, 700,  500, 2, 2, 2);
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/vpdVzvsRefMult.png", dir.Data()));
  
  //  c2->cd(1);
  //  gPad->SetLogy(0);
  //  gPad->SetLogz(1);
  //  hVzDiffvsTpcVz->Draw("colz");
  //  c2->cd(2);
  //  gPad->SetLogy(0);
  //  ddVpdVzResvsTpcVz->Draw("c");
  //  setGraph(grVpdVzResvsTpcVz, 20, 0.6, 1, 1, 2);
  //  grVpdVzResvsTpcVz->Draw("pzsame");
  //  c2->cd(3);
  //  gPad->SetLogy(0);
  //  gPad->SetLogz(1);
  //  hVzDiffvsRefMult->Draw("colz");
  //  c2->cd(4);
  //  gPad->SetLogy(0);
  //  ddVpdVzResvsRefMult->Draw("c");
  //  setGraph(grVpdVzResvsRefMult, 20, 0.6, 1, 1, 2);
  //  grVpdVzResvsRefMult->Draw("pzsame");
  //  pdfAction(c2, ps);
  //  c2->SaveAs(Form("%s/vpdVzRes_selVtx.png", dir.Data()));
  
  //  c1->cd();
  //  gPad->SetLogy(1);
  //  hRefMult->GetXaxis()->SetRangeUser(0, 800);
  //  hRefMult->SetLineColor(1);
  //  hRefMult_VzVrCut->SetLineColor(2);
  //  hRefMult_EvtCut->SetLineColor(4);
  //  hRefMult->Draw();
  //  hRefMult_VzVrCut->Draw("same");
  //  hRefMult_EvtCut->Draw("same");
  //  leg->Draw("same");
  //  drawLatex(0.20, 0.48, "Selected Vertex", 22, 0.06, 1);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/refMult_selVtx.png", dir.Data()));
  
  //  hRawVpdVzvsRefMult->RebinX(rebRefMult);
  //  hVtxIdxvsRefMult->RebinX(rebRefMult);
  //  hVtxIdxvsRefMult_VzDiffCut->RebinX(rebRefMult);
  //  Int_t nPoints = hVtxIdxvsRefMult->GetXaxis()->FindBin(650-1.e-3) - hVtxIdxvsRefMult->GetXaxis()->FindBin(0+1.e-3) + 1;
  //  TGraphErrors *grVpdEffvsRefMult     = new TGraphErrors(nPoints);
  //  TGraphErrors *grDefVtxFracvsRefMult = new TGraphErrors(nPoints);
  //  TGraphErrors *grDefVtxFracvsRefMult_VzDiffCut = new TGraphErrors(nPoints);
  //  for(Int_t i=0; i<nPoints; i++){
  //    Double_t nNumErr, nDenErr;
  //
  //    Int_t vpdVzBinLow = hRawVpdVzvsRefMult->GetYaxis()->FindBin(-500+1.e-5);
  //    Int_t vpdVzBinHi  = hRawVpdVzvsRefMult->GetYaxis()->FindBin(500-1.e-5);
  //    Int_t nYBins      = hRawVpdVzvsRefMult->GetNbinsY();
  //    Double_t nNum = hRawVpdVzvsRefMult->IntegralAndError(i+1,i+1,vpdVzBinLow,vpdVzBinHi,nNumErr);
  //    Double_t nDen = hRawVpdVzvsRefMult->IntegralAndError(i+1,i+1,0,nYBins+1,nDenErr);
  //    Double_t ratio    = nNum/nDen;
  //    Double_t ratioErr = TMath::Abs( ((1.-2.*ratio)*pow(nNumErr,2)+pow(ratio*nDenErr,2))/pow(nDen,2) );
  //    ratioErr = sqrt(ratioErr);
  //    grVpdEffvsRefMult->SetPoint(i, hRawVpdVzvsRefMult->GetXaxis()->GetBinCenter(i+1), ratio);
  //    grVpdEffvsRefMult->SetPointError(i, hRawVpdVzvsRefMult->GetXaxis()->GetBinWidth(i+1)/2., ratioErr);
  //
  //    nYBins = hVtxIdxvsRefMult->GetNbinsY();
  //    nNum = hVtxIdxvsRefMult->IntegralAndError(i+1,i+1,1,1,nNumErr);
  //    nDen = hVtxIdxvsRefMult->IntegralAndError(i+1,i+1,1,nYBins,nDenErr);
  //    ratio    = nNum/nDen;
  //    ratioErr = TMath::Abs( ((1.-2.*ratio)*pow(nNumErr,2)+pow(ratio*nDenErr,2))/pow(nDen,2) );
  //    ratioErr = sqrt(ratioErr);
  //    grDefVtxFracvsRefMult->SetPoint(i, hVtxIdxvsRefMult->GetXaxis()->GetBinCenter(i+1), ratio);
  //    grDefVtxFracvsRefMult->SetPointError(i, hVtxIdxvsRefMult->GetXaxis()->GetBinWidth(i+1)/2., ratioErr);
  //
  //    nYBins = hVtxIdxvsRefMult_VzDiffCut->GetNbinsY();
  //    nNum = hVtxIdxvsRefMult_VzDiffCut->IntegralAndError(i+1,i+1,1,1,nNumErr);
  //    nDen = hVtxIdxvsRefMult_VzDiffCut->IntegralAndError(i+1,i+1,1,nYBins,nDenErr);
  //    ratio    = nNum/nDen;
  //    ratioErr = TMath::Abs( ((1.-2.*ratio)*pow(nNumErr,2)+pow(ratio*nDenErr,2))/pow(nDen,2) );
  //    ratioErr = sqrt(ratioErr);
  //    grDefVtxFracvsRefMult_VzDiffCut->SetPoint(i, hVtxIdxvsRefMult_VzDiffCut->GetXaxis()->GetBinCenter(i+1), ratio);
  //    grDefVtxFracvsRefMult_VzDiffCut->SetPointError(i, hVtxIdxvsRefMult_VzDiffCut->GetXaxis()->GetBinWidth(i+1)/2., ratioErr);
  //  }
  //
  //  TH2D *ddFrac = (TH2D *)histo("ddFrac", 0, 650, 0.2, 1.05,"refMult","Default Vertex Fraction");
  //  ddFrac->GetXaxis()->SetTitleSize(0.06);
  //  ddFrac->GetYaxis()->SetTitleSize(0.055);
  //  c1->cd();
  //  gPad->SetLogy(0);
  //  ddFrac->Draw("c");
  //  //setGraph(grDefVtxFracvsRefMult, 20, 1., 1, 1, 2);
  //  setGraph(grDefVtxFracvsRefMult_VzDiffCut, 20, 1., 1, 1, 2);
  //  //grDefVtxFracvsRefMult->Draw("pzsame");
  //  grDefVtxFracvsRefMult_VzDiffCut->Draw("pzsame");
  //  drawLatex(0.24, 0.32, Form("|TPC V_{z} - VPD V_{z}|<%1.1f cm", mVzDiffCut), 22, 0.05, 1);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/defVtxFraction_selVtx.png", dir.Data()));
  
  //  c1->cd();
  //  gPad->SetLogy(0);
  //  TH1D *hRefMultRatio = (TH1D *)hRefMult_EvtCut->Clone("hRefMultRatio");
  //  hRefMultRatio->Reset();
  //  hRefMultRatio->Divide(hRefMult_EvtCut, hRefMult_VzVrCut, 1, 1, "B");
  //  setHisto(hRefMultRatio, 20, 1., 1, 1);
  //  setGraph(grVpdEffvsRefMult, 29, 1.5, 4, 4, 2);
  //  hRefMultRatio->SetAxisRange(0, 200, "X");
  //  hRefMultRatio->SetAxisRange(0, 1.05, "Y");
  //  hRefMultRatio->GetYaxis()->SetTitle("Ratio (Efficiency)");
  //  hRefMultRatio->Draw("p");
  //  //grVpdEffvsRefMult->Draw("pzsame");
  //  drawLine(10, 0, 10, 1, 2, 2, 2);
  //  drawLine(22, 0, 22, 1, 2, 2, 2);
  //  drawLine(43, 0, 43, 1, 2, 2, 2);
  //  drawLatex(0.42, 0.64, "Selected Vertex", 22, 0.06, 1);
  //  drawLatex(0.195, 0.16, "70-80%", 22, 0.05, 2, 90);
  //  drawLatex(0.260, 0.16, "60-70%", 22, 0.05, 2, 90);
  //  setLegend(leg, 0.32, 0.2, 0.8, 0.4, 0.05);
  //  leg->AddEntry(hRefMultRatio, "#frac{Vz+Vr+VzDiff Cuts}{Vz+Vr Cuts}", "pl" );
  //  leg->AddEntry(grVpdEffvsRefMult, "VPD Efficiency", "pl");
  //  //leg->Draw("same");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/refMultRatio_selVtx.png", dir.Data()));
  //
  //  TH1D *hSelOverDef_RefMultRatio = (TH1D *)hRefMult_EvtCut->Clone("hSelOverDef_RefMultRatio");
  //  hSelOverDef_RefMultRatio->Reset();
  //  hSelOverDef_RefMultRatio->Divide(hRefMult_EvtCut, hRefMult_EvtCut_DefVtx, 1, 1, "B");
  //  setHisto(hSelOverDef_RefMultRatio, 20, 0.6, 1, 1);
  //  hSelOverDef_RefMultRatio->SetAxisRange(0, 600, "X");
  //  hSelOverDef_RefMultRatio->SetAxisRange(0.8, 1.2, "Y");
  //  hSelOverDef_RefMultRatio->GetYaxis()->SetTitle("SelectedVtx / DefaultVtx");
  //  hSelOverDef_RefMultRatio->Draw("p");
  //  drawLatex(0.2, 0.32, "w/ event level cuts", 22, 0.06, 1);
  //  //drawLine(10, 0, 10, 1, 2, 2, 2);
  //  //drawLine(22, 0, 22, 1, 2, 2, 2);
  //  //drawLine(43, 0, 43, 1, 2, 2, 2);
  //  //drawLatex(0.42, 0.48, "Selected Vertex", 22, 0.06, 1);
  //  //drawLatex(0.195, 0.16, "70-80%", 22, 0.05, 2, 90);
  //  //drawLatex(0.260, 0.16, "60-70%", 22, 0.05, 2, 90);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/refMultRatio_selVtxOverdelVtx.png", dir.Data()));
  
  //*********   import histograms of inclusive QA - track level  *********//
  TH1D *hNHitsFit       = (TH1D *)f->Get("hNHitsFit");
  TH1D *hNHitsPoss      = (TH1D *)f->Get("hNHitsPoss");
  TH1D *hNHitsDedx      = (TH1D *)f->Get("hNHitsDedx");
  TH2D *hEtavsPt        = (TH2D *)f->Get("hEtavsPt");
  TH2D *hFFPhivsPt      = (TH2D *)f->Get("hFFPhivsPt");
  TH2D *hRFFPhivsPt     = (TH2D *)f->Get("hRFFPhivsPt");
  TH2D *hPosTrkEtavsPhi = (TH2D *)f->Get("hPosTrkEtavsPhi");
  TH2D *hNegTrkEtavsPhi = (TH2D *)f->Get("hNegTrkEtavsPhi");
  TH2D *hDcavsPt        = (TH2D *)f->Get("hDcavsPt");
  TH2D *hdEdxvsP        = (TH2D *)f->Get("hdEdxvsP");
  TH2D *hdEdxvsPhi      = (TH2D *)f->Get("hdEdxvsPhi");
  TH2D *hdEdxvsEta      = (TH2D *)f->Get("hdEdxvsEta");
  TH2D *hBetavsP        = (TH2D *)f->Get("hBetavsP");
  TH2D *hBEMCeEtavsPt   = (TH2D *)f->Get("hBEMCeEtavsPt");
  TH2D *hBEMCePhivsPt   = (TH2D *)f->Get("hBEMCePhivsPt");
  TH2D *hBEMCeEtavsPhi  = (TH2D *)f->Get("hBEMCeEtavsPhi");
  
  c1->cd();
  gPad->SetLogy(1);
  hNHitsFit->Draw();
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsFit.png", dir.Data()));
  
  c1->cd();
  hNHitsPoss->Draw();
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsPoss.png", dir.Data()));
  
  c1->cd();
  hNHitsDedx->SetAxisRange(-40, 40, "X");
  hNHitsDedx->Draw();
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsDedx.png", dir.Data()));
  
  c1->cd();
  gPad->SetLogy(0);
  gPad->SetLogz(1);
  hEtavsPt->RebinY(2);
  hEtavsPt->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/etavsPt.png", dir.Data()));
  
  //  c1->cd();
  //  hFFPhivsPt->Rebin2D();
  //  hFFPhivsPt->Draw("colz");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/phivsPt_FF.png", dir.Data()));
  
  c1->cd();
  hRFFPhivsPt->Rebin2D();
  hRFFPhivsPt->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/phivsPt_RFF.png", dir.Data()));
  
  c1->cd();
  hPosTrkEtavsPhi->Rebin2D(4, 2);
  hPosTrkEtavsPhi->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/posTrkEtavsPhi.png", dir.Data()));
  
  c1->cd();
  hNegTrkEtavsPhi->Rebin2D(4, 2);
  hNegTrkEtavsPhi->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/negTrkEtavsPhi.png", dir.Data()));
  
  c1->cd();
  hDcavsPt->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dcavsPt.png", dir.Data()));
  
  c1->cd();
  hdEdxvsP->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dEdxvsP.png", dir.Data()));
  
  c1->cd();
  hdEdxvsPhi->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dEdxvsPhi.png", dir.Data()));
  
  c1->cd();
  hdEdxvsEta->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dEdxvsEta.png", dir.Data()));
  
  c1->cd();
  hBetavsP->Draw("colz");
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/betavsP.png", dir.Data()));
  
  //  c1->cd();
  //  hBEMCeEtavsPt->Draw("colz");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigBEMCe_etavsPt.png", dir.Data()));
  //
  //  c1->cd();
  //  hBEMCePhivsPt->Draw("colz");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigBEMCe_phivsPt.png", dir.Data()));
  //
  //  c1->cd();
  //  hBEMCeEtavsPhi->Draw("colz");
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigBEMCe_etavsPhi.png", dir.Data()));
  
  //*********   import histograms of run by  run QA  *********//
  hnEvtsvsRun              = (TH1D *)f->Get("hnEvtsvsRun");
  TProfile *hBFieldvsRun         = (TProfile *)f->Get("hBFieldvsRun");
  TProfile *hZdcXvsRun           = (TProfile *)f->Get("hZdcXvsRun");
  TProfile *hBbcXvsRun           = (TProfile *)f->Get("hBbcXvsRun");
  TProfile *hZdcXoverBbcXvsRun   = (TProfile *)f->Get("hZdcXoverBbcXvsRun");
  TProfile *hTpcVxvsRun          = (TProfile *)f->Get("hTpcVxvsRun");
  TProfile *hTpcVyvsRun          = (TProfile *)f->Get("hTpcVyvsRun");
  TProfile *hTpcVzvsRun          = (TProfile *)f->Get("hTpcVzvsRun");
  TProfile *hVpdVzvsRun          = (TProfile *)f->Get("hVpdVzvsRun");
  TProfile *hVzDiffvsRun         = (TProfile *)f->Get("hVzDiffvsRun");
  TProfile *hRefMultvsRun        = (TProfile *)f->Get("hRefMultvsRun");
  TProfile *hNTrksvsRun          = (TProfile *)f->Get("hNTrksvsRun");
  TProfile *hPtvsRun             = (TProfile *)f->Get("hPtvsRun");
  TProfile *hEtavsRun            = (TProfile *)f->Get("hEtavsRun");
  TProfile *hPhivsRun            = (TProfile *)f->Get("hPhivsRun");
  TProfile *hDcavsRun            = (TProfile *)f->Get("hDcavsRun");
  TProfile *hNHitsFitvsRun       = (TProfile *)f->Get("hNHitsFitvsRun");
  TProfile *hNHitsPossvsRun      = (TProfile *)f->Get("hNHitsPossvsRun");
  TProfile *hNHitsDedxvsRun      = (TProfile *)f->Get("hNHitsDedxvsRun");
  TProfile *hDedxvsRun           = (TProfile *)f->Get("hDedxvsRun");
  TProfile *hNSigmaEvsRun        = (TProfile *)f->Get("hNSigmaEvsRun");
  TProfile *hNSigmaPivsRun       = (TProfile *)f->Get("hNSigmaPivsRun");
  TProfile *hNSigmaKvsRun        = (TProfile *)f->Get("hNSigmaKvsRun");
  TProfile *hNSigmaPvsRun        = (TProfile *)f->Get("hNSigmaPvsRun");
  TProfile *hBetavsRun           = (TProfile *)f->Get("hBetavsRun");
  TProfile *hNMthTrksvsRun       = (TProfile *)f->Get("hNMthTrksvsRun");
  TProfile *hMthTrkPtvsRun       = (TProfile *)f->Get("hMthTrkPtvsRun");
  TProfile *hMthTrkEtavsRun      = (TProfile *)f->Get("hMthTrkEtavsRun");
  TProfile *hMthTrkPhivsRun      = (TProfile *)f->Get("hMthTrkPhivsRun");
  TProfile *hMthTrkNSigmaEvsRun  = (TProfile *)f->Get("hMthTrkNSigmaEvsRun");
  TProfile *hMthTrkBetavsRun     = (TProfile *)f->Get("hMthTrkBetavsRun");
  TProfile *hMthTrkAdc0vsRun     = (TProfile *)f->Get("hMthTrkAdc0vsRun");
  TProfile *hMthTrkE0vsRun       = (TProfile *)f->Get("hMthTrkE0vsRun");
  TProfile *hMthTrkEvsRun        = (TProfile *)f->Get("hMthTrkEvsRun");
  TProfile *hMthTrkZDistvsRun    = (TProfile *)f->Get("hMthTrkZDistvsRun");
  TProfile *hMthTrkPhiDistvsRun  = (TProfile *)f->Get("hMthTrkPhiDistvsRun");
  TProfile *hMthTrkNEtavsRun     = (TProfile *)f->Get("hMthTrkNEtavsRun");
  TProfile *hMthTrkNPhivsRun     = (TProfile *)f->Get("hMthTrkNPhivsRun");
  TProfile *hNTrigTrksvsRun      = (TProfile *)f->Get("hNTrigTrksvsRun");
  TProfile *hTrigTrkPtvsRun      = (TProfile *)f->Get("hTrigTrkPtvsRun");
  TProfile *hTrigTrkEtavsRun     = (TProfile *)f->Get("hTrigTrkEtavsRun");
  TProfile *hTrigTrkPhivsRun     = (TProfile *)f->Get("hTrigTrkPhivsRun");
  TProfile *hTrigTrkNSigmaEvsRun = (TProfile *)f->Get("hTrigTrkNSigmaEvsRun");
  TProfile *hTrigTrkAdc0vsRun    = (TProfile *)f->Get("hTrigTrkAdc0vsRun");
  TProfile *hTrigTrkE0vsRun      = (TProfile *)f->Get("hTrigTrkE0vsRun");
  TProfile *hTrigTrkEvsRun       = (TProfile *)f->Get("hTrigTrkEvsRun");
  TProfile *hTrigTrkZDistvsRun   = (TProfile *)f->Get("hTrigTrkZDistvsRun");
  TProfile *hTrigTrkPhiDistvsRun = (TProfile *)f->Get("hTrigTrkPhiDistvsRun");
  TProfile *hTrigTrkNEtavsRun    = (TProfile *)f->Get("hTrigTrkNEtavsRun");
  TProfile *hTrigTrkNPhivsRun    = (TProfile *)f->Get("hTrigTrkNPhivsRun");
  TProfile *hNBemcEsvsRun        = (TProfile *)f->Get("hNBemcEsvsRun");
  TProfile *hMthTrkDeltaYvsRun     = (TProfile *)f->Get("hMthTrkDeltaYvsRun");
  TProfile *hMthTrkDeltaZvsRun     = (TProfile *)f->Get("hMthTrkDeltaZvsRun");
  TProfile *hMthTrkDeltaTOFvsRun     = (TProfile *)f->Get("hMthTrkDeltaTOFvsRun");
  
  
  Int_t binLow[nTrgIds],binHi[nTrgIds];
  memset(binLow,0,sizeof(binLow));
  memset(binHi,0,sizeof(binHi));
  for(Int_t i=0; i<nTrgIds; i++){
    IntMap::iterator iter = mRunId2RunIdx.find(runStart[i]);
    if(iter != mRunId2RunIdx.end()){
      binLow[i] = iter->second+1;
    }
    else{
      cout<<"Can not find the start runId for trgId["<<i<<"] in the run list !"<<endl;
      return;
    }
    
    iter = mRunId2RunIdx.find(runStop[i]);
    if(iter != mRunId2RunIdx.end()){
      binHi[i] = iter->second+1;
    }
    else{
      cout<<"Can not find the stop runId for trgId["<<i<<"] in the run list !"<<endl;
      return;
    }
  }
  for(Int_t i=0; i<nTrgIds; i++)
  cout<<"start run bin:"<<binLow[i]<<" \t stop run bin"<<binHi[i]<<endl;
  
  Int_t binLowT = binLow[0];
  Int_t binHiT = binHi[0];
  Int_t nDaySize;
  nDaySize = drawDay(runStartT, runStopT, runlist, binLowT, binHiT);
  TText *dayText = new TText(0.01,0.93,"Day");
  dayText->SetNDC();
  dayText->SetX(0.03);
  dayText->SetY(0.97);
  dayText->SetTextSize(0.025);
  dayText->SetY(0.95);
  
  //==============================================================================
  c1->cd();
  setHisto(hnEvtsvsRun, MStyle, MSize, MColor, MColor);
  hnEvtsvsRun->Draw("p");
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/hnEvtsvsRun.png", dir.Data()));
  
  if(dumpEmptyHT2Runs){
    ofstream outData("run11_zeroHT2EvtsInDD_runnumber.dat");
    outData<<"Run list with zero HT2 evts in Distributed Disk:"<<endl;
    Int_t nRuns = 0;
    for(Int_t i=grabRunIdx(runStart[0]); i<=grabRunIdx(runStop[nTrgIds-1]); i++){
      if(hnEvtsvsRun->GetBinContent(i+1)==0){
        outData<<grabRunId(i)<<endl;
      }
      else{
        nRuns++;
      }
    }
    outData.close();
    cout<<"# of none-zero statistic runs from HT2 events counter: "<<nRuns<<endl;
  }
  
  Int_t varIdx = 0;
  currentBadRuns.clear();
  currentEmptyRuns.clear();
  newBadRuns.clear();
  newEmptyRuns.clear();
  totalBadRuns.clear();
  totalEmptyRuns.clear();
  totalRejectRuns.clear();
  
  c1->cd();
  setProfile(hZdcXvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hZdcXvsRun->SetAxisRange(-1, 1, "Y");
  hZdcXvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hZdcXvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "zdcX", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/zdcXvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hBbcXvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hBbcXvsRun->SetAxisRange(-10, 60, "Y");
  hBbcXvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hBbcXvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "bbcX", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/bbcXvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hZdcXoverBbcXvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hZdcXoverBbcXvsRun->SetAxisRange(-0.1, 0.1, "Y");
  hZdcXoverBbcXvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hZdcXoverBbcXvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "zdcXOverBbcX", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/zdcXOverBbcXvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hTpcVxvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hTpcVxvsRun->SetAxisRange(-0.6, 0.6, "Y");
  hTpcVxvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hTpcVxvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "tpcVx", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/tpcVxvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hTpcVyvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hTpcVyvsRun->SetAxisRange(-0.8, 0.8, "Y");
  hTpcVyvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hTpcVyvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "tpcVy", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/tpcVyvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hTpcVzvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hTpcVzvsRun->SetAxisRange(-100, 100, "Y");
  hTpcVzvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hTpcVzvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "tpcVz", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/tpcVzvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hVpdVzvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hVpdVzvsRun->SetAxisRange(-500, 500, "Y");
  hVpdVzvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hVpdVzvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "vpdVz", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/vpdVzvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hVzDiffvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hVzDiffvsRun->SetAxisRange(-500, 500, "Y");
  hVzDiffvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hVzDiffvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "vzDiff", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/vzDiffvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hRefMultvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hRefMultvsRun->SetAxisRange(0, 100, "Y");
  hRefMultvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hRefMultvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "refMult", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/refMultvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hNTrksvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hNTrksvsRun->SetAxisRange(0, 250, "Y");
  hNTrksvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hNTrksvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "npTrks", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/npTrksvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hPtvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hPtvsRun->SetAxisRange(0.2, 0.8, "Y");
  hPtvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hPtvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "pt", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/ptvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hEtavsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hEtavsRun->SetAxisRange(-0.4, 0.4, "Y");
  hEtavsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hEtavsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "eta", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/etavsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hPhivsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hPhivsRun->SetAxisRange(-0.5, 0.5, "Y");
  hPhivsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hPhivsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "phi", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/phivsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hDcavsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hDcavsRun->SetAxisRange(0., 1.5, "Y");
  hDcavsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hDcavsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "dca", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dcavsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hNHitsFitvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hNHitsFitvsRun->SetAxisRange(40, 70, "Y");
  hNHitsFitvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hNHitsFitvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "nHitsFit", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsFitvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hNHitsPossvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hNHitsPossvsRun->SetAxisRange(50, 80, "Y");
  hNHitsPossvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hNHitsPossvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "nHitsPoss", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsPossvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hNHitsDedxvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hNHitsDedxvsRun->SetAxisRange(40, 70, "Y");
  hNHitsDedxvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hNHitsDedxvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "nHitsDedx", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nHitsDedxvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hDedxvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hDedxvsRun->SetAxisRange(2, 5, "Y");
  hDedxvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hDedxvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "dEdx", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/dEdxvsRun.png", dir.Data()));
  
  c1->cd();
  setProfile(hNSigmaEvsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hNSigmaEvsRun->SetAxisRange(-7, -1, "Y");
  hNSigmaEvsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hNSigmaEvsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "nSigmaE", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/nSigmaEvsRun.png", dir.Data()));
  
  
  c1->cd();
  setProfile(hBetavsRun, MStyle, MSize, MColor, MColor);
  if(setProfileRange) hBetavsRun->SetAxisRange(0.5, 5, "Y");
  hBetavsRun->Draw("p");
  for(Int_t i=0; i<nTrgIds; i++){
    addShade(hBetavsRun, binLow[i], binHi[i], i);
  }
  dumpBadRuns(badRunDir, "beta", varIdx);
  dayText->Draw("same");
  for (int l=0; l<nDaySize; l++){
    linesGlobal.at(l)->Draw("same");
    textsGlobal.at(l)->Draw("same"); //l+1?
  }
  pdfAction(c1, ps);
  c1->SaveAs(Form("%s/betavsRun.png", dir.Data()));
  
  //  c1->cd();
  //  setProfile(hNMthTrksvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hNMthTrksvsRun->SetAxisRange(0, 25, "Y");
  //  hNMthTrksvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hNMthTrksvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "nBemcMthTrks", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/nBemcMthTrksvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkPtvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkPtvsRun->SetAxisRange(0.6, 1.0, "Y");
  //  hMthTrkPtvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkPtvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkPt", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkPtvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkEtavsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkEtavsRun->SetAxisRange(-0.7, 0.7, "Y");
  //  hMthTrkEtavsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkEtavsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkEta", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkEtavsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkPhivsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkPhivsRun->SetAxisRange(-0.3, 0.6, "Y");
  //  hMthTrkPhivsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkPhivsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkPhi", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkPhivsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkNSigmaEvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkNSigmaEvsRun->SetAxisRange(-7, 0., "Y");
  //  hMthTrkNSigmaEvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkNSigmaEvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkNSigmaE", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkNSigmaEvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkAdc0vsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkAdc0vsRun->SetAxisRange(80, 240, "Y");
  //  hMthTrkAdc0vsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkAdc0vsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkAdc0", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkAdc0vsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkEvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkEvsRun->SetAxisRange(1.2, 2.2, "Y");
  //  hMthTrkEvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkEvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkE", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkEvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkZDistvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkZDistvsRun->SetAxisRange(-0.20, 20, "Y");
  //  hMthTrkZDistvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkZDistvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkZDist", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkZDistvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkPhiDistvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkPhiDistvsRun->SetAxisRange(-0.07, 0.02, "Y");
  //  hMthTrkPhiDistvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkPhiDistvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkPhiDist", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkPhiDistvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkNEtavsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkNEtavsRun->SetAxisRange(-1.2, 0.8, "Y");
  //  hMthTrkNEtavsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkNEtavsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkNEta", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkNEtavsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkNPhivsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkNPhivsRun->SetAxisRange(-1.2, 0.8, "Y");
  //  hMthTrkNPhivsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkNPhivsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "bemcMthTrkNPhi", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/bemcMthTrkNPhivsRun.png", dir.Data()));
  //  //===========================================================
  //  // MTD
  //  c1->cd();
  //  setProfile(hMthTrkDeltaYvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkDeltaYvsRun->SetAxisRange(-20, 20, "Y");
  //  hMthTrkDeltaYvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkDeltaYvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "mtdMthTrkDeltaY", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/mtdMthTrkDeltaYvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkDeltaZvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkDeltaZvsRun->SetAxisRange(-40, 60, "Y");
  //  hMthTrkDeltaZvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkDeltaZvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "mtdMthTrkDeltaZ", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/mtdMthTrkDeltaZvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hMthTrkDeltaTOFvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hMthTrkDeltaTOFvsRun->SetAxisRange(-1200, 0, "Y");
  //  hMthTrkDeltaTOFvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hMthTrkDeltaTOFvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "mtdMthTrkDeltaTOF", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/mtdMthTrkDeltaTOFvsRun.png", dir.Data()));
  
  //
  //  c1->cd();
  //  setProfile(hNTrigTrksvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hNTrigTrksvsRun->SetAxisRange(0, 0.1, "Y");
  //  hNTrigTrksvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hNTrigTrksvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "nTrigTrks", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/nTrigTrksvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkPtvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkPtvsRun->SetAxisRange(0, 30, "Y");
  //  hTrigTrkPtvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkPtvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkPt", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkPtvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkEtavsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkEtavsRun->SetAxisRange(-0.2, 0.2, "Y");
  //  hTrigTrkEtavsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkEtavsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkEta", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkEtavsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkPhivsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkPhivsRun->SetAxisRange(-0.6, 0.6, "Y");
  //  hTrigTrkPhivsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkPhivsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkPhi", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkPhivsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkNSigmaEvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkNSigmaEvsRun->SetAxisRange(-4, -2.5, "Y");
  //  hTrigTrkNSigmaEvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkNSigmaEvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkNSigmaE", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkNSigmaEvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkAdc0vsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkAdc0vsRun->SetAxisRange(340, 400, "Y");
  //  hTrigTrkAdc0vsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkAdc0vsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkAdc0", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkAdc0vsRun.png", dir.Data()));
  
  //c1->cd();
  //setProfile(hTrigTrkE0vsRun, MStyle, MSize, MColor, MColor);
  //if(setProfileRange) hTrigTrkE0vsRun->SetAxisRange(4, 8, "Y");
  //hTrigTrkE0vsRun->Draw("p");
  //for(Int_t i=0; i<nTrgIds; i++){
  //	addShade(hTrigTrkE0vsRun, binLow[i], binHi[i], i);
  //}
  //dumpBadRuns(badRunDir, "trigTrkE0", varIdx);
  //pdfAction(c1, ps);
  //c1->SaveAs(Form("%s/trigTrkE0vsRun.png", dir.Data()));
  
  //  c1->cd();
  //  setProfile(hTrigTrkEvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkEvsRun->SetAxisRange(6, 8, "Y");
  //  hTrigTrkEvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkEvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkE", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkEvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkZDistvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkZDistvsRun->SetAxisRange(-2, 2, "Y");
  //  hTrigTrkZDistvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkZDistvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkZDist", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkZDistvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkPhiDistvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkPhiDistvsRun->SetAxisRange(-0.01, 0.01, "Y");
  //  hTrigTrkPhiDistvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkPhiDistvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkPhiDist", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkPhiDistvsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkNEtavsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkNEtavsRun->SetAxisRange(-1.5, 2.5, "Y");
  //  hTrigTrkNEtavsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkNEtavsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkNEta", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkNEtavsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hTrigTrkNPhivsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hTrigTrkNPhivsRun->SetAxisRange(-1.5, 2.5, "Y");
  //  hTrigTrkNPhivsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hTrigTrkNPhivsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "trigTrkNPhi", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/trigTrkNPhivsRun.png", dir.Data()));
  //
  //  c1->cd();
  //  setProfile(hNBemcEsvsRun, MStyle, MSize, MColor, MColor);
  //  if(setProfileRange) hNBemcEsvsRun->SetAxisRange(0, 0.012, "Y");
  //  hNBemcEsvsRun->Draw("p");
  //  for(Int_t i=0; i<nTrgIds; i++){
  //    addShade(hNBemcEsvsRun, binLow[i], binHi[i], i);
  //  }
  //  dumpBadRuns(badRunDir, "nTrigEs", varIdx);
  //  pdfAction(c1, ps);
  //  c1->SaveAs(Form("%s/nTrigEsvsRun.png", dir.Data()));
  
  ps->On();
  ps->Close();
  
  cout<<"End of program !"<<endl;
  return;
}
//___________________________________________________________________
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
//___________________________________________________________________
Bool_t produceShade(TProfile *h, Int_t binLow, Int_t binHi, Double_t &mean, Double_t &rms, TGraphErrors *shade)
{
  mean=0.;
  rms=0.;
  
  Int_t nRuns = 0;
  for(Int_t i=binLow;i<=binHi;i++){
    Double_t binContent = h->GetBinContent(i);
    Double_t binError = h->GetBinError(i);
    if(TMath::Abs(binContent)<epsilon && TMath::Abs(binError)<epsilon) continue; // zero good variable event
    
    mean += h->GetBinContent(i);
    nRuns++;
  }
  mean = mean/nRuns;
  
  cout<<endl;
  //  cout<<"# of none-zero good variable statistic runs:"<<nRuns<< " mean: "<<mean<<" rms: "<<rms<<endl;
  
  for(Int_t i=binLow;i<=binHi;i++){
    Double_t binContent = h->GetBinContent(i);
    Double_t binError = h->GetBinError(i);
    if(TMath::Abs(binContent)<epsilon && TMath::Abs(binError)<epsilon) continue; // zero good variable event
    
    rms += pow(h->GetBinContent(i)-mean,2);
  }
  rms = sqrt(rms/nRuns);
  
  cout<<"# of none-zero good variable statistic runs:"<<nRuns<< " mean: "<<mean<<" rms: "<<rms<<endl;
  
  if(rms<epsilon){
    cout<<"the rms is zero !"<<endl;
    return kFALSE;
  }
  
  Int_t NPoints = binHi-binLow+1;
  for(Int_t i=0;i<NPoints;i++){
    shade->SetPoint(i, h->GetBinCenter(binLow+i), mean);
    shade->SetPointError(i, 0, nSigma*rms);
  }
  
  return kTRUE;
}
//___________________________________________________________________
void addShade(TProfile *h, Int_t binLow, Int_t binHi, Int_t trgIdx)
{
  if(trgIdx==0){
    currentBadRuns.clear();
    currentEmptyRuns.clear();
    newBadRuns.clear();
    newEmptyRuns.clear();
  }
  
  TGraphErrors *shade = new TGraphErrors(binHi-binLow+1);
  
  Double_t mean,rms;
  if( !produceShade(h,binLow,binHi,mean,rms,shade) ){
    shade->Delete();
    return;
  }
  shade->SetFillStyle(3003);
  shade->SetFillColor(kBlue-4);
  shade->DrawClone("e3same");
  shade->Delete();
  //cout<<mean<<"  "<<rms<<endl;
  
  Double_t empty_x[5000],    empty_y[5000];
  Double_t empty_xErr[5000], empty_yErr[5000];
  memset(empty_x, -1, sizeof(empty_x));
  memset(empty_y, -1, sizeof(empty_y));
  memset(empty_xErr, -1, sizeof(empty_xErr));
  memset(empty_yErr, -1, sizeof(empty_yErr));
  
  Double_t bad_x[5000],    bad_y[5000];
  Double_t bad_xErr[5000], bad_yErr[5000];
  memset(bad_x, -1, sizeof(bad_x));
  memset(bad_y, -1, sizeof(bad_y));
  memset(bad_xErr, -1, sizeof(bad_xErr));
  memset(bad_yErr, -1, sizeof(bad_yErr));
  
  Int_t nEmptyRuns = 0;
  Int_t nBadRuns   = 0;
  for(Int_t i=binLow;i<=binHi;i++){
    Int_t    nTrigEvts  = hnEvtsvsRun->GetBinContent(i);
    Double_t binCenter  = h->GetBinCenter(i);
    Double_t binWidth   = h->GetBinWidth(i);
    Double_t binContent = h->GetBinContent(i);
    Double_t binError   = h->GetBinError(i);
    
    if(nTrigEvts==0) continue; // zero HT2 triggered event
    
    if(TMath::Abs(binContent)<epsilon && TMath::Abs(binError)<epsilon){
      empty_x[nEmptyRuns]    = binCenter;
      empty_xErr[nEmptyRuns] = binWidth/2;
      empty_y[nEmptyRuns]    = binContent;
      empty_yErr[nEmptyRuns] = binError;
      nEmptyRuns++;
      
      currentEmptyRuns[i-1] = grabRunId(i-1);
      
      //cout<< setiosflags(ios::left) <<"runIdx: "<< setw(8) <<i-1<< "runId: "<< setw(12) <<grabRunId(i-1)<< "HT2 triggered events: "<< setw(10) <<nTrigEvts<<"bad_y(var): "<< setw(10) <<binContent<<"bad_yErr(varErr): "<< setw(10) <<binError<<endl;
    }
    else if(TMath::Abs(binContent-mean)>nSigma*rms && TMath::Abs(binContent-mean)>nSigma*binError){
      bad_x[nBadRuns]    = binCenter;
      bad_xErr[nBadRuns] = binWidth/2;
      bad_y[nBadRuns]    = binContent;
      bad_yErr[nBadRuns] = binError;
      nBadRuns++;
      
      currentBadRuns[i-1]  = grabRunId(i-1);
      
      //cout<< setiosflags(ios::left) <<"runIdx: "<< setw(8) <<i-1<< "runId: "<< setw(12) <<grabRunId(i-1)<< "HT2 triggered events: "<< setw(10) <<nTrigEvts<<"bad_y(var): "<< setw(10) <<binContent<<"bad_yErr(varErr): "<< setw(10) <<binError<<endl;
    }
  }
  
  TGraphErrors *grEmptyRuns = new TGraphErrors(nEmptyRuns, empty_x, empty_y, empty_xErr, empty_yErr);
  setGraph(grEmptyRuns, emptyMStyle, emptyMSize, emptyMColor, emptyMColor);
  if(nEmptyRuns>0) grEmptyRuns->DrawClone("pzsame");
  grEmptyRuns->Delete();
  
  TGraphErrors *grBadRuns = new TGraphErrors(nBadRuns, bad_x, bad_y, bad_xErr, bad_yErr);
  setGraph(grBadRuns, badMStyle, badMSize, badMColor, badMColor);
  if(nBadRuns>0) grBadRuns->DrawClone("pzsame");
  grBadRuns->Delete();
  
  for(IntMap::iterator curBadIter=currentBadRuns.begin();curBadIter!=currentBadRuns.end();curBadIter++){
    IntMap::iterator totBadIter = totalBadRuns.find(curBadIter->first);
    if(totBadIter == totalBadRuns.end()){
      newBadRuns[curBadIter->first]   = curBadIter->second;
      totalBadRuns[curBadIter->first] = curBadIter->second;
    }
    
    IntMap::iterator totRejIter = totalRejectRuns.find(curBadIter->first);
    if(totRejIter == totalRejectRuns.end()){
      totalRejectRuns[curBadIter->first] = curBadIter->second;
    }
  }
  
  for(IntMap::iterator curEmpIter=currentEmptyRuns.begin();curEmpIter!=currentEmptyRuns.end();curEmpIter++){
    IntMap::iterator totEmpIter = totalEmptyRuns.find(curEmpIter->first);
    if(totEmpIter == totalEmptyRuns.end()){
      newEmptyRuns[curEmpIter->first]   = curEmpIter->second;
      totalEmptyRuns[curEmpIter->first] = curEmpIter->second;
    }
    
    IntMap::iterator totRejIter = totalRejectRuns.find(curEmpIter->first);
    if(totRejIter == totalRejectRuns.end()){
      totalRejectRuns[curEmpIter->first] = curEmpIter->second;
    }
  }
}
//___________________________________________________________________
Bool_t dumpBadRuns(TString dirName, TString varName, Int_t &varIdx){
  ofstream outData(Form("%s/var%d_%s_badruns.dat", dirName.Data(), varIdx++, varName.Data()));
  if(!outData.is_open()){
    cout<<Form("Failed to dump bad runs for variable: %s", varName.Data())<<endl;
    return kFALSE;
  }
  
  outData<<"Bad runs (Non-zero good variable events) determined by current variable:"<<endl;
  for(IntMap::iterator curBadIter=currentBadRuns.begin();curBadIter!=currentBadRuns.end();curBadIter++){
    if(detailBadRunInfo) outData<< setiosflags(ios::left) <<"currentBadRuns - runIdx: "<< setw(8) <<curBadIter->first<<"runId: "<<curBadIter->second<<endl;
    else                 outData<<curBadIter->second<<endl;
  }
  outData<<endl;
  outData<<endl;
  outData<<"Empty runs (Zero good variable events) determined by current variable:"<<endl;
  for(IntMap::iterator curEmpIter=currentEmptyRuns.begin();curEmpIter!=currentEmptyRuns.end();curEmpIter++){
    if(detailBadRunInfo) outData<< setiosflags(ios::left) <<"currentEmptyRuns     - runIdx: "<< setw(8) <<curEmpIter->first<<"   runId: "<<curEmpIter->second<<endl;
    else                 outData<<curEmpIter->second<<endl;
  }
  outData<<endl;
  outData<<endl;
  outData<<"New bad runs determined by current variable:"<<endl;
  for(IntMap::iterator newBadIter=newBadRuns.begin();newBadIter!=newBadRuns.end();newBadIter++){
    if(detailBadRunInfo) outData<< setiosflags(ios::left) <<"newBadRuns     - runIdx: "<< setw(8) <<newBadIter->first<<"   runId: "<<newBadIter->second<<endl;
    else                 outData<<newBadIter->second<<endl;
  }
  outData<<endl;
  outData<<endl;
  outData<<"New empty runs determined by current variable:"<<endl;
  for(IntMap::iterator newEmpIter=newEmptyRuns.begin();newEmpIter!=newEmptyRuns.end();newEmpIter++){
    if(detailBadRunInfo) outData<< setiosflags(ios::left) <<"newEmptyRuns     - runIdx: "<< setw(8) <<newEmpIter->first<<"   runId: "<<newEmpIter->second<<endl;
    else                 outData<<newEmpIter->second<<endl;
  }
  outData<<endl;
  outData<<endl;
  outData<<"Total reject runs so far (totalBadRuns+totalEmptyRuns: won't count overlap runs between these two categories):"<<endl;
  for(IntMap::iterator totRejIter=totalRejectRuns.begin();totRejIter!=totalRejectRuns.end();totRejIter++){
    if(detailBadRunInfo) outData<< setiosflags(ios::left) <<"totalRejectRuns     - runIdx: "<< setw(8) <<totRejIter->first<<"   runId: "<<totRejIter->second<<endl;
    else                 outData<<totRejIter->second<<endl;
  }
  
  outData.close();
  return kTRUE;
}
//___________________________________________________________________
Int_t drawDay(Int_t runStartT, Int_t runStopT, TString runlist,Int_t binLow, Int_t binHi){
  
  vector<double> runNum;
  
  vector<double> runCount;
  vector<double> runCount_err;
  vector<double> newDay;
  vector<TString> dayN;
  
  int run_count = 0;
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
  for(int i=grabRunIdx(runStartT); i<=grabRunIdx(runStopT); i++){
    int runNumberI = grabRunId(i);
    //    //    if (runNumberI >= 20366000 && runNumberI < 21001001) runNumberI = 21001001;
    //    //    if (runNumberI >= 21246012 && runNumberI < 21251016) runNumberI = 21251016;
    TString runNumber(Form("%d",runNumberI));
    TSubString day = runNumber(2,3);
    TSubString run = runNumber(5,3);
    //    cout<<"day: "<<day<<"          run: "<<run<<endl;
    if(runNumberI<21000000) dayI = (runNumberI-20000000)/1000;
    else dayI = (runNumberI-21000000)/1000;
    if(runNumberI==runStartT){
      dayStart = day;
    }
    dayEnd = day;
    
    
    if ((dayI != dayPrev) || (runNumberI == runStartT)){
      newDay.push_back(run_count);
      dayN.push_back(day);
    }
    
    
    //    cout << "working on run " << runNumber << endl;
    runCount.push_back(run_count);
    runNum.push_back(runNumberI);
    runNum.push_back(runNumberI);
    
    
    run_count++;
    dayPrev = dayI;
  }
  int nRuns = runCount.size();
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
      
      lineGlobal->SetX1(0.09+0.71*newDay.at(l)/nCount);
      lineGlobal->SetX2(0.09+0.71*newDay.at(l)/nCount);
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
  return newDay.size();
}
//___________________________________________________________________
void pdfAction(TCanvas *c, TPDF *ps)
{
  ps->On();
  c->Update();
  c->cd();
  ps->NewPage();
  ps->Off();
};
