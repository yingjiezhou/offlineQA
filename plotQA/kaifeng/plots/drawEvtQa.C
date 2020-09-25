//#include "rootlogon.h"
void addpdf(TPDF* pdf)
{
  pdf->On();
  pdf->NewPage();
  gPad->Update();
  pdf->Off(); 
}
void drawtitle(TPDF* pdf,TCanvas* c,string s){
  c->cd();
  c->Draw();
  // setPad(0.1,0.1,0.05,0.12);
  TLatex t;
  t.SetTextSize(0.05);
  t.DrawText(0.2,0.5,s.c_str());
  TLatex la;
  la.SetTextSize(0.035);
  la.DrawText(0.1,0.3,(new TDatime())->AsSQLString());
  la.DrawText(0.1,0.2,"by Kaifeng");
  pdf->On();
  pdf->NewPage();
  gPad->Update();
  pdf->Off();
}

void drawEvtQa(){
//  myStyle();
TFile* f = TFile::Open("test_out_.QA.root");
TCanvas* c = new TCanvas("c","c");
TPDF* pdf = new TPDF("QAplots.pdf");
pdf->Off();
drawtitle(pdf,c,"QA plots");

TH1F* hgDca = (TH1F*)f->Get("hgDca");;
TH1F* hVz = (TH1F*)f->Get("hVz");
TH1F* hVzVpdVz = (TH1F*)f->Get("hVzVpdVz");
TH2F* hnEvsEtavsVz = (TH2F*)f->Get("hnEvsEtavsVz");
TH2F* hnTofMulvsRef = (TH2F*)f->Get("hnTofMulvsRef");
TH2F* hnTofMatvsRef = (TH2F*)f->Get("hnTofMatvsRef");
TH1F* hnHitsFit = (TH1F*)f->Get("hnHitsFit");
TH2F* hinvBetavsP = (TH2F*)f->Get("hinvBetavsP");
TH2F* hdEdx = (TH2F*)f->Get("hdEdx");
TH1F* hEta = (TH1F*)f->Get("hEta");
TH1F* hPhi = (TH1F*)f->Get("hPhi");
// TH1F* hevt = (TH1F*)

//draw


hVz->Draw();
addpdf(pdf);
hVzVpdVz->Draw();
hgDca->Draw();

hnEvsEtavsVz->Draw("colz");
addpdf(pdf);
hnTofMulvsRef->Draw("colz");
TF1* fmulcut = new TF1("fmulcut","x*0.3405+48",0,1600);
fmulcut->SetLineColor(kRed);
fmulcut->Draw("same");
TLine* lmulcut = new TLine(250,20,1600,20);
lmulcut->Draw("same");
lmulcut->SetLineColor(kRed);
gPad->SetLogz(1);
addpdf(pdf);
hnTofMatvsRef->Draw("colz");
addpdf(pdf);
hinvBetavsP->Draw("colz");

addpdf(pdf);
hnHitsFit->Draw();
addpdf(pdf);
hdEdx->Draw("colz");
addpdf(pdf);
TProfile* pdEdx  = (TProfile*)f->Get("Dedx")->Clone("pdedx");
TH1F* hdEdxerr  = (TH1F*)pdEdx->ProjectionX("hdedxerr");
for (int i = 0;i<hdEdxerr->GetNbinsX();i++)
{
  hdEdxerr->SetBinContent(i+1,pdEdx->GetBinError(i+1)*sqrt(pdEdx->GetBinEntries(i+1)));
}
hdEdxerr->GetYaxis()->SetRangeUser(0.22,0.26);
hdEdxerr->Draw();
addpdf(pdf);
TProfile* ptof  = (TProfile*)f->Get("Tof")->Clone("ptof");
TH1F* htoferr  = (TH1F*)ptof->ProjectionX("htoferr");
for (int i = 0;i<htoferr->GetNbinsX();i++)
{
  htoferr->SetBinContent(i+1,ptof->GetBinError(i+1)*sqrt(ptof->GetBinEntries(i+1)));
}
htoferr->GetYaxis()->SetRangeUser(0.18,0.24);
htoferr->Draw();
addpdf(pdf);
hEta->Draw();
addpdf(pdf);
hPhi->Draw();
addpdf(pdf);
hgDca->Draw();
hgDca->GetXaxis()->SetTitle("gDCA (cm)");
gPad->SetLogy(1);
addpdf(pdf);



pdf->On();
pdf->Close();
}
// TH3F	hVxVyVz;1	VxVyVz
// TH1F	hVz;1	Vz
// TH1F	hVr;1	Vr
// TH1F	hVzVpdVz;1	Vz-VpdVz(cm)
// TH2F	hnEvsEtavsVz;1	nElectron
// TH2F	hnTofMulvsRef;1	hnTofMul vs Ref
// TH2F	hnTofMatvsRef;1	nTofmatch VS Refmult
// TH2F	hnTofHitvsRef;1	hnTofHit vs Ref
// TH1D	hevt;1	hevt
// TH1D	hevtcut;1	hevtcut
// TH1D	hevtbadcut;1	Events after remove bad run
// TH1F	hnHitsFit;1	nHitsFit
// TH1F	hgDca;1	gDca
// TH2F	hinvBetavsP;1	#frac{1}{#beta} vs p
// TH2F	hdEdx;1	dEdx vs p*charge
// TH1F	hpt;1	hpt
// TH1F	hEta;1	Eta
// TH1F	hPhi;1	Phi
