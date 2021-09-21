#include "plot.h"

void plotjpsi(
	TString inFile            	= "jpsi_test.root"
){

gStyle->SetOptStat(0);
// load tree
TChain *const tt_event = new TChain("T");

    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
 
if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}

Int_t bins1 =100;
TH2F *kinep  = new TH2F("kinep", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinem  = new TH2F("kinem", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinjpsi  = new TH2F("kinjpsi", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinel  = new TH2F("kinel", "", bins1, 0., 20, bins1, -20, 20);
TH2F *kinp  = new TH2F("kinp", "", bins1, 100., 300., bins1, -10, 20);

TCut cut = "Q2MClepton>1";

tt_event->Draw("epetaMC:eppMC>>kinep",cut,"goff");
tt_event->Draw("emetaMC:empMC>>kinem",cut,"goff");
tt_event->Draw("JpsietaMC:JpsipMC>>kinjpsi",cut,"goff");
tt_event->Draw("eletaMC:elpMC>>kinel",cut,"goff");
tt_event->Draw("petaMC:ppMC>>kinp",cut,"goff");

//========================================================================//
TCanvas *c1 = new TCanvas("c1","",0,0,1100,1000);
DrawGammaCanvasSettings( c1, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(kinep, "e^{-} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinem, "e^{+} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinjpsi, "J/#psi", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinel, "Electron", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinp, "Proton", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
c1->Divide(3,2);
c1->cd(1);
kinep->Draw("colz");
c1->cd(2);
kinem->Draw("colz");
c1->cd(3);
kinjpsi->Draw("colz");
c1->cd(4);
kinel->Draw("colz");
c1->cd(5);
kinp->Draw("colz");


TH2F *Q2xbjhadron  = new TH2F("Q2xbjhadron", "", bins1,0,1, bins1,0,100);
TH2F *Q2xbjlepton  = new TH2F("Q2xbjlepton", "", bins1,0,1, bins1,0,100);



tt_event->Draw("Q2MChadron:xbjMChadron>>Q2xbjhadron", "", "goff");
tt_event->Draw("Q2MClepton:xbjMClepton>>Q2xbjlepton", "", "goff");

TCanvas *c2 = new TCanvas("c2","",0,0,1000,500);
DrawGammaCanvasSettings( c2, 0.095, 0.01, 0.01, 0.105);
c2->Divide(2,1);
SetStyleHistoTH1ForGraphs(Q2xbjlepton, "Using only electrons information", "x_{bj}", "Q^{2} [GeV^{2}]", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(Q2xbjhadron, "Jacquet-Blondel method", "x_{bj}", "Q^{2} [GeV^{2}]", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
c2->cd(1);
Q2xbjhadron->Draw("colz");
c2->cd(2);
Q2xbjlepton->Draw("colz");

TH1F *Q2hadron  = new TH1F("Q2hadron", "", bins1,0,10);
TH1F *Q2lepton  = new TH1F("Q2lepton", "", bins1,0,10);

tt_event->Draw("Q2MChadron>>Q2hadron", "", "goff");
tt_event->Draw("Q2MClepton>>Q2lepton", "", "goff");

TCanvas *c3 = new TCanvas("c3","",0,0,500,500);
DrawGammaCanvasSettings( c3, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(Q2hadron, "Using only electrons information", "x_{bj}", "Q^{2} [GeV^{2}]",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(Q2lepton, "Jacquet-Blondel method", "x_{bj}", "Q^{2} [GeV^{2}]", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
Q2hadron->SetLineColor(kRed);
Q2hadron->Draw();
Q2lepton->Draw("same");

TCut cut2 = "Q2_track3_lepton > 1.";
TH1F *histJpsi_track3_M1  = new TH1F("histJpsi_track3_M1", "", 100,2,4.55);
TH1F *histJpsi_track3_M1_cut  = new TH1F("histJpsi_track3_M1_cut", "", 100,2,4.55);
tt_event->Draw("Jpsi_track3_M1>>histJpsi_track3_M1", "", "goff");
tt_event->Draw("Jpsi_track3_M1>>histJpsi_track3_M1_cut", cut2, "goff");

TCanvas *c4 = new TCanvas("c4","",0,0,500,500);
DrawGammaCanvasSettings( c4, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(histJpsi_track3_M1, "3 Tracks", "M(ee)", "",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
histJpsi_track3_M1->SetLineColor(kBlack);
histJpsi_track3_M1->Draw();
histJpsi_track3_M1_cut->SetLineColor(kRed);
histJpsi_track3_M1_cut->Draw("same");

}