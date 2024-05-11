#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <TStyle.h>
#include "TMath.h"
#include "TH3.h"
#include "TGraph2D.h"



using namespace std;
using namespace fastjet;

const int nEvents = 100;

int j = 0;
TTree* tree;
TFile* file;
//TCanvas* c;
TH2F* hist[nEvents];
TH1I *hNJetsPerEvent, *hNJetsPerEventCut;
TH2F* hRapPhiCorr;
TH1D *hJetPt, *hJetPtCut, *hDijetPhi, *hDijetPhi2;
TH3F *hRapPhiPtCorr[nEvents];
TGraph2D *h3DGraph;

vector<TH1*> histograms1D;

//settings for plots
const int textFont = 42;
const Double_t labelSize = 0.05;


//initial setting for jet clusters
//possible algorithms: kt_algorithm, antikt_algorithm, cambridge_algorithm
double R = 0.4;
JetDefinition jet_def(kt_algorithm, R);
double minPt = 10; //GeV

bool ConnectInput(const char* fileName, const char* treeName){
  //open .root file with tree that contains info about the events and load the tree
  file = new TFile(fileName);
  if(!file){
    cout << "Could not open input file. Returning." << endl;
    return false;
  }
  tree = (TTree*)file->Get(treeName);
  if(!tree){
    cout << "Could not get tree. Returning." << endl;
    return false;
  }

  hNJetsPerEvent = new TH1I("hNJetsPerEvent", "Number of jets per event", 20, 0, 200);
  histograms1D.push_back(hNJetsPerEvent);
  hNJetsPerEventCut = new TH1I("hNJetsPerEventCut", "Number of jets per event", 10, 0, 10);
  histograms1D.push_back(hNJetsPerEventCut);
  hJetPt = new TH1D("hJetPt", "hJetPt", 50, 0, 50);
  histograms1D.push_back(hJetPt);
  hJetPtCut = new TH1D("hJetPtCut", "hJetPtCut", 50, 0, 50);
  histograms1D.push_back(hJetPtCut);
  hDijetPhi = new TH1D("hDijetPhi", "hDijetPhi", 36,0,360);
  histograms1D.push_back(hDijetPhi);
  hDijetPhi2 = new TH1D("hDijetPhi2", "hDijetPhi2", 36, -90, 270);
  histograms1D.push_back(hDijetPhi2);

  h3DGraph = new TGraph2D();

  hRapPhiCorr = new TH2F("hRapPhiCorr", "hRapPhiCorr",5 , -4, 4, 5, 0,6.28);


  return true;
}



void Make(vector<PseudoJet> particles, int EvtNum){

  // run the clustering, extract the jets
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  vector<PseudoJet> jetsFiltered;
  for (int i = 0; i < jets.size(); ++i){
    hJetPt->Fill( jets[i].pt() );
    if(jets[i].pt() > minPt){
      jetsFiltered.push_back(jets[i]);
      hJetPtCut->Fill( jets[i].pt() );
    }
  }
  hNJetsPerEvent->Fill( jets.size() );
  hNJetsPerEventCut->Fill( jetsFiltered.size() );

  // fill histogram with number of jets per event
  //cout << "event: " << EvtNum << " jets originally: " << jets.size() << " after pT cut: " << jetsFiltered.size() << endl;

  if(jetsFiltered.size() == 0 )
    return;

  // run over all the jets and save info about constituent particles
  vector<PseudoJet> constituents;
  for (int iJet = 0; iJet < jetsFiltered.size() ; ++iJet) {
    constituents = jetsFiltered[iJet].constituents();
    hist[EvtNum] = new TH2F(Form("hist%iJet%i", EvtNum, iJet), "histogram of rapidity and azimuthal angle",8 , -4, 4, 36, -180,180);
    hRapPhiPtCorr[EvtNum] = new TH3F(Form("hRPP%iJet%i", EvtNum, iJet), "hRapPhiPtCorr", 32, -4,4, 36, -180, 180, 50, 0, 100 );
    for (int iTrk = 0; iTrk < constituents.size() ; ++iTrk) {
      hist[EvtNum]->Fill( constituents[iTrk].rap(), constituents[iTrk].phi()*180/TMath::Pi() - 180);
      //hRapPhiCorr->Fill(constituents[iTrk].rap(), constituents[iTrk].phi() );
      h3DGraph->SetPoint(constituents[iTrk].rap() , constituents[iTrk].phi() *180/TMath::Pi(), constituents[iTrk].pt() , j);
      ++j;
    }// tracks in jet
    hRapPhiPtCorr[EvtNum]->Fill( jetsFiltered[iJet].rap(), jetsFiltered[iJet].phi()*180/TMath::Pi() - 180, jetsFiltered[iJet].pt() );
    hDijetPhi2->Fill( (jetsFiltered[iJet].phi() - jetsFiltered[0].phi() )*180/TMath::Pi() );
  }// jets

  if(jetsFiltered.size() == 2){
    hDijetPhi->Fill( (jetsFiltered[0].phi() - jetsFiltered[1].phi() )*180/TMath::Pi() );
  }

}


vector<PseudoJet> SetUpEvent(int EvtNum){
  // setup addresses and load variables from tree
  int nEntries = tree->GetEntries();
  vector<PseudoJet> particles;
  Float_t track_Px, track_Py, track_Pz, track_E;
  Int_t track_Evt, track_mult;
  tree->SetBranchAddress("mEvt", &track_Evt);
  tree->SetBranchAddress("mEnergy", &track_E);
  tree->SetBranchAddress("mPx", &track_Px);
  tree->SetBranchAddress("mPy", &track_Py);
  tree->SetBranchAddress("mPz", &track_Pz);
  tree->SetBranchAddress("mMultiplicity", &track_mult);


  // go over all the tracks stored in TTree and save those which have the specified event number in particles vector
  for (int iTrk = 0; iTrk < nEntries; ++iTrk) {
    // if particle belongs to the certain event, save its info in particles
    tree->GetEntry(iTrk);
    if(track_Evt == EvtNum){
      particles.push_back( PseudoJet(track_Px, track_Py, track_Pz, track_E) );
    }//if
  }//track loop

  return particles;

}//SetUpEvent

void DrawText(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(0.027);
   textSTAR -> SetFillColor(0);
   textSTAR -> SetTextFont(62);
   textSTAR -> SetTextAlign(13);
   textSTAR->AddText("JetsAtLHC");
   textSTAR->AddText(" Pythia8 #sqrt{s}= 5 TeV");

   textSTAR -> Draw("same");
}


void hist1D(){
  TCanvas *c;
  for(int i = 0; i < histograms1D.size(); ++i){
    c = new TCanvas("c", "c", 800, 800);
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.14,0.07,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLeftMargin(0.15);

    if(i == 0){
      histograms1D[i]->GetXaxis()->SetTitle("N_{jets} [-]");
    } else if(i ==1){
      histograms1D[i]->GetXaxis()->SetTitle("N_{jets}^{after Cut} [-]");
    } else if(i ==2){
      histograms1D[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      c->SetLogy();
    } else if(i == 3){
      histograms1D[i]->GetXaxis()->SetTitle("p_{T}^{after Cut} [GeV/c]");
      c->SetLogy();

    }
    histograms1D[i]->GetYaxis()->SetTitle("counts");
    histograms1D[i]->GetXaxis()->SetTitleFont(textFont);
    histograms1D[i]->GetYaxis()->SetTitleFont(textFont);
    histograms1D[i]->GetXaxis()->SetLabelFont(textFont);
    histograms1D[i]->GetYaxis()->SetLabelFont(textFont);
    histograms1D[i]->GetXaxis()->SetLabelSize(labelSize);
    histograms1D[i]->GetYaxis()->SetLabelSize(labelSize);
    histograms1D[i]->GetXaxis()->SetTitleSize(labelSize);
    histograms1D[i]->GetYaxis()->SetTitleSize(labelSize);
    histograms1D[i]->GetXaxis()->SetTitleOffset(1.);
    histograms1D[i]->GetYaxis()->SetTitleOffset(1.5);
    histograms1D[i]->SetTitle("");
    histograms1D[i]->Draw("same");


    DrawText(0.65, 0.88, 0.85,0.93);

    TString namePos = "figures/" + TString(histograms1D[i]->GetName()) + ".pdf";
    //cerr << "Saving as " << namePos << endl;
    c->SaveAs(namePos);
    c->Close();
  }

}

void drawGraph() {
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.11); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  gPad->SetTickx();
  gPad->SetTicky(); 

  h3DGraph->GetXaxis()->SetTitle("y [-]");
  h3DGraph->GetYaxis()->SetTitle("#varphi [#circ]");
  h3DGraph->GetXaxis()->SetTitleFont(textFont);
  h3DGraph->GetYaxis()->SetTitleFont(textFont);
  h3DGraph->GetXaxis()->SetLabelFont(textFont);
  h3DGraph->GetYaxis()->SetLabelFont(textFont);
  h3DGraph->GetXaxis()->SetLabelSize(labelSize);
  h3DGraph->GetYaxis()->SetLabelSize(labelSize);
  h3DGraph->GetXaxis()->SetTitleSize(labelSize);
  h3DGraph->GetYaxis()->SetTitleSize(labelSize);
  h3DGraph->GetXaxis()->SetTitleOffset(0.9);
  h3DGraph->GetYaxis()->SetTitleOffset(1.3);
  h3DGraph->Draw("same LEGO");

  c->SaveAs("figures/hRapPhiPtGraph.pdf");
  c->Clear();

}

void RapPhiPlot(){
  // print all the histograms from different events to canvas with different colors
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.11); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  gPad->SetTickx();
  gPad->SetTicky(); 

  //cout << "here" << endl;
  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
    hist[iEvt]->SetStats(0); // Here 'hist' is your histogram object, and '0' means false.
    hist[iEvt]->GetXaxis()->SetTitle("y [-]");
    hist[iEvt]->GetYaxis()->SetTitle("#varphi [#circ]");
    hist[iEvt]->GetXaxis()->SetTitleFont(textFont);
    hist[iEvt]->GetYaxis()->SetTitleFont(textFont);
    hist[iEvt]->GetXaxis()->SetLabelFont(textFont);
    hist[iEvt]->GetYaxis()->SetLabelFont(textFont);
    hist[iEvt]->GetXaxis()->SetLabelSize(labelSize);
    hist[iEvt]->GetYaxis()->SetLabelSize(labelSize);
    hist[iEvt]->GetXaxis()->SetTitleSize(labelSize);
    hist[iEvt]->GetYaxis()->SetTitleSize(labelSize);
    hist[iEvt]->GetXaxis()->SetTitleOffset(0.9);
    hist[iEvt]->GetYaxis()->SetTitleOffset(1.3);
    hist[iEvt]->SetMarkerColor(iEvt);
    hist[iEvt]->SetMarkerStyle(iEvt);
    hist[iEvt]->SetMarkerSize(1.);
    hist[iEvt]->Draw("same");
  }
  //hRapPhiCorr->Draw("same");
  cout << "Saving rapidity and azimuthal angle correlation plot as: figures/hRapPhi.pdf" << endl;

  c->SaveAs("figures/hRapPhi.pdf");
  c->Clear();

}

void rapPhiPtPlot(){
// print all the histograms from different events to canvas with different colors
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.11); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  //gPad->SetTickx();
  //gPad->SetTicky(); 

  cout << "here" << endl;
  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
    hRapPhiPtCorr[iEvt]->SetStats(0); // Here 'hist' is your histogram object, and '0' means false.
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetTitle("y [-]");
    hRapPhiPtCorr[iEvt]->GetYaxis()->SetTitle("#varphi [#circ]");
    hRapPhiPtCorr[iEvt]->GetZaxis()->SetTitle("p_{T} [GeV]");
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetTitleFont(textFont);
    hRapPhiPtCorr[iEvt]->GetYaxis()->SetTitleFont(textFont);
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetLabelFont(textFont);
    hRapPhiPtCorr[iEvt]->GetZaxis()->SetLabelFont(textFont);
    hRapPhiPtCorr[iEvt]->GetZaxis()->SetTitleFont(textFont);
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetLabelFont(textFont);
    //hRapPhiPtCorr[iEvt]->GetXaxis()->SetLabelSize(labelSize);
    //hRapPhiPtCorr[iEvt]->GetYaxis()->SetLabelSize(labelSize);
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetTitleSize(labelSize);
    hRapPhiPtCorr[iEvt]->GetYaxis()->SetTitleSize(labelSize);
    //hRapPhiPtCorr[iEvt]->GetZaxis()->SetLabelSize(labelSize);
    hRapPhiPtCorr[iEvt]->GetZaxis()->SetTitleSize(labelSize);
    hRapPhiPtCorr[iEvt]->GetXaxis()->SetTitleOffset(1.4);
    hRapPhiPtCorr[iEvt]->GetYaxis()->SetTitleOffset(1.4);
    hRapPhiPtCorr[iEvt]->GetZaxis()->SetTitleOffset(1.4);
    hRapPhiPtCorr[iEvt]->SetMarkerColor(iEvt);
    hRapPhiPtCorr[iEvt]->SetMarkerStyle(20);
    hRapPhiPtCorr[iEvt]->SetMarkerSize(1.);
    hRapPhiPtCorr[iEvt]->Draw("same SURF");
  }
  //hRapPhiCorr->Draw("same");
  cout << "Saving 3D plot as: figures/hRapPhiPt.pdf" << endl;

  c->SaveAs("figures/hRapPhiPt.pdf");
  c->Clear();


}


void finish(){
  if(file){
    file->Close();
    delete file;
  }

  //if(tree){
    //delete file;
  //}
}


int main () {

  cout << "Starting fastjet analysis..." << endl;


  const char* fileName = "jetGen.root";
  const char* treeName = "recTree";
  cout << "Connecting input..." << endl;
  if ( !ConnectInput(fileName, treeName) )
    return 1;
  cout << "Input root file connected..." << endl;
  
  // print out some infos
  cout << "Clustering with " << jet_def.description() << "..." << endl;

  // for loop over all the events
  vector<PseudoJet> particles;
  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
    particles = SetUpEvent(iEvt);
    if(particles.size() == 0)
      continue;
    Make(particles, iEvt);
  }

  // create a plot with number of jets per event
  cout << "Creating 1D histograms..." << endl;
  hist1D();

  // create a 2D plot of rapidity and phi 
  cout << "Creating a 2D correlation histogram of rapidity and azimuthal angle..." << endl;
  RapPhiPlot();

  // create a 3D plot of rapidity, phi and pT
  cout << "Creating a 2D correlation histogram of rapidity and azimuthal angle..." << endl;
  rapPhiPtPlot();

  // create a 3D plot of rapidity, phi and pT
  cout << "Creating a 3D correlation graph of rapidity, azimuthal angle and pT..." << endl;
  drawGraph();

  cout << "Ending fastjet code. Goodbye..." << endl;

  finish();

  return 0;
} 