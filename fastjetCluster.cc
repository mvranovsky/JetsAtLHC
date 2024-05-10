

#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace fastjet;

const int nEvents = 100;

TTree* tree;
TFile* file;
//TCanvas* c;
TH2F* hist[nEvents];
TH1I *hNJetsPerEvent, *hNJetsPerEventCut;
TH2F* hRapPhiCorr;
TH1D *hJetPt, *hJetPtCut;

vector<TH1*> histograms1D;

//settings for plots
const int textFont = 42;
const Double_t labelSize = 0.05;


//initial setting for jet clusters
//possible algorithms: kt_algorithm, antikt_algorithm, cambridge_algorithm
double R = 0.4;
JetDefinition jet_def(kt_algorithm, R);
vector<PseudoJet> particles;
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

  hNJetsPerEvent = new TH1I("hNJetsPerEvent", "Number of jets per event", 50, 0, 50);
  histograms1D.push_back(hNJetsPerEvent);
  hNJetsPerEventCut = new TH1I("hNJetsPerEventCut", "Number of jets per event", 10, 0, 10);
  histograms1D.push_back(hNJetsPerEventCut);
  hJetPt = new TH1D("hJetPt", "hJetPt", 50, 0, 50);
  histograms1D.push_back(hJetPt);
  hJetPtCut = new TH1D("hJetPtCut", "hJetPtCut", 50, 0, 50);
  histograms1D.push_back(hJetPtCut);

  hRapPhiCorr = new TH2F("hRapPhiCorr", "hRapPhiCorr",40 , -4, 4, 50, 0,6.28);



  return true;
}



void Make(int EvtNum){

  if(particles.size() == 0)
    return;
  // run the clustering, extract the jets
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  vector<PseudoJet> jetsFiltered;
  for (int i = 0; i < jets.size(); ++i){
    hJetPt->Fill( jets[i].pt() );
    if(jets[i].pt() > minPt)
      jetsFiltered.push_back(jets[i]);
      hJetPtCut->Fill( jets[i].pt() );
  }
  hNJetsPerEvent->Fill( jets.size() );
  hNJetsPerEventCut->Fill( jetsFiltered.size() );

  // fill histogram with number of jets per event
  cout << "event: " << EvtNum << " jets originally: " << jets.size() << " after pT cut: " << jetsFiltered.size() << endl;


  // run over all the jets and save info about constituent particles
  vector<PseudoJet> constituents;
  for (int iJet = 0; iJet < jetsFiltered.size() ; ++iJet) {
    constituents = jetsFiltered[iJet].constituents();
    hist[EvtNum] = new TH2F(Form("hist%iJet%i", EvtNum, iJet), "histogram of rapidity and azimuthal angle",40 , -4, 4, 50, 0,6.28);
    for (int iTrk = 0; iTrk < constituents.size() ; ++iTrk) {
      hist[EvtNum]->Fill( constituents[iTrk].rap(), constituents[iTrk].phi() );
      hRapPhiCorr->Fill(constituents[iTrk].rap(), constituents[iTrk].phi() );
    }// tracks in jet
  }// jets

  // clear particles vector
  particles.clear();
}

void SetUpEvent(int EvtNum){
  // setup addresses and load variables from tree
  int nEntries = tree->GetEntries();
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

}//SetUpEvent

void hist1D(){
  cout << "here" << endl;
  for(int i; i < histograms1D.size(); ++i){
    TCanvas* c = new TCanvas("c", "c", 800, 800);
    gPad->SetMargin(0.14,0.07,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 

    if(histograms1D[i]->GetName() == "hNJetsPerEvent"){
      histograms1D[i]->GetXaxis()->SetTitle("N_{jets} [-]");
    } else if(histograms1D[i]->GetName() == "hNJetsPerEventCut"){
      histograms1D[i]->GetXaxis()->SetTitle("N_{jets}^{after Cut} [-]");
    } else if(histograms1D[i]->GetName() == "hJetPt"){
      histograms1D[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    } else if(histograms1D[i]->GetName() == "hJetPtCut"){
      histograms1D[i]->GetXaxis()->SetTitle("p_{T}^{after Cut} [GeV/c]");
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
    histograms1D[i]->GetXaxis()->SetTitleOffset(0.8);
    histograms1D[i]->GetYaxis()->SetTitleOffset(0.9);
    histograms1D[i]->Draw("same");

    TString namePos = "figures/" + TString(histograms1D[i]->GetName()) + ".pdf";
    cerr << "Saving as " << namePos << endl;
    c->SaveAs(namePos);
    c->Close();
  }

}

void RapPhiPlot(){
  // print all the histograms from different events to canvas with different colors
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.11); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  gPad->SetTickx();
  gPad->SetTicky(); 

  cout << "here" << endl;
  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
    hist[iEvt]->SetStats(0); // Here 'hist' is your histogram object, and '0' means false.
    hist[iEvt]->GetXaxis()->SetTitle("y [-]");
    hist[iEvt]->GetYaxis()->SetTitle("#varphi [-]");
    hist[iEvt]->GetXaxis()->SetTitleFont(textFont);
    hist[iEvt]->GetYaxis()->SetTitleFont(textFont);
    hist[iEvt]->GetXaxis()->SetLabelFont(textFont);
    hist[iEvt]->GetYaxis()->SetLabelFont(textFont);
    hist[iEvt]->GetXaxis()->SetLabelSize(labelSize);
    hist[iEvt]->GetYaxis()->SetLabelSize(labelSize);
    hist[iEvt]->GetXaxis()->SetTitleSize(labelSize);
    hist[iEvt]->GetYaxis()->SetTitleSize(labelSize);
    hist[iEvt]->GetXaxis()->SetTitleOffset(0.9);
    hist[iEvt]->GetYaxis()->SetTitleOffset(0.9);
    hist[iEvt]->SetMarkerColor(iEvt);
    hist[iEvt]->SetMarkerSize(1.5);
    hist[iEvt]->SetLineColor(iEvt);
    hist[iEvt]->Draw("same");
  }
  //hRapPhiCorr->Draw("same colz");
  cout << "down here" << endl;

  c->SaveAs("figures/hRapPhi.pdf");
  c->Clear();

}

void finish(){
  if(file){
    file->Close();
    delete file;
  }

  if(tree){
    delete file;
  }
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
  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
    SetUpEvent(iEvt);
    if(particles.size() == 0)
      continue;
    Make(iEvt);
  }

  // create a plot with number of jets per event
  cout << "Creating 1D histograms..." << endl;
  hist1D();

  // create a 2D plot of rapidity and phi 
  cout << "Creating a 2D correlation histogram of rapidity and azimuthal angle..." << endl;
  RapPhiPlot();

  cout << "Ending fastjet code. Goodbye..." << endl;

  //finish();

  return 0;
} 