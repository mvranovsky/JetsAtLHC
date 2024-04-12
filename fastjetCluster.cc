

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
TCanvas* c;
TH2F* hist[nEvents];
TH1I* hNJetsPerEvent;

//settings for plots
const int textFont = 42;
const Double_t labelSize = 0.05;


//initial setting for jet clusters
//possible algorithms: kt_algorithm, antikt_algorithm, cambridge_algorithm
double R = 0.4;
JetDefinition jet_def(antikt_algorithm, R);
vector<PseudoJet> particles;


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

  hNJetsPerEvent = new TH1I("hNJetsPerEvent", "Number of jets per event", 20, 0, 20);



  return true;
}



void Make(int EvtNum){


  // run the clustering, extract the jets
  ClusterSequence cs(particles, jet_def);
  cout << "here" << endl; //problem here in defining ClusterSequence 
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  // fill histogram with number of jets per event
  hNJetsPerEvent->Fill( jets.size() );

  // run over all the jets and save info about constituent particles
  vector<PseudoJet> constituents;
  for (int iJet = 0; iJet < jets.size() ; ++iJet) {
    constituents = jets[iJet].constituents();
    hist[EvtNum] = new TH2F(Form("hist%i", EvtNum), "histogram of rapidity and azimuthal angle in event",100, -2, 2, 100, -3.14,3.14);
    for (int iTrk = 0; iTrk < constituents.size() ; ++iTrk) {
      hist[EvtNum]->Fill( constituents[iTrk].rap(), constituents[iTrk].phi() );
    }// tracks in jet
  }// jets

  // clear particles vector
  particles.clear();
}

void SetUpEvent(int EvtNum){
  // setup addresses and load variables from tree
  int nEntries = tree->GetEntries();
  Float_t track_Px[nEntries], track_Py[nEntries], track_Pz[nEntries], track_E[nEntries];
  Int_t track_Evt[nEntries];
  tree->SetBranchAddress("mEvt", track_Evt);
  tree->SetBranchAddress("mEnergy", track_E);
  tree->SetBranchAddress("mPx", track_Px);
  tree->SetBranchAddress("mPy", track_Py);
  tree->SetBranchAddress("mPz", track_Pz);

  // go over all the tracks stored in TTree and save those which have the specified event number in particles vector
  for (int iTrk = 0; iTrk < nEntries; ++iTrk) {
    // if particle belongs to the certain event, save its info in particles
    if(track_Evt[iTrk] == EvtNum){
      particles.push_back( PseudoJet(track_Px[iTrk], track_Py[iTrk], track_Pz[iTrk], track_E[iTrk] ) );
    }//if
  }//track loop
}//SetUpEvent

void NumOfJetsHist(){
  // save histogram with number of jets per event
  c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  gPad->SetTickx();
  gPad->SetTicky(); 


  hNJetsPerEvent->GetXaxis()->SetTitle("jets per event [-]");
  hNJetsPerEvent->GetYaxis()->SetTitle("counts");
  hNJetsPerEvent->GetXaxis()->SetTitleFont(textFont);
  hNJetsPerEvent->GetYaxis()->SetTitleFont(textFont);
  hNJetsPerEvent->GetXaxis()->SetLabelFont(textFont);
  hNJetsPerEvent->GetYaxis()->SetLabelFont(textFont);
  hNJetsPerEvent->GetXaxis()->SetLabelSize(labelSize);
  hNJetsPerEvent->GetYaxis()->SetLabelSize(labelSize);
  hNJetsPerEvent->GetXaxis()->SetTitleSize(labelSize);
  hNJetsPerEvent->GetYaxis()->SetTitleSize(labelSize);
  hNJetsPerEvent->GetXaxis()->SetTitleOffset(0.8);
  hNJetsPerEvent->GetYaxis()->SetTitleOffset(0.9);
  hNJetsPerEvent->Draw();


  c->SaveAs("figures/hNJets.pdf");
  c->Close();

}

void RapPhiPlot(){
  // print all the histograms from different events to canvas with different colors
  c = new TCanvas("c", "c", 800, 800);
  gPad->SetMargin(0.14,0.07,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
  gPad->SetTickx();
  gPad->SetTicky(); 


  for (int iEvt = 0; iEvt < nEvents; ++iEvt) {
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
    hist[iEvt]->GetXaxis()->SetTitleOffset(0.8);
    hist[iEvt]->GetYaxis()->SetTitleOffset(0.9);
    hist[iEvt]->SetMarkerColor(iEvt);
    hist[iEvt]->SetLineColor(iEvt);
    hist[iEvt]->Draw("same colz");
  }

  c->SaveAs("figures/hRapPhi.pdf");
  c->Close();
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
  cout << "Creating histogram with number of jets per event..." << endl;
  NumOfJetsHist();

  // create a 2D plot of rapidity and phi 
  cout << "Creating a 2D correlation histogram of rapidity and azimuthal angle..." << endl;
  RapPhiPlot();

  cout << "Ending fastjet code. Goodbye..." << endl;

  return 0;
} 