// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; charged multiplicity
// using 91 because it includes ROOT library to save output to root files
// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "fastjet/ClusterSequence.hh"
#include "Pythia8/Pythia.h"
#include <iostream>
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
/*
#include "/home/michal/root/root/include/ROOT/RConfig.hxx"
#include "/home/michal/root/root/include/TH2F.h"
#include "/home/michal/root/root/include/TTree.h"
#include "/home/michal/root/root/include/TString.h"
#include "/home/michal/root/root/include/TFile.h"
*/
using namespace std;
using namespace Pythia8;

TTree *tree;
TFile *outFile;

TH2F *hEtaPhi; 

Float_t mEta, mPhi, mCharge, mPx, mPy, mPz;
Int_t mEvt, mTrk, mMult;

bool ConnectOutput(const char* output, const char* treeName){

  outFile = new TFile(output, "RECREATE");
  if(!outFile)
    return false;

  outFile->cd();
  tree = new TTree(treeName, treeName);
  if(!tree)
    return false;

  tree->Branch("mEvt", &mEvt, "mEvt/I");
  tree->Branch("mTrk", &mTrk, "trackID/I");
  tree->Branch("mMultiplicity", &mMult, "mMult/I");
  tree->Branch("mPx", &mPx, "mPx/F");
  tree->Branch("mPy", &mPy, "mPy/F");
  tree->Branch("mPz", &mPz, "mPz/F");
  tree->Branch("mCharge", &mCharge, "mCharge/F");
  tree->Branch("mEta", &mEta, "mEta/F");
  tree->Branch("mPhi", &mPhi, "mPhi/F");


  hEtaPhi = new TH2F("hEtaPhi", "Correlation distribution of eta and phi; #eta [-];#varphi [-]", 100, -2, 2, 100, -4,4);
  if(!hEtaPhi)
    return false;

  return true;
}

void fillInfo(Event& event,int iEvt, int iTrk){
  
  mEvt = iEvt;
  mTrk = event[iTrk].id();
  mPx = event[iTrk].px();
  mPy = event[iTrk].py();
  mPz = event[iTrk].pz();  
  mCharge = event[iTrk].charge();
  mPhi = event[iTrk].phi();
  mEta = event[iTrk].eta();
  mMult = event.size();

  hEtaPhi->Fill(mEta, mPhi);

}


int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  /*
  if (argc != 2){
    cout << "Incorrect number of arguments. Required only one argument and that is the position of outputfile. Returning." << endl;
    return 1;
  }
  string outputName = argv[1];
  */
  
  //need to initialize TApplication only if i am going to be using root graphical features
  //TApplication theApp("hist", &argc, argv);

  const char* outputName = "/home/michal/Pythia/jetGenerator/jetGen.root";
  if(!ConnectOutput(outputName, "recTree")){
    cout << "Could not create a proper outfile in " << outputName << ". Returning." << endl;
    return 1;
  }


  Pythia pythia;

  pythia.readFile("jetGenerator.cmnd");

  //read-in important settings from the input file
  int nAbort = pythia.settings.mode("Main:timesAllowErrors");
  double eCM = pythia.settings.parm("Beams:eCM");
  int nEvent = pythia.settings.mode("Main:numberOfEvents");

  //initialize pythia
  if(!pythia.init()){
    return 1;
  }

  int iAbort = 0;
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()){ 
      ++iAbort;
      if (iAbort == nAbort)
        return 1;
      continue;
    }

    // Find number of all final charged particles and fill histogram.
    for (int iTrk = 0; iTrk < pythia.event.size(); ++iTrk){
      if (pythia.event[iTrk].isFinal() && pythia.event[iTrk].isCharged()){
        fillInfo(pythia.event, iEvent, iTrk);
        tree->Fill();
      } //if 
    } //track loop
  } // End of event loop. Statistics. Histogram. Done.


  //saving the output
  outFile->cd();
  tree->Write();
  hEtaPhi->Write();
  outFile->Close();

  cout << "Events have been generated and saved to file: " << outputName << endl;

  return 0;
}
