#include <TH3F.h>
#include <TCanvas.h>

void draw3DHistograms() {
    TCanvas *c1 = new TCanvas("c1", "3D Histogram Draw Options", 800, 600);
    TH3F *h3 = new TH3F("h3", "3D Histogram", 10, 0, 10, 10, 0, 10, 10, 0, 10);

    // Fill the histogram with some random data
    for (int i = 0; i < 1000; ++i) {
        h3->Fill(gRandom->Uniform(10), gRandom->Uniform(10), gRandom->Uniform(10));
    }

    c1->Divide(3, 2);  // Create a canvas with a 3x2 grid

    c1->cd(1);
    h3->Draw("LEGO");

    c1->cd(2);
    h3->Draw("LEGO1");

    c1->cd(3);
    h3->Draw("SURF");

    c1->cd(4);
    h3->Draw("SURF2");

    c1->cd(5);
    h3->Draw("COLZ");

    c1->cd(6);
    h3->Draw("BOX");

    c1->Update();
}