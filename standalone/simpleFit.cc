// A program to read in data from a file, do a fit and make a plot
// Glen Cowan, Ben Sowden
// RHUL Physics
// November 2013

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TROOT.h>
#include <TPostScript.h>
#include <TStyle.h>

using namespace std;

// Simple polynomial for test fit; modify as needed.

const int npar = 3;    // number of parameters.  Needs to be global.
TCanvas* canvas;       // global so main can clean up after running.

double fitFunc (double xArray[], double par[]){
  double x = xArray[0];
  double f = 0.;
  for (int i=0; i<npar; i++){
    f += par[i] * pow(x,i);
  }
  return f;
}

void simpleFit() {

// first make a canvas and a 2D histogram for the axes

  canvas = new TCanvas("canvas", "canvas", 10, 10, 500, 500);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameBorderMode(0);   // need this to turn off red hist frame!

  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();

  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.07);
  gPad->SetBottomMargin(0.17);

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);

  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.04);
  gStyle->SetTitleFont(42, "hxy");    // for histogram and axis title
  gStyle->SetLabelFont(42, "xyz");    // for axis labels (values)

  gStyle->SetTitleOffset(0.8, "h");        // what does this do?
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.99);

  gROOT->ForceStyle();

  // can make histogram or alternatively use the histograms automatically
  // connected to the TF1 or TGraph objects

  double xMin = 0.;
  double xMax = 10.;
  double yMin = 0.;
  double yMax = 6.;
  TH2F* axhist = new TH2F("axhist", "title", 10, xMin, xMax, 10, yMin, yMax);
  axhist->SetTitle("");
  axhist->SetXTitle("x");
  axhist->SetYTitle("y");

  const int maxCol = 10;
  const int maxRow = 500;

  double u[maxCol];
  double dataArray[maxCol][maxRow];

  // Read in data sets from file and insert into an array.
  // For file with n colums, data will go into dataArray[i][j]
  // where i = column number and j = rolw (line) number.
  // To access a 1-d array for a single column i use dataArray[i].

  TString fileName;
  cout << "Enter file name: ";
  cin >> fileName;

  ifstream inFile;
  inFile.open(fileName);
  if (inFile.fail()) {
    cout << "Couldn't open file!" << endl;
    exit(1);
  }

  bool readLine = true;
  int lineNum = 0;
  int ncol;
  while ( readLine ){

    TString line;
    stringstream ss;
    line.ReadLine(inFile);
    readLine = inFile.good();

    if ( readLine ) {

      TString firstChar = line(0,1);
      bool useLine = firstChar != "#";       // # for comment line

      if ( useLine ){

        int i = 0;
        stringstream ss;
        ss << line;              // put whole line into ss
        TString token;
        bool getToken = true;

        while ( getToken ) {
          ss >> token;           // extracts one token
          if ( token.Length() > 0 ) {
            u[i] = token.Atof();
            i++;
          }
          else {
            getToken = false;
          }
        }             // getToken
        ncol = i;

        for (int i=0; i<ncol; i++){
          dataArray[i][lineNum] = u[i];
        }
        lineNum++;

      }               // useLine

    }         // readLine

  }           // readLine

  int numPoints = lineNum;
  inFile.close();

  // Create TGraphErrors object for data

  double* x = dataArray[0];
  double* y = dataArray[1];
  double* dx = 0;       // null pointer means errors set to zero
  double* dy = dataArray[2];
  TGraphErrors* tge = new TGraphErrors(numPoints, x, y, dx, dy);

  // Create TF1 object and fit to the data

  double xmin = 0.;
  double xmax = 10.;
  TF1* f = new TF1("f", fitFunc, xmin, xmax, npar);
  TFitResultPtr r = tge->Fit(f, "S0");    // needed for covariances

  // Access fit results

  double chi2 = f->GetChisquare();
  int ndof = f->GetNDF();
  double pval = TMath::Prob(chi2, ndof);
  cout << "chi2 = " << chi2 << " for " << ndof
       << " degrees of freedom" << endl;
  cout << "p-value = " << pval << endl;
  cout << "par number, value, error..." << endl;
  double thetaHat[npar];
  double sigmaThetaHat[npar];
  for (int i=0; i<npar; i++){
    thetaHat[i] = f->GetParameter(i);
    sigmaThetaHat[i] = f->GetParError(i);
    cout << i << "  " << thetaHat[i] << " +- " << sigmaThetaHat[i] << endl;
  }
  cout << endl;

  cout << "Covariance matrix and its inverse..." << endl;
  TMatrixDSym V = r->GetCovarianceMatrix();
  TMatrixDSym Vinv(V);           // initialize Vinv with V
  Vinv.Invert();                 // writes inverse into Vinv
  for (int i=0; i<npar; i++){
    for (int j=0; j<npar; j++){
      cout << i << "  " << j << "  " << V[i][j] << "  " << Vinv[i][j] << endl;
     }
  }
  cout << endl;

  // Make plot

  TAxis* xa = axhist->GetXaxis();
  TAxis* ya = axhist->GetYaxis();

  xa->SetTitleOffset(1.4);    //  factor multiplies default offset
  ya->SetTitleOffset(1.4);

  xa->SetLabelOffset(0.015);
  ya->SetLabelOffset(0.015);

  xa->SetTickLength(0.015);  // default  = 0.03
  ya->SetTickLength(0.015);  // default  = 0.03

  xa->SetTitleSize(0.04);
  ya->SetTitleSize(0.04);

  xa->SetNdivisions(-5); // negative value should force number of divisions?
  ya->SetNdivisions(-6);

  xa->SetLabelSize(0.04);
  ya->SetLabelSize(0.04);

  // Draw axes and then add stuff

  // kDot=1, kPlus, kStar, kCircle=4, kMultiply=5,
  // kFullDotSmall=6, kFullDotMedium=7, kFullDotLarge=8,
  // kFullCircle=20, kFullSquare=21, kFullTriangleUp=22,
  // kFullTriangleDown=23, kOpenCircle=24, kOpenSquare=25,
  // kOpenTriangleUp=26, kOpenDiamond=27, kOpenCross=28,
  // kFullStar=29, kOpenStar=30

  axhist->Draw();

  tge->SetLineColor(1);
  tge->SetLineWidth(2);
  tge->SetLineStyle(1);
  tge->SetMarkerColor(1);
  tge->SetMarkerSize(0.8);
  tge->SetMarkerStyle(20);
  tge->Draw("PZ");              // Z turns off error bar tick marks

  f->SetRange(xmin, xmax);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  f->SetLineStyle(1);
  f->Draw("same");

  TLegend* leg = new TLegend(0.2, 0.7, 0.5, 0.9); // x1, y1, x2, y2
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(tge, "test data", "p");  // l for line, f for box, p for point
  leg->AddEntry(f, "fit result", "l");
  leg->Draw();

  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.04);
  tl->SetTextFont(42);
  tl->SetNDC();
  tl->DrawLatex(.65, 0.78, "A simple fit");

  TPostScript psfile("plot.eps", 113);     // 113 makes eps
  canvas->Draw();
  psfile.Close();
  canvas->Print("plot.gif", "gif");

}

int main(int argc, char* argv[]) {

  TApplication theApp("App", &argc, argv);
  simpleFit();
  theApp.Run(true);

  canvas->Close();
  delete canvas;

  return 0;
}


