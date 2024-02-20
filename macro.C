#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
#include <numeric>
#include <bits/stdc++.h>
#include <TMultiGraph.h>
#include <cmath>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TAttMarker.h>
#include <THStack.h>
#include <TLegend.h>
#include <string>
#include <ctime>

// ------------------------------------------------------------------------------------------------------------------------
// functions

int calo_track_corresponder (int calo_column, int track_layer){
  if (calo_column == 0 && track_layer >= 0 && track_layer <= 6) return 1;
  if (calo_column == 1 && track_layer >= 2 && track_layer <= 11) return 1;
  if (calo_column == 2 && track_layer >= 8 && track_layer <= 17) return 1;
  if (calo_column == 3 && track_layer >= 14 && track_layer <= 24) return 1;
  if (calo_column == 4 && track_layer >= 19 && track_layer <= 29) return 1;
  if (calo_column == 5 && track_layer >= 25 && track_layer <= 34) return 1;
  if (calo_column == 6 && track_layer >= 31 && track_layer <= 40) return 1;
  if (calo_column == 7 && track_layer >= 37 && track_layer <= 46) return 1;
  if (calo_column == 8 && track_layer >= 43 && track_layer <= 52) return 1;
  if (calo_column == 9 && track_layer >= 49 && track_layer <= 58) return 1;
  if (calo_column == 10 && track_layer >= 55 && track_layer <= 64) return 1;
  if (calo_column == 11 && track_layer >= 61 && track_layer <= 70) return 1;
  if (calo_column == 12 && track_layer >= 66 && track_layer <= 75) return 1;
  if (calo_column == 13 && track_layer >= 72 && track_layer <= 81) return 1;
  if (calo_column == 14 && track_layer >= 78 && track_layer <= 87) return 1;
  if (calo_column == 15 && track_layer >= 84 && track_layer <= 93) return 1;
  if (calo_column == 16 && track_layer >= 90 && track_layer <= 99) return 1;
  if (calo_column == 17 && track_layer >= 96 && track_layer <= 105) return 1;
  if (calo_column == 18 && track_layer >= 101 && track_layer <= 110) return 1;
  if (calo_column == 19 && track_layer >= 107 && track_layer <= 112) return 1;
  else return 0;
}

int get_cell (int bin_no, int fcr) {
  /*
  function to convert waveform bin to memory cell by shifting by the fcr
  */
  int cell = bin_no;
  if (cell + fcr < 1024){ //order bins by cell
    cell = cell + fcr;
  }
  else {
    cell = cell + fcr - 1024;
  }
  return cell;
} 


void make_hist(TH1D* graph, const char* title, const char* xtitle, const char* ytitle, const char* save_title){
  /*
  function to make and save a canvas with histogram plot
  */
  TCanvas* canv = new TCanvas(title, title, 1250, 900);
  gStyle->SetOptStat(0);
  graph->SetTitle("");
  graph->SetLabelSize(0.035, "xy");
  graph->SetTitleSize(0.035, "xy");
  graph->SetXTitle(xtitle);
  graph->SetYTitle(ytitle);
  graph->SetFillColor(598);
  graph->Draw("BAR");
  canv->SaveAs(save_title);
}

void make_tgraph(int n, double* x, double* y, int upper_limit, const char* title, const char* xtitle, const char* ytitle, const char* save_title, const char* display_title){
  /*
  function to make and save a canvas with TGraph plot
  */
  TGraph *graph = new TGraph(n, x, y);
  
  TCanvas* canv = new TCanvas(title, title, 1250, 900);
  graph->SetLineColor(9);
  graph->Draw();
  graph->GetXaxis()->SetLimits(0, upper_limit);
  //graph->SetTitle(display_title);
  graph->SetTitle("");
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetTitleSize(0.035);
  graph->GetYaxis()->SetTitleSize(0.035);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  canv->SaveAs(save_title);
}

double* slice(std::vector<double>& vec, int a, int b){
 
    // Starting and Ending iterators
    auto start = vec.begin() + a - 1;
    auto end = vec.begin() + b;
 
    // To store the sliced vector
    double* result = new double[b - a + 1];
 
    // Copy vector using copy function()
    std::copy(start, end, result);
 
    // Return the final sliced array
    return result;
}

void make_tgraph_with_removal(Int_t n, Double_t *x, Double_t *y, const char* title, const char* xtitle, const char* ytitle, const char* save_title){
  /*
  function to make and save a canvas with TGraph plot
  */

  //for splitting of y-array into two arrays - removing 0 values:
  std::vector<Double_t> vecy1;
  std::vector<Double_t> vecy2;
  std::vector<Double_t> vecx1;
  std::vector<Double_t> vecx2;
  int first_zero;

  for (int i = 0; i < n; i++){
    if (y[i] != 0){
      vecy1.push_back(y[i]);
      vecx1.push_back(x[i]);
    } else{
      first_zero = i;
      break;
    }
  }
  
  for (int j = first_zero+64; j < n; j++){
    vecy2.push_back(y[j]);
    vecx2.push_back(x[j]);
  }

  Double_t y1[vecy1.size()];
  Double_t y2[vecy2.size()];
  Double_t x1[vecx1.size()];
  Double_t x2[vecx2.size()];

  for (int i = 0; i < vecy1.size(); i++){
    y1[i] = vecy1.at(i);
    x1[i] = vecx1.at(i);
  }

  for (int i = 0; i < vecy2.size(); i++){
    y2[i] = vecy2.at(i);
    x2[i] = vecx2.at(i);
  }

  TCanvas* canv = new TCanvas(title, title, 1250, 900);

  gStyle->SetOptFit(100);

  TGraph *g1 = new TGraph(vecy1.size(), x1, y1);
  g1->SetLineColor(9);
  TGraph *g2 = new TGraph(vecy2.size(), x2, y2);
  g2->SetLineColor(9);
   
  auto graph = new TMultiGraph();
  graph->Add(g1);
  graph->Add(g2);

  graph->GetXaxis()->SetLimits(0, 1024);
  graph->SetMinimum(3538.5); // only checked for event 1, OM 124
  graph->SetMaximum(3557.5);
  graph->SetTitle("");
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetTitleSize(0.035);
  graph->GetYaxis()->SetTitleSize(0.035);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  graph->Draw("AL");
  canv->SaveAs(save_title);
}

void make_tgraph_single_wave(Int_t n, Double_t *x, Double_t *y, Double_t *xerr, Double_t *yerr, const char* title, const char* xtitle, const char* ytitle, const char *fitfunc, const char* save_title){
  /*
  function to make and save a canvas with TGraph plot
  */

  //for splitting of y-array into two arrays - removing 0 values:
  std::vector<Double_t> vecy1;
  std::vector<Double_t> vecy2;
  std::vector<Double_t> vecyerr1;
  std::vector<Double_t> vecyerr2;
  std::vector<Double_t> vecx1;
  std::vector<Double_t> vecx2;
  std::vector<Double_t> vecxerr1;
  std::vector<Double_t> vecxerr2;

  for (int i = 0; i<n; i++){
    if (y[i] != 0){
      vecy1.push_back(y[i]);
      vecyerr1.push_back(yerr[i]);
      vecx1.push_back(x[i]);
      vecxerr1.push_back(xerr[i]);
    }
    else {
      i=i+64;
      for (i; i<n; i++){
        vecy2.push_back(y[i]);
        vecyerr2.push_back(yerr[i]);
        vecx2.push_back(x[i]);
        vecxerr2.push_back(xerr[i]);
      }
      break;
    }
  }

  Double_t y1[vecy1.size()];
  Double_t y2[vecy2.size()];
  Double_t yerr1[vecy1.size()];
  Double_t yerr2[vecy2.size()];
  Double_t x1[vecx1.size()];
  Double_t x2[vecx2.size()];
  Double_t xerr1[vecx1.size()];
  Double_t xerr2[vecx2.size()];

  for (int i=0; i<vecy1.size(); i++){
    y1[i] = vecy1.at(i);
    yerr1[i] = vecyerr1.at(i);
    x1[i] = vecx1.at(i);
    xerr1[i] = vecxerr1.at(i);
  }

  for (int i=0; i<vecy2.size(); i++){
    y2[i] = vecy2.at(i);
    yerr2[i] = vecyerr2.at(i);
    x2[i] = vecx2.at(i);
    xerr2[i] = vecxerr2.at(i);
  }

  TCanvas* canv = new TCanvas(title, title, 1250, 900);

  //gStyle->SetOptFit(0000);

  TGraphErrors *g1 = new TGraphErrors(vecy1.size(), x1, y1, xerr1, yerr1);
  TGraphErrors *g2 = new TGraphErrors(vecy2.size(), x2, y2, xerr2, yerr2);
  g1->SetLineColor(9);
  g2->SetLineColor(9);
   
  auto graph = new TMultiGraph();
  graph->GetHistogram()->SetStats(0);
  graph->Add(g1);
  graph->Add(g2);
  
  graph->Fit(fitfunc);

  graph->GetXaxis()->SetLimits(0, 1024);
  graph->SetMinimum(3542);
  graph->SetMaximum(3554);
  graph->SetTitle("");
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetTitleSize(0.035);
  graph->GetYaxis()->SetTitleSize(0.035);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  graph->Draw("AP");
  canv->SaveAs(save_title);
}

void make_cal_tgraph_single_wave(Int_t n, Double_t *x, Double_t *y, Double_t *xerr, Double_t *yerr, Double_t *cal, const char* title, const char* xtitle, const char* ytitle, const char *fitfunc, const char* save_title){
  /*
  function to make and save a canvas with TGraph plot
  */

  //for splitting of y-array into two arrays - removing 0 values:
  std::vector<Double_t> vecy1;
  std::vector<Double_t> vecy2;
  std::vector<Double_t> vecyerr1;
  std::vector<Double_t> vecyerr2;
  std::vector<Double_t> vecx1;
  std::vector<Double_t> vecx2;
  std::vector<Double_t> vecxerr1;
  std::vector<Double_t> vecxerr2;
  std::vector<Double_t> veccal1;
  std::vector<Double_t> veccal2;

  for (int i = 0; i<n; i++){
    if (y[i] != 0){
      vecy1.push_back(y[i]);
      vecyerr1.push_back(yerr[i]);
      vecx1.push_back(x[i]);
      vecxerr1.push_back(xerr[i]);
      veccal1.push_back(cal[i]);
    }
    else {
      i=i+64;
      for (i; i<n; i++){
        vecy2.push_back(y[i]);
        vecyerr2.push_back(yerr[i]);
        vecx2.push_back(x[i]);
        vecxerr2.push_back(xerr[i]);
        veccal2.push_back(cal[i]);
      }
      break;
    }
  }

  Double_t y1[vecy1.size()];
  Double_t y2[vecy2.size()];
  Double_t yerr1[vecy1.size()];
  Double_t yerr2[vecy2.size()];
  Double_t x1[vecx1.size()];
  Double_t x2[vecx2.size()];
  Double_t xerr1[vecx1.size()];
  Double_t xerr2[vecx2.size()];

  for (int i=0; i<vecy1.size(); i++){
    y1[i] = vecy1.at(i) - veccal1.at(i);
    yerr1[i] = vecyerr1.at(i);
    x1[i] = vecx1.at(i);
    xerr1[i] = vecxerr1.at(i);
  }

  for (int i=0; i<vecy2.size(); i++){
    y2[i] = vecy2.at(i) - veccal2.at(i);
    yerr2[i] = vecyerr2.at(i);
    x2[i] = vecx2.at(i);
    xerr2[i] = vecxerr2.at(i);
  }

  TCanvas* canv = new TCanvas(title, title, 1250, 900);

  gStyle->SetOptFit(100);

  TGraphErrors *g1 = new TGraphErrors(vecy1.size(),x1,y1,xerr1,yerr1);
  TGraphErrors *g2 = new TGraphErrors(vecy2.size(),x2,y2,xerr2,yerr2);
  g1->SetLineColor(9);
  g2->SetLineColor(9);
   
  auto graph = new TMultiGraph();
  graph->Add(g1);
  graph->Add(g2);

  graph->GetXaxis()->SetLimits(0, 1024);
  graph->Fit(fitfunc);
  graph->SetTitle("");
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetTitleSize(0.035);
  graph->GetYaxis()->SetTitleSize(0.035);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  graph->Draw("AP");
  canv->SaveAs(save_title);
}

void calc_standard_dev(std::vector<std::vector<Double_t>> vec, Double_t *standard_dev_array, Double_t *calibration_const, const char *cal_check){

  //check whether to calibrate
  if (cal_check == "yes"){
    //calibrate
    for (int i = 0; i<vec.size(); i++) {
      int j = 0;
      for (auto it = vec.at(i).begin(); it != vec.at(i).end(); it++){
        *it = *it - calibration_const[j];
        j++;
      }
    }
    //removing null sections (with calibration)
    for (int i = 0; i<vec.size(); i++) {
      int j = 0;
      for (auto it = vec.at(i).begin(); it != vec.at(i).end(); it++){
        if (*it == (0 - calibration_const[j])){
          vec.at(i).erase(it);
          it--;
        }
        j++;
      }
    }
  }
  else{
    //removing null sections (without calibration)
    for (int i = 0; i<vec.size(); i++) {
      for (auto it = vec.at(i).begin(); it != vec.at(i).end(); it++){
        if (*it == 0){
          vec.at(i).erase(it);
          it--;
        }
      }
    }
  }

  //array of mean values
  Double_t mean_vals[vec.size()];
  for (int i=0; i<vec.size(); i++){
    Double_t sum_of_elems = 0;
    for(auto it = vec.at(i).begin(); it != vec.at(i).end(); it++){
      sum_of_elems = sum_of_elems + *it;
    }
    mean_vals[i] = sum_of_elems/960;
  }

  //array of std deviation
  Double_t standard_dev_from_vec[vec.size()];
  for (int i=0; i<vec.size(); i++){
    Double_t Nvariance = 0;
    for(auto it = vec.at(i).begin(); it != vec.at(i).end(); it++){
      Nvariance = Nvariance + ((*it - mean_vals[i]) * (*it - mean_vals[i]));
    }
    standard_dev_array[i] = (sqrt(Nvariance/960));
  }
}

void plot_standard_dev(Double_t *standard_dev_pre, Double_t *standard_dev_post, int length, int om_chosen, const char* title_pre, const char* title_post, int no_constants){
  std::string om_num_str = std::to_string(om_chosen);
  std::string no_constants_str = std::to_string(no_constants);
  std::string title1 = "calibration standard deviation values - OM: ";
  std::string title2 = " - No. calibration constants: ";
  std::string title0_pre = "Pre-";
  std::string title0_post = "Post-";
  std::string title_str_pre = title0_pre + title1 + om_num_str;
  std::string title_str_post = title0_post + title1 + om_num_str + title2 + no_constants_str;
  const char *title_char_pre = title_str_pre.c_str();
  const char *title_char_post = title_str_post.c_str();
  
  std::vector<double> vec_pre_cal;
  std::vector<double> vec_post_cal;

  for (int i=0; i<length; i++){
    vec_pre_cal.push_back(standard_dev_pre[i]);
    vec_post_cal.push_back(standard_dev_post[i]);
  }

  int nbins = 100; // redefine bins e.g. 100 bins between 0 nd 3
  
  std::sort (vec_pre_cal.begin(), vec_pre_cal.end());
  double xlow_pre = vec_pre_cal.front();
  double xup_pre = vec_pre_cal.back();
  std::sort (vec_post_cal.begin(), vec_post_cal.end());
  double xlow_post = vec_post_cal.front();
  double xup_post = vec_post_cal.back();

  TH1D *pre_cal = new TH1D(title_char_pre, title_char_pre, 10, xlow_pre, xup_pre);
    for (auto it = vec_pre_cal.begin(); it != vec_pre_cal.end(); it++) {
      pre_cal -> Fill(*it);
    }

  TCanvas* pre_cal_canv = new TCanvas(title_char_pre, title_char_pre, 900, 600);
    gStyle->SetOptStat("eMR");
    pre_cal->SetXTitle("Standard Deviation");
    pre_cal->SetYTitle("Frequency");
    pre_cal->Fit("gaus");
    pre_cal->Draw("E1");
    pre_cal_canv->SaveAs(title_pre);

  TH1D *post_cal = new TH1D(title_char_post, title_char_post, 10, xlow_post, xup_post);
    for (auto it = vec_post_cal.begin(); it != vec_post_cal.end(); it++) {
      post_cal -> Fill(*it);
    }

  TCanvas* post_cal_canv = new TCanvas(title_char_post, title_char_post, 900, 600);
    gStyle->SetOptStat("eMR");
    post_cal->SetXTitle("Standard Deviation");
    post_cal->SetYTitle("Frequency");
    post_cal->Fit("gaus");
    post_cal->Draw("E1");
    post_cal_canv->SaveAs(title_post);
}

void plot_combined_standard_dev(Double_t *standard_dev_pre, Double_t *standard_dev_post1, Double_t *standard_dev_post2, Double_t *standard_dev_post3, int length, int chosen_om){
  // removed variables:  float *pre_cal_mean_spread_array, float *post_cal_mean_spread_array, int initial_om
  
  std::string om_num_str = std::to_string(chosen_om);
  std::string title1 = "Standard deviation distributions - OM: ";
  std::string title_str =  title1 + om_num_str;
  std::string save_str0 = "calibration_plots/combined_dev_dist_OM_";
  std::string save_str_png = ".png";
  std::string save_str = save_str0 + om_num_str + save_str_png;
  const char *title_char = title_str.c_str();
  const char *save_char = save_str.c_str();
  
  std::vector<double> vec_pre_cal;
  std::vector<double> vec_post_cal1;
  std::vector<double> vec_post_cal2;
  std::vector<double> vec_post_cal3;

  for (int i = 0; i < length; i++){
    vec_pre_cal.push_back(standard_dev_pre[i]);
    vec_post_cal1.push_back(standard_dev_post1[i]);
    vec_post_cal2.push_back(standard_dev_post2[i]);
    vec_post_cal3.push_back(standard_dev_post3[i]);
  }

  int nbins = 65;
  
  std::sort (vec_pre_cal.begin(), vec_pre_cal.end());
  double xlow_pre = vec_pre_cal.front();
  double xup_pre = vec_pre_cal.back();
  std::sort (vec_post_cal1.begin(), vec_post_cal1.end());
  double xlow_post1 = vec_post_cal1.front();
  double xup_post1 = vec_post_cal1.back();
  std::sort (vec_post_cal2.begin(), vec_post_cal2.end());
  double xlow_post2 = vec_post_cal2.front();
  double xup_post2 = vec_post_cal2.back();
  std::sort (vec_post_cal3.begin(), vec_post_cal3.end());
  double xlow_post3 = vec_post_cal3.front();
  double xup_post3 = vec_post_cal3.back();

  TH1D *pre_cal = new TH1D("Pre-calibration", title_char, nbins, 1, 2.3);
  for (auto it = vec_pre_cal.begin(); it != vec_pre_cal.end(); it++){
    pre_cal->Fill(*it);
  }

  TH1D *post_cal1 = new TH1D("Post-calibration, 1024 constants", title_char, nbins, 1, 2.3);
  for (auto it = vec_post_cal1.begin(); it != vec_post_cal1.end(); it++){
    post_cal1->Fill(*it);
  }
    
  TH1D *post_cal3 = new TH1D("Post-calibration, 64 constants", title_char, nbins, 1, 2.3);
  for (auto it = vec_post_cal3.begin(); it != vec_post_cal3.end(); it++){
    post_cal3->Fill(*it);
  }
  
  TH1D *post_cal2 = new TH1D("Post-calibration, 16 constants", title_char, nbins, 1, 2.3);
  for (auto it = vec_post_cal2.begin(); it != vec_post_cal2.end(); it++){
    post_cal2->Fill(*it);
  }
  

  pre_cal->SetLineColor(kRed);
  pre_cal->SetTitle("");
  post_cal1->SetLineColor(kBlue);
  post_cal1->SetTitle("");
  post_cal2->SetLineColor(417);
  post_cal2->SetTitle("");
  post_cal3->SetLineColor(kViolet);
  post_cal2->SetTitle("");

  TCanvas *dev_canv = new TCanvas("dev_canv", "dev_canv", 1250, 900);
  
  // PRE CAL HIST
  
  pre_cal->Draw("E1");
  pre_cal->SetMarkerStyle(8);
  pre_cal->SetMarkerColor(kRed);
  pre_cal->Fit("gaus");
  pre_cal->SetMaximum(220);
  pre_cal->GetXaxis()->SetTitle("Standard deviation");
  pre_cal->GetYaxis()->SetTitle("Frequency");
  gStyle->SetOptStat("ne");
  TF1 *fit_pre = pre_cal->GetFunction("gaus");
  //pre_cal_mean_spread_array[chosen_om-initial_om] = fit_pre->GetParameter(1);
  fit_pre->SetLineColor(kBlack);

  dev_canv->Update();
  // Retrieve the stat box
  TPaveStats *ps_pre = (TPaveStats*)dev_canv->GetPrimitive("stats");
  ps_pre->SetName("mystats_pre");
  ps_pre->SetY1NDC(.8);
  ps_pre->SetY2NDC(.9);
  ps_pre->SetX1NDC(.65);
  ps_pre->SetX2NDC(.95);
  ps_pre->SetTextColor(kRed);
  TList *listOfLines_pre = ps_pre->GetListOfLines();

  // Remove the std dev line
  TText *tconst0_pre = ps_pre->GetLineWith("Std Dev");
  listOfLines_pre->Remove(tconst0_pre);
  
  TText *tconst1_pre = ps_pre->GetLineWith("Mean");
  listOfLines_pre->Remove(tconst1_pre);
  
  TText *tconst2_pre = ps_pre->GetLineWith("Entries");
  listOfLines_pre->Remove(tconst2_pre);
  
  TText *tconst3_pre = ps_pre->GetLineWith("chi");
  listOfLines_pre->Remove(tconst3_pre);
    
  // Add a new line in the stat box.
  // Note that "=" is a control character
  TLatex *myt_mean_pre = new TLatex(0,0, TString::Format("Mean: #mu = %.3f #pm %.3f", fit_pre->GetParameter(1), fit_pre->GetParError(1)));
  myt_mean_pre ->SetTextFont(42);
  myt_mean_pre ->SetTextColor(kRed);
  myt_mean_pre ->SetTextSize(0.03);

  TLatex *myt_chi_pre = new TLatex(0, 0.2, TString::Format("#chi^{2}#/Ndof = %.3f #/ %d", fit_pre->GetChisquare(), fit_pre->GetNDF()));
  myt_chi_pre ->SetTextFont(42);
  myt_chi_pre ->SetTextColor(kRed);
  myt_chi_pre ->SetTextSize(0.03);
  listOfLines_pre->Add(myt_mean_pre);
  //listOfLines_pre->Add(myt_chi_pre);
  ps_pre->SetTextSize(0.03);
  pre_cal->SetStats(0);
  
  // POST CAL HIST 1
  
  post_cal1->Draw("E1 sames");
  post_cal1->SetMarkerStyle(8);
  post_cal1->SetMarkerColor(kBlue);
  post_cal1->SetXTitle("Standard Deviation");
  post_cal1->SetYTitle("Frequency");
  post_cal1->Fit("gaus","", "sames");
  TF1 *fit_post1 = post_cal1->GetFunction("gaus");
  //post_cal_mean_spread_array[chosen_om-initial_om] = fit_post->GetParameter(1);
  fit_post1->SetLineColor(kBlack);
  
  dev_canv->Modified();
  dev_canv->Update();
  // Retrieve the stat box
  TPaveStats *ps_post1 = (TPaveStats*)dev_canv->GetPrimitive("stats");
  ps_post1->SetName("mystats_post1");
  ps_post1->SetY1NDC(.675);
  ps_post1->SetY2NDC(.775);
  ps_post1->SetX1NDC(.65);
  ps_post1->SetX2NDC(.95);
  ps_post1->SetTextColor(kBlue);
  TList *listOfLines_post1 = ps_post1->GetListOfLines();
  
  // Remove the std dev line
  TText *tconst0_post1 = ps_post1->GetLineWith("Std Dev");
  listOfLines_post1->Remove(tconst0_post1);
  
  TText *tconst1_post1 = ps_post1->GetLineWith("Mean");
  listOfLines_post1->Remove(tconst1_post1);
  
  TText *tconst2_post1 = ps_post1->GetLineWith("Entries");
  listOfLines_post1->Remove(tconst2_post1);
  
  TText *tconst3_post1 = ps_post1->GetLineWith("chi");
  listOfLines_post1->Remove(tconst3_post1);
  
  // Add a new line in the stat box.
  // Note that "=" is a control character
  TLatex *myt_mean_post1 = new TLatex(0,0, TString::Format("Mean: #mu = %.3f #pm %.3f", fit_post1->GetParameter(1), fit_post1->GetParError(1)));
  myt_mean_post1 ->SetTextFont(42);
  myt_mean_post1 ->SetTextColor(kBlue);
  myt_mean_post1 ->SetTextSize(0.03);

  TLatex *myt_chi_post1 = new TLatex(0, 0.2, TString::Format("#chi^{2}#/Ndof = %.3f #/ %d", fit_post1->GetChisquare(), fit_post1->GetNDF()));
  myt_chi_post1 ->SetTextFont(42);
  myt_chi_post1 ->SetTextColor(kBlue);
  myt_chi_post1 ->SetTextSize(0.03);
  listOfLines_post1->Add(myt_mean_post1);
  //listOfLines_post1->Add(myt_chi_post1);
  ps_post1->SetTextSize(0.03);
  post_cal1->SetStats(0);
  
  // POST CAL HIST 2
  
  post_cal2->Draw("E1 sames");
  post_cal2->SetMarkerStyle(8);
  post_cal2->SetMarkerColor(417);
  post_cal2->SetXTitle("Standard Deviation");
  post_cal2->SetYTitle("Frequency");
  gStyle->SetOptStat("ne");
  post_cal2->Fit("gaus","", "sames");
  TF1 *fit_post2 = post_cal2->GetFunction("gaus");
  //post_cal_mean_spread_array[chosen_om-initial_om] = fit_post->GetParameter(1);
  fit_post2->SetLineColor(kBlack);
  
  dev_canv->Modified();
  dev_canv->Update();
  // Retrieve the stat box
  TPaveStats *ps_post2 = (TPaveStats*)dev_canv->GetPrimitive("stats");
  ps_post2->SetName("mystats_post2");
  ps_post2->SetY1NDC(.55);
  ps_post2->SetY2NDC(.65);
  ps_post2->SetX1NDC(.65);
  ps_post2->SetX2NDC(.95);
  ps_post2->SetTextColor(417);
  TList *listOfLines_post2 = ps_post2->GetListOfLines();

  // Remove the std dev line
  TText *tconst0_post2 = ps_post2->GetLineWith("Std Dev");
  listOfLines_post2->Remove(tconst0_post2);
  
  TText *tconst1_post2 = ps_post2->GetLineWith("Mean");
  listOfLines_post2->Remove(tconst1_post2);
  
  TText *tconst2_post2 = ps_post2->GetLineWith("Entries");
  listOfLines_post2->Remove(tconst2_post2);
  
  TText *tconst3_post2 = ps_post2->GetLineWith("chi");
  listOfLines_post2->Remove(tconst3_post2);
    
  // Add a new line in the stat box.
  // Note that "=" is a control character
  TLatex *myt_mean_post2 = new TLatex(0,0, TString::Format("Mean: #mu = %.3f #pm %.3f", fit_post2->GetParameter(1), fit_post2->GetParError(1)));
  myt_mean_post2 ->SetTextFont(42);
  myt_mean_post2 ->SetTextColor(417);
  myt_mean_post2 ->SetTextSize(0.03);

  TLatex *myt_chi_post2 = new TLatex(0, 0.2, TString::Format("#chi^{2}#/Ndof = %.3f #/ %d", fit_post2->GetChisquare(), fit_post2->GetNDF()));
  myt_chi_post2 ->SetTextFont(42);
  myt_chi_post2 ->SetTextColor(417);
  myt_chi_post2 ->SetTextSize(0.03);
  listOfLines_post2->Add(myt_mean_post2);
  //listOfLines_post2->Add(myt_chi_post);
  ps_post2->SetTextSize(0.03);
  post_cal2->SetStats(0);
  
  // POST CAL HIST 3
  
  post_cal3->Draw("E1 sames");
  post_cal3->SetMarkerStyle(8);
  post_cal3->SetMarkerColor(kViolet);
  post_cal3->SetXTitle("Standard Deviation");
  post_cal3->SetYTitle("Frequency");
  gStyle->SetOptStat("ne");
  post_cal3->Fit("gaus","", "sames");
  TF1 *fit_post3 = post_cal3->GetFunction("gaus");
  //post_cal_mean_spread_array[chosen_om-initial_om] = fit_post->GetParameter(1);
  fit_post3->SetLineColor(kBlack);
  
  dev_canv->Modified();
  dev_canv->Update();
  // Retrieve the stat box
  TPaveStats *ps_post3 = (TPaveStats*)dev_canv->GetPrimitive("stats");
  ps_post3->SetName("mystats_post2");
  ps_post3->SetY1NDC(.425);
  ps_post3->SetY2NDC(.525);
  ps_post3->SetX1NDC(.65);
  ps_post3->SetX2NDC(.95);
  ps_post3->SetTextColor(kViolet);
  TList *listOfLines_post3 = ps_post3->GetListOfLines();

  // Remove the std dev line
  TText *tconst0_post3 = ps_post3->GetLineWith("Std Dev");
  listOfLines_post3->Remove(tconst0_post3);
  
  TText *tconst1_post3 = ps_post3->GetLineWith("Mean");
  listOfLines_post3->Remove(tconst1_post3);
  
  TText *tconst2_post3 = ps_post3->GetLineWith("Entries");
  listOfLines_post3->Remove(tconst2_post3);
  
  TText *tconst3_post3 = ps_post3->GetLineWith("chi");
  listOfLines_post3->Remove(tconst3_post3);
    
  // Add a new line in the stat box.
  // Note that "=" is a control character
  TLatex *myt_mean_post3 = new TLatex(0,0, TString::Format("Mean: #mu = %.3f #pm %.3f", fit_post3->GetParameter(1), fit_post3->GetParError(1)));
  myt_mean_post3 ->SetTextFont(42);
  myt_mean_post3 ->SetTextColor(kViolet);
  myt_mean_post3 ->SetTextSize(0.03);

  TLatex *myt_chi_post3 = new TLatex(0, 0.2, TString::Format("#chi^{2}#/Ndof = %.3f #/ %d", fit_post3->GetChisquare(), fit_post3->GetNDF()));
  myt_chi_post3 ->SetTextFont(42);
  myt_chi_post3 ->SetTextColor(kViolet);
  myt_chi_post3 ->SetTextSize(0.03);
  listOfLines_post3->Add(myt_mean_post3);
  //listOfLines_post3->Add(myt_chi_post);
  ps_post3->SetTextSize(0.03);
  post_cal3->SetStats(0);
  
  //Float_t rightmax = 1.1*post_cal->GetMaximum();
  //Float_t scale = dev_canv->GetUymax()/rightmax;
  //post_cal->Scale(scale);
  dev_canv->SaveAs(save_char);
}

std::string get_date(){

  time_t now = time(0);
  tm *ltm = localtime(&now);
  
  std::string day = std::to_string(ltm->tm_mday);
  std::string month = std::to_string(ltm->tm_mon + 1);
  std::string year = std::to_string(ltm->tm_year + 1900);
  
  std::string todays_date = day + "/" + month + "/" + year;
  
  return todays_date;
}

int snemo_run_time(int run_number)
{
  const char *cbd_base_path = "/sps/nemo/snemo/snemo_data/raw_data/CBD";

  std::vector<std::string> log_paths;
  log_paths.push_back(Form("%s/run-%d/snemo_trigger_run-%d.log", cbd_base_path, run_number, run_number));

  for (int crate=6; crate>=0; crate--)
    log_paths.push_back(Form("%s/run-%d/snemo_crate-%d_run-%d.log", cbd_base_path, run_number, crate, run_number));

  int unixtime = -1;

  for (const std::string & log_path : log_paths)
    {
      std::ifstream log_file (log_path);
      if (!log_file.is_open()) continue;

      std::string log_line;

      while (getline(log_file, log_line))
	{
	  // look if the line contains the wished prefix
	  size_t unixtime_index = log_line.find("run.run_unixtime_ms=");

	  if (unixtime_index == std::string::npos)
	    continue;

	  // extract the unix time
	  int an_unixtime = std::stoi(log_line.substr(20));

	  // keep the latest one
	  if (an_unixtime > unixtime)
	    unixtime = an_unixtime;
	}
    }

  return unixtime;
}

void write_csv(int run_number, Int_t run_date, std::vector<int> optical_module_numbers, std::vector<double> cells, std::vector<double> calibration_constants, std::vector<int> number_entries, std::vector<double> reduced_chi_squares, std::vector<double> standard_deviations, std::string filename) {

  //std::string date = get_date();
  int date = snemo_run_time(run_number);
  TDatime root_datime (date);
  std:: cout << "Today's date: " << date << std::endl;

  // create an ofstream for the file output
  std::fstream output_file;
  
  // create and open the .csv file
  output_file.open(filename, std::ios::out | std::ios::app);
    
  // write the file headers
  output_file << "Run number" << "," << "Date" << "," << "OM number" << "," << "Memory cell" << "," << "Calibration constant (ADC)" << ","  << "Number of entries" << "," << "Reduced chi square" << "," << "Standard deviation (ADC)" << std::endl;
  
  for (int cell = 0; cell <  1024; cell++) {  
    // write the data to the output file
    output_file << run_number << "," << root_datime.AsSQLString() << "," << optical_module_numbers.at(cell) << "," << cells.at(cell) << "," << calibration_constants.at(cell) << "," << number_entries.at(cell) << "," << reduced_chi_squares.at(cell) << "," << standard_deviations.at(cell) << std::endl;
  }
  
  output_file.close();
}

// ------------------------------------------------------------------------------------------------------------------------

int main() {
  TFile *file = new TFile("../data/snemo_run-1143_udd.root", "READ");
  std::vector<std::vector<short>> *wave = new std::vector<std::vector<short>>;
  std::vector<int> *calo_wall   = new std::vector<int>;
  std::vector<int> *calo_side   = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row    = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;
  std::vector<int> *calo_ampl   = new std::vector<int>;
  std::vector<int> *calo_type   = new std::vector<int>;
  std::vector<int> *fcr         = new std::vector<int>;
  std::vector<int> *tracker_side   = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_row    = new std::vector<int>;
  std::vector<int> *tracker_layer  = new std::vector<int>;

  int ncalo_main_wall = 520;
  int counter = 0;

  int eventnumber = 0;
  int calo_nohits = 0;
  int tracker_nohits = 0;
  Int_t date;

  TTree* tree = (TTree*)file->Get("SimData");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("header.eventnumber",1);
  tree->SetBranchAddress("header.eventnumber", &eventnumber);
  tree->SetBranchStatus("header.date", 1);
  tree->SetBranchAddress("header.date", &date);
  tree->SetBranchStatus("digicalo.nohits",1);
  tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
  tree->SetBranchStatus("digicalo.waveform",1);
  tree->SetBranchAddress("digicalo.waveform", &wave);
  tree->SetBranchStatus("digicalo.wall",1);
  tree->SetBranchAddress("digicalo.wall", &calo_wall);
  tree->SetBranchStatus("digicalo.side",1);
  tree->SetBranchAddress("digicalo.side", &calo_side);
  tree->SetBranchStatus("digicalo.column",1);
  tree->SetBranchAddress("digicalo.column", &calo_column);
  tree->SetBranchStatus("digicalo.row",1);
  tree->SetBranchAddress("digicalo.row", &calo_row);
  tree->SetBranchStatus("digicalo.charge",1);
  tree->SetBranchAddress("digicalo.charge", &calo_charge);
  tree->SetBranchStatus("digicalo.peakamplitude",1);
  tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
  tree->SetBranchStatus("digicalo.type", 1);
  tree->SetBranchAddress("digicalo.type", &calo_type);
  tree->SetBranchStatus("digicalo.fcr",1);
  tree->SetBranchAddress("digicalo.fcr", &fcr);
  tree->SetBranchStatus("digitracker.nohits",1);
  tree->SetBranchAddress("digitracker.nohits", &tracker_nohits);
  tree->SetBranchStatus("digitracker.side",1);
  tree->SetBranchAddress("digitracker.side", &tracker_side);
  tree->SetBranchStatus("digitracker.layer",1);
  tree->SetBranchAddress("digitracker.layer", &tracker_layer);
  tree->SetBranchStatus("digitracker.column",1);
  tree->SetBranchAddress("digitracker.column", &tracker_column);
  
  int single_waveform_event_no = 10001; // number of event to plot single waveform
  int no_cells = 1024;
  int max_entries = tree->GetEntries();
  std::cout << "Number of events: " << max_entries << "---------------------------------------------\n";
  
  int chosen_om = 124; // choose optical module
  std::vector<Double_t> tracking_event;
  std::vector<std::vector<double> > all_cells(no_cells); // memory cell * event vector, so we have the value for each event at each memory cell
  std::cout << "Date: " << date << std::endl;
  
  // ------------------------------------------------------------------------------------------------------------------------
  // build waveforms
  
  // setup waveform histograms to be filled
  std::string full_title_str = "Waveform summed over all events, OM: ";
  full_title_str += std::to_string(chosen_om);
  const char* full_title = full_title_str.c_str();
  std::string single_bin_title_str = "Single waveform as a function of bin, OM: ";
  single_bin_title_str += std::to_string(chosen_om);
  const char* single_bin_title = single_bin_title_str.c_str();
  std::string single_title_str = "Single waveform, OM: ";
  single_title_str += std::to_string(chosen_om);
  const char* single_title = single_title_str.c_str();
  std::string single_removed_title_str = "Single waveform with bins removed, OM: ";
  single_removed_title_str += std::to_string(chosen_om);
  const char* single_removed_title = single_removed_title_str.c_str();
  std::string av_title_str = "Average (mean) waveform over all events, OM: ";
  av_title_str += std::to_string(chosen_om);
  const char* av_title = av_title_str.c_str();
  
  TH1D *waveform = new TH1D(full_title, full_title, no_cells, 0, no_cells);
  TH1D *swaveform_bin = new TH1D(single_bin_title, single_bin_title, no_cells, 0, no_cells);
  TH1D *swaveform = new TH1D(single_title, single_title, no_cells, 0, no_cells);
  TH1D *swaveform_removed = new TH1D(single_removed_title, single_removed_title, no_cells, 0, no_cells);
  TH1D *avwaveform = new TH1D(av_title, av_title, no_cells, 0, no_cells);
  
  // setup arrays to plot waveforms as TGraphs
  Double_t waveformT[no_cells];
  for (int bin = 0; bin < no_cells; bin++){
    waveformT[bin] = 0;
  }
  Double_t swaveformT[no_cells];
  Double_t swaveform_binT[no_cells];
  Double_t swaveform_removedT[no_cells];
  Double_t cell_numsT[no_cells];
  Double_t timeT[no_cells];
  
  //for calibrating single waveform
  Double_t single_waveT[no_cells];
  std::vector<std::vector<Double_t>> cell_ordered_waveforms; // vector containing readings from all memory cells for each event

  for (int event = 0; event < max_entries; event++) {  // event: loop on event number
    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; k++) {  // k: loop on calo hit number
      
      int om_num;
      
      // get OM number based on calorimeter wall
      if (calo_type->at(k) == 0){
        om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);
      } else if (calo_type->at(k) == 1){
        om_num = 520 + calo_side->at(k)*64 + calo_wall->at(k)*32 + calo_column->at(k)*16 + calo_row->at(k);
      } else if (calo_type->at(k) == 2) {
        om_num = 648 + calo_side->at(k)*32 + calo_wall->at(k)*16 + calo_column->at(k);
      }
      
      //int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);

      if (om_num == chosen_om){ // select only events for chosen OM

        tracking_event.push_back(event);
        std::vector<double> temp_ADC_vec(no_cells); // temporary vector to append to cell_ordered_waveforms
        
        std::vector<short> wave_k = wave->at(k); // get the waveform corresponding to calo hit k
        int fcr_k = fcr->at(k); // get the fcr corresponding to k
        
        //std::cout << "event no: " << event << "\nOM no: " << om_num << "\nFCR: " << fcr_k << "\n------------" << std::endl;
        
        for (int bin = 0; bin < no_cells; bin++){ // bin: loop on waveform bin
          
          double wave_at_bin = wave_k.at(bin); // wave value corresponding to calo hit number and waveform bin
          int cell = get_cell(bin, fcr_k); // convert waveform bin to memory cell
          cell_numsT[cell] = cell+1;
          timeT[bin] = bin * (400.0/1024.0);
            
          waveform->Fill(bin, wave_at_bin); // fill sum waveform without fcr shift
          waveformT[bin] += wave_at_bin;
          
          if (event == single_waveform_event_no){
              swaveform->Fill(cell, wave_at_bin);
              swaveformT[cell] = wave_at_bin;
              swaveform_bin->Fill(bin, wave_at_bin);
              swaveform_binT[bin] = wave_at_bin;
          }
          
          if (bin >= 960){
            temp_ADC_vec.at(cell) = 0;
            
            if (event == single_waveform_event_no){
              swaveform_removedT[cell] = 0;
            }
            
          } else{
            temp_ADC_vec.at(cell) = wave_at_bin;
            all_cells[cell].push_back(wave_at_bin); // input the wave value into xth vector of all_events (xth memory cell)
            
            if (event == single_waveform_event_no){
              swaveform_removed->Fill(cell, wave_at_bin);
              swaveform_removedT[cell] = wave_at_bin;
              single_waveT[cell] = wave_at_bin;
            }
          }
        }
      cell_ordered_waveforms.push_back(temp_ADC_vec);
      }
    counter++;
    }
  }
  
  // make waveform histograms
  make_hist(swaveform, single_title, "Cell", "ADC", "waveform_histograms/single_waveform.png");
  make_hist(swaveform_bin, single_bin_title, "Cell", "ADC", "waveform_histograms/single_waveform_against_bin.png");
  make_hist(swaveform_removed, single_removed_title, "Cell", "ADC", "waveform_histograms/single_waveform_with_removed.png");
  make_hist(waveform, full_title, "Bin", "ADC", "waveform_histograms/waveform.png");
  
  // make waveform TGraphs
  make_tgraph(no_cells, timeT, waveformT, 400, "Summed waveform", "Time (ns)", "ADC", "waveform_histograms/waveform_tgraph.png", "");
  make_tgraph(no_cells, cell_numsT, swaveformT, 1024, "Single waveform", "Memory cell", "ADC", "waveform_histograms/swaveform_tgraph.png", "");
  make_tgraph(no_cells, timeT, swaveform_binT, 400, "Single waveform (bin)", "Time (ns)", "ADC", "waveform_histograms/swaveform_bin_tgraph.png", "");
  make_tgraph_with_removal(no_cells, cell_numsT, swaveform_removedT, "Single waveform (removed)", "Memory cell", "ADC", "waveform_histograms/swaveform_removed_tgraph.png");
  
  // ------------------------------------------------------------------------------------------------------------------------
  // find average wavefom and calculate calibration constant for each memory cell

  std::vector<double> averages; // vector to store the average value over all events for each memory cell
  double chi_sum = 0; // to find average chi square
  
  // define variables to plot against memory cell
  std::vector<double> means;
  std::vector<double> sigmas;
  std::vector<double> red_chi2s;
  std::vector<double> cell_nums;
  std::vector<int> no_events;
  
  std::vector<int> om_numbers;
  
  // also define as arrays
  Double_t empty_x_errT[1024];
  Double_t meansT[1024];
  Double_t sigmasT[1024];
  Double_t red_chi2sT[1024];
  Double_t no_eventsT[1024];
  
  for (int cell = 0; cell < no_cells; cell++){ // cell: loop on memory cell
    std::vector<double> ADC_vals = all_cells.at(cell); // select the memory cell from all_cells
    int no_vals = ADC_vals.size(); // number of ADC readings/entries
    
    double average = std::accumulate(ADC_vals.begin(), ADC_vals.end(), 0.0) / no_vals; // calculate average ADC across all events for this cell
    averages.push_back(average); // add average to averages vector
    avwaveform->Fill(cell, average);
    
    // plot a graph of ADC against event number to check for time dependence of the mean
    
    std::string temp2 = "ADC by Event number, OM: ";
    temp2 += std::to_string(chosen_om);
    temp2 += ",   Memory cell: ";
    temp2 += std::to_string(cell+1);
    const char* plot_title = temp2.c_str();

    std::string temp3 = "event_num_plots/ADCbyevent_cell_";
    temp3 += std::to_string(cell+1);
    temp3 += ".png";
    const char* plot_name = temp3.c_str();
    
    auto tprof = new TProfile(plot_title, plot_title, 100, 1, no_vals);
    
    for (int event = 0; event < no_vals; event++){
      tprof->Fill(event+1, ADC_vals.at(event));
    }
    
    TCanvas* canv_en = new TCanvas(plot_title, plot_title, 900, 600);
    //gStyle->SetOptStat(10);
    tprof->SetXTitle("Event number");
    tprof->SetYTitle("ADC");
    tprof->Draw();
    canv_en->SaveAs(plot_name);

    // create histogram title
    std::string temp = "OM: ";
    temp += std::to_string(chosen_om);
    temp += ",   Memory cell: ";
    temp += std::to_string(cell+1);
    const char* title = temp.c_str();
    
    // create histogram path
    std::string temp1 = "mem_cell_histograms/cell_";
    temp1 += std::to_string(cell+1);
    temp1 += ".png";
    const char* save_name = temp1.c_str();
    
    // create histogram
    std::sort (ADC_vals.begin(), ADC_vals.end());
    int min_val = ADC_vals.front();
    int max_val = ADC_vals.back();
    TH1D *hist = new TH1D(title, title, no_vals, min_val, max_val);
    for (int val = 0; val < no_vals; val++){
      hist -> Fill(ADC_vals[val]);
    }

    // save in a canvas
    TCanvas* canv_mc = new TCanvas(title, title, 1250, 900);
    gStyle->SetOptStat("eMR");
    //gStyle->SetOptFit();

    hist->SetTitle("");
    hist->SetTitleSize(0.035, "xy");
    hist->SetLabelSize(0.035, "xy");
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus", "Q");
    hist->Draw("E1");
    
    canv_mc->Update();

    // get fit parameters
    TF1 *fit = hist->GetFunction("gaus");
    double mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double chi2 = fit->GetChisquare();
    double ndof = fit->GetNDF();
    
    // Retrieve the stat box
    TPaveStats *ps = (TPaveStats*)canv_mc->GetPrimitive("stats");
    ps->SetName("mystats");
    TList *listOfLines = ps->GetListOfLines();

    // Remove unwanted lines
    TText *tconst0 = ps->GetLineWith("Mean");
    listOfLines->Remove(tconst0);
    TText *tconst1 = ps->GetLineWith("Std Dev");
    listOfLines->Remove(tconst1);
    TText *tconst2 = ps->GetLineWith("chi");
    listOfLines->Remove(tconst2);
    
    // Add a new line in the stat box, note that "=" is a control character
    TLatex *myt0 = new TLatex(0, 0, TString::Format("#mu = %.2f #pm %.2f", mean, fit->GetParError(1)));
    myt0->SetTextFont(42);
    myt0->SetTextSize(0.045);

    TLatex *myt1 = new TLatex(0, 0.2, TString::Format("#sigma = %.2f #pm %.2f", sigma, fit->GetParError(2)));
    myt1->SetTextFont(42);
    myt1->SetTextSize(0.045);
    
    TLatex *myt2 = new TLatex(0, 0.2, TString::Format("#chi^{2}#/ndof = %.2f #/ %d", chi2, fit->GetNDF()));
    myt2->SetTextFont(42);
    myt2->SetTextSize(0.045);
    
    listOfLines->Add(myt0);
    listOfLines->Add(myt1);
    listOfLines->Add(myt2);
  
    // The following line is needed to avoid that the automatic redrawing of stats
    hist->SetStats(0);

    canv_mc->Modified();
    canv_mc->SaveAs(save_name);
    
    // fill variables to plot against memory cell
    cell_nums.push_back(cell+1);
    means.push_back(mean);
    sigmas.push_back(sigma);
    red_chi2s.push_back(chi2/ndof);
    no_events.push_back(no_vals);
    
    om_numbers.push_back(chosen_om);
    
    empty_x_errT[cell] = 0;
    meansT[cell] = mean;
    sigmasT[cell] = sigma;
    red_chi2sT[cell] = chi2 / ndof;
    no_eventsT[cell] = no_vals;
    
    // add chi2 to chi_sum
    chi_sum += chi2 / ndof;
  }
  
  // find average chi square
  double av_chi = chi_sum / 1024;
  
  // plot graphs against memory cell
  make_tgraph(no_cells, cell_numsT, meansT, 1024, "Mean ADC readings", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC_tgraph.png", "");
  make_tgraph(no_cells, cell_numsT, sigmasT, 1024, "Standard deviation", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas_tgraph.png", "");
  make_tgraph(no_cells, cell_numsT, red_chi2sT, 1024, "Chi squares", "Memory cell", "Reduced chi square", "mem_cell_plots/chi2_tgraph.png", "");
  make_tgraph(no_cells, cell_numsT, no_eventsT, 1024, "Number of events", "Memory cell", "Number of events", "mem_cell_plots/no_events_tgraph.png", "");
  
  // make average waveform histogram
  make_hist(avwaveform, av_title, "Cell", "ADC", "waveform_histograms/average_waveform.png");
  
  // plot tgraphs for calibrated and uncalibrated single wave
  make_tgraph_single_wave(1024, cell_numsT, swaveform_removedT, empty_x_errT, sigmasT, "ADC by cell - single waveform", "Cell", "ADC", "pol1", "calibration_plots/singular_waveform_pre_cal.png");
  make_cal_tgraph_single_wave(1024, cell_numsT, swaveform_removedT, empty_x_errT, sigmasT, meansT, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal.png");
  
  // calulate standard deviations
  Double_t pre_cal_std_dev [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev, meansT, "no");

  Double_t post_cal_std_dev [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev, meansT, "yes");
  
  // plot standard deviations
  plot_standard_dev(pre_cal_std_dev, post_cal_std_dev, sizeof(post_cal_std_dev)/sizeof(post_cal_std_dev[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation.png", "calibration_plots/post_cal_standard_deviation.png", 1024);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // calculate calibration constant for every 64 memory cells
  
  std::vector<std::vector<double> > all_cells64(64); // vector to replace all_cells, with data from every 64 cells combined

  Double_t meansT_64[64];
  Double_t sigmasT_64[64];
  Double_t red_chi2sT_64[64];
  Double_t cell_numsT_64[64];
  Double_t no_eventsT_64[64];
  
  // combine data from every 64 cells
  for (int i = 0; i < 64; i++){
    for (int j = 0; j < 16; j++){
      all_cells64.at(i).insert(all_cells64.at(i).end(), all_cells.at((64*j)+i).begin(), all_cells.at((64*j)+i).end());
    }
  }
  
  for (int cell = 0; cell < 64; cell++){ // cell: loop on memory cell (but in multiples of 64)
    std::vector<double> ADC_vals = all_cells64.at(cell); // select set of 16 memory cells from all_cells_64
    int no_vals = ADC_vals.size(); // number of ADC readings/entries
    
    // create histogram title
    std::string temp = "OM: ";
    temp += std::to_string(chosen_om);
    temp += ",   Memory cells: 64n + ";
    temp += std::to_string(cell+1);
    const char* title = temp.c_str();
    
    // create histogram path
    std::string temp1 = "mem_cell_histograms/cell_64n_";
    temp1 += std::to_string(cell+1);
    temp1 += ".png";
    const char* save_name = temp1.c_str();
    
    // create histogram
    std::sort (ADC_vals.begin(), ADC_vals.end());
    int min_val = ADC_vals.front();
    int max_val = ADC_vals.back();
    TH1D *hist = new TH1D(title, title, 100, min_val, max_val);
    for (int val = 0; val < no_vals; val++){
      hist -> Fill(ADC_vals[val]);
    }

    // save in a canvas
    TCanvas* canv_mc = new TCanvas(title, title, 900, 600);
    gStyle->SetOptStat("eMR");
    gStyle->SetOptFit();

    hist->SetTitleSize(0.035, "xy");
    hist->SetLabelSize(0.035, "xy");
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus", "Q");
    hist->Draw("E1");
    
    canv_mc->SaveAs(save_name);

    // get fit parameters
    TF1 *fit = hist->GetFunction("gaus");
    double mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double chi2 = fit->GetChisquare();
    double ndof = fit->GetNDF();

    meansT_64[cell] = mean;
    sigmasT_64[cell] = sigma;
    red_chi2sT_64[cell] = chi2 / ndof;
    cell_numsT_64[cell] = cell + 1;
    no_eventsT_64[cell] = no_vals;
  }
  
  // plot graphs against memory cell
  make_tgraph(64, cell_numsT_64, meansT_64, 64, "Mean ADC readings overlayed", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC64_tgraph.png", "");
  make_tgraph(64, cell_numsT_64, sigmasT_64, 64, "Standard deviation overlayed", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas64_tgraph.png", "");
  make_tgraph(64, cell_numsT_64, red_chi2sT_64, 64, "Chi squares overlayed", "Memory cell", "Reduced chi square", "mem_cell_plots/chi264_tgraph.png", "");
  make_tgraph(64, cell_numsT_64, no_eventsT_64, 64, "Number of events overlayed", "Memory cell", "Number of events", "mem_cell_plots/no_events64_tgraph.png", "");
  
  // create a meansT_64 of size 1024
  double meansT_64_full[1024];
  double sigmasT_64_full[1024];
  for (int i = 0; i < 64; i++){
    for (int j = 0; j < 16; j++){
      meansT_64_full[i + 64*j] = meansT_64[i];
      sigmasT_64_full[i + 64*j] = sigmasT_64[i];
    }
  }
  
  // plot tgraph for calibrated single wave
  make_cal_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT_64_full, meansT_64_full, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal64.png");
  
  // calculate standard deviations
  Double_t pre_cal_std_dev64 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev64, meansT_64_full, "no");

  Double_t post_cal_std_dev64 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev64, meansT_64_full, "yes");

  // plot standard deviations
  plot_standard_dev(pre_cal_std_dev64, post_cal_std_dev64, sizeof(post_cal_std_dev64)/sizeof(post_cal_std_dev64[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation64.png", "calibration_plots/post_cal_standard_deviation64.png", 64);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // calculate calibration constant for every 16 memory cells
  
  std::vector<std::vector<double> > all_cells16(16); // vector to replace all_cells, with data from every 16 cells combined

  Double_t meansT_16[16];
  Double_t sigmasT_16[16];
  Double_t red_chi2sT_16[16];
  Double_t cell_numsT_16[16];
  Double_t no_eventsT_16[16];
  
  // combine data from every 16 cells
  for (int i = 0; i < 16; i++){
    for (int j = 0; j < 64; j++){
      all_cells16.at(i).insert(all_cells16.at(i).end(), all_cells.at((16*j)+i).begin(), all_cells.at((16*j)+i).end());
    }
  }
  
  for (int cell = 0; cell < 16; cell++){ // cell: loop on memory cell (but in multiples of 16)
    std::vector<double> ADC_vals = all_cells16.at(cell); // select set of 64 memory cells from all_cells_16
    int no_vals = ADC_vals.size(); // number of ADC readings/entries

    // create histogram title
    std::string temp = "OM: ";
    temp += std::to_string(chosen_om);
    temp += ",   Memory cells: 16n + ";
    temp += std::to_string(cell+1);
    const char* title = temp.c_str();
    
    // create histogram path
    std::string temp1 = "mem_cell_histograms/cell_16n_";
    temp1 += std::to_string(cell+1);
    temp1 += ".png";
    const char* save_name = temp1.c_str();
    
    // create histogram
    std::sort (ADC_vals.begin(), ADC_vals.end());
    int min_val = ADC_vals.front();
    int max_val = ADC_vals.back();
    TH1D *hist = new TH1D(title, title, no_vals, min_val, max_val);
    for (int val = 0; val < no_vals; val++){
      hist -> Fill(ADC_vals[val]);
    }

    // save in a canvas
    TCanvas* canv_mc = new TCanvas(title, title, 900, 600);
    gStyle->SetOptStat("eMR");
    gStyle->SetOptFit();

    hist->SetTitleSize(0.035, "xy");
    hist->SetLabelSize(0.035, "xy");
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus", "Q");
    hist->Draw("E1");
    
    canv_mc->SaveAs(save_name);

    // get fit parameters
    TF1 *fit = hist->GetFunction("gaus");
    double mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double chi2 = fit->GetChisquare();
    double ndof = fit->GetNDF();
    
    meansT_16[cell] = mean;
    sigmasT_16[cell] = sigma;
    red_chi2sT_16[cell] = chi2 / ndof;
    cell_numsT_16[cell] = cell + 1;
    no_eventsT_16[cell] = no_vals;
  }

  // plot graphs against memory cell
  make_tgraph(16, cell_numsT_16, meansT_16, 16, "Mean ADC readings overlayed", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC16_tgraph.png", "");
  make_tgraph(16, cell_numsT_16, sigmasT_16, 16, "Standard deviation overlayed", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas16_tgraph.png", "");
  make_tgraph(16, cell_numsT_16, red_chi2sT_16, 16, "Chi squares overlayed", "Memory cell", "Reduced chi square", "mem_cell_plots/chi216_tgraph.png", "");
  make_tgraph(16, cell_numsT_16, no_eventsT_16, 16, "Number of events overlayed", "Memory cell", "Number of events", "mem_cell_plots/no_events16_tgraph.png", "");
  
  // create a means_arr16 of size 1024
  double meansT_16_full[1024];
  double sigmasT_16_full[1024];
  for (int i = 0; i < 16; i++){
    for (int j = 0; j < 64; j++){
      meansT_16_full[i + 16*j] = meansT_16[i];
      sigmasT_16_full[i + 16*j] = sigmasT_16[i];
    }
  }
  
  // plot tgraph for calibrated single wave
  make_cal_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT_16_full, meansT_16_full, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal16.png");
  
  // calculate standard deviations
  Double_t pre_cal_std_dev16 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev16, meansT_16_full, "no");

  Double_t post_cal_std_dev16 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev16, meansT_16_full, "yes");

  // plot standard deviations
  plot_standard_dev(pre_cal_std_dev16, post_cal_std_dev16, sizeof(post_cal_std_dev16)/sizeof(post_cal_std_dev16[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation16.png", "calibration_plots/post_cal_standard_deviation16.png", 16);
  
  // plot standard deviations for 1024 and 16 together
  plot_combined_standard_dev(pre_cal_std_dev, post_cal_std_dev, post_cal_std_dev64, post_cal_std_dev16, sizeof(post_cal_std_dev16)/sizeof(post_cal_std_dev16[0]), chosen_om);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // overlayed plots
  
  TCanvas* canv_overlayed_64 = new TCanvas("Every 64 means overlaid", "Every 64 means overlaid", 1250, 900);
  
  double* means_64 = slice(means, 1, 64);
  double* cell_nums_64 = slice(cell_nums, 1, 64);
  
  TGraph *gr64 = new TGraph(64, cell_nums_64, means_64);
  gr64->Draw("AL");
  gr64->GetXaxis()->SetLimits(1, 64);
  gr64->SetTitle("");
  gr64->GetXaxis()->SetLabelSize(0.035);
  gr64->GetYaxis()->SetLabelSize(0.035);
  gr64->GetXaxis()->SetTitleSize(0.035);
  gr64->GetYaxis()->SetTitleSize(0.035);
  gr64->GetXaxis()->SetTitle("Memory cell");
  gr64->GetYaxis()->SetTitle("Mean ADC");
  
  for (int i = 1; i < 16; i++){
    double* means_slice = slice(means, 1+64*i, 64*(i+1));
    TGraph *gr = new TGraph(64, cell_nums_64, means_slice);
    gr->SetLineColor(i+1);
    gr->Draw("LSame");
  }
  
  canv_overlayed_64->SaveAs("mem_cell_plots/overlaid_means_every_64.png");
  
  TCanvas* canv_overlayed_16_legend = new TCanvas("Every 16 means overlaid with legend", "Every 16 means overlaid with legend", 1250, 900);
  
  double* means_16 = slice(means, 1, 16);
  double* cell_nums_16 = slice(cell_nums, 1, 16);
  int colours[3] = {632, 600, 417};
  
  TGraph *gr16_legend = new TGraph(16, cell_nums_16, means_16);
  auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  
  gr16_legend->SetLineColor(colours[0]);
  gr16_legend->SetLineWidth(2);
  gr16_legend->Draw("AL");
  gr16_legend->GetXaxis()->SetLimits(1, 16);
  gr16_legend->SetTitle("");
  gr16_legend->GetXaxis()->SetLabelSize(0.035);
  gr16_legend->GetYaxis()->SetLabelSize(0.035);
  gr16_legend->GetXaxis()->SetTitleSize(0.035);
  gr16_legend->GetYaxis()->SetTitleSize(0.035);
  gr16_legend->GetXaxis()->SetTitle("Memory cell");
  gr16_legend->GetYaxis()->SetTitle("Mean ADC");
  legend->AddEntry(gr16_legend, "cells 1 to 16", "l");

  for (int i = 1; i < 3; i++){
    double* means_slice = slice(means, 1+16*i, 16*(i+1));
    TGraph *gr = new TGraph(16, cell_nums_16, means_slice);
    gr->SetLineColor(colours[i]);
    gr->SetLineWidth(2);
    gr->Draw("LSame");
    
    std::string first_str = std::to_string(1+16*i);
    std::string last_str = std::to_string(16*(i+1));
    std::string legend_str = "cells ";
    legend_str += first_str;
    legend_str += " to ";
    legend_str += last_str;
    const char* legend_char = legend_str.c_str();
    legend->AddEntry(gr, legend_char, "l");
  }

  legend->SetTextSize(0.04);
  legend->Draw();
  
  canv_overlayed_16_legend->SaveAs("mem_cell_plots/overlaid_means_every_16_legend.png");
  
  TCanvas* canv_overlayed_16 = new TCanvas("Every 16 means overlaid", "Every 16 means overlaid", 1250, 900);
  
  TGraph *gr16 = new TGraph(16, cell_nums_16, means_16);
  gr16->Draw("AL");
  gr16->GetXaxis()->SetLimits(1, 16);
  gr16->SetTitle("");
  gr16->GetXaxis()->SetLabelSize(0.035);
  gr16->GetYaxis()->SetLabelSize(0.035);
  gr16->GetXaxis()->SetTitleSize(0.035);
  gr16->GetYaxis()->SetTitleSize(0.035);
  gr16->GetXaxis()->SetTitle("Memory cell");
  gr16->GetYaxis()->SetTitle("Mean ADC");
  
  for (int i = 1; i < 64; i++){
    double* means_slice = slice(means, 1+16*i, 16*(i+1));
    TGraph *gr = new TGraph(16, cell_nums_16, means_slice);
    gr->SetLineColor(i+1);
    gr->Draw("LSame");
  }
  
  canv_overlayed_16->SaveAs("mem_cell_plots/overlaid_means_every_16.png");
  
  // ------------------------------------------------------------------------------------------------------------------------

  std::cout << "Number of events: " << counter << std::endl;
  TFile *newfile = new TFile("output_from_macro.root", "RECREATE");
  newfile->cd();
  waveform->Write();
  swaveform->Write();
  newfile->Close();
  
  write_csv(1143, date, om_numbers, cell_nums, means, no_events, red_chi2s, sigmas, "calibration_output.csv");
  
  // make histograms of data in csv
  
  std::sort (means.begin(), means.end());
  int min_mean = means.front();
  int max_mean = means.back();
  TH1D *means_hist = new TH1D("histogram of mean baselines", "histogram of mean baselines", 30, min_mean, max_mean);
  std::sort (no_events.begin(), no_events.end());
  int min_no_events = no_events.front();
  int max_no_events = no_events.back();
  TH1D *no_events_hist = new TH1D("histogram of number of events", "histogram of number of events", 30, min_no_events, max_no_events);
  std::sort (red_chi2s.begin(), red_chi2s.end());
  int min_red_chi2 = red_chi2s.front();
  int max_red_chi2 = red_chi2s.back();
  TH1D *red_chi2s_hist = new TH1D("histogram of chi squares", "histogram of chi squares", 30, min_red_chi2, max_red_chi2);
  std::sort (sigmas.begin(), sigmas.end());
  int min_sigma = sigmas.front();
  int max_sigma = sigmas.back();
  TH1D *sigmas_hist = new TH1D("histogram of standard deviations", "histogram of standard deviations", 30, min_sigma, max_sigma);
  
  for (int cell = 0; cell < no_cells; cell++){
    means_hist->Fill(means.at(cell));
    no_events_hist->Fill(no_events.at(cell));
    red_chi2s_hist->Fill(red_chi2s.at(cell));
    sigmas_hist->Fill(sigmas.at(cell));
  }
  
  make_hist(means_hist, "histogram of mean baselines", "Calibration constant (ADC)", "Frequency", "calibration_plots/means_histogram.png");
  make_hist(no_events_hist, "histogram of number of events", "Number of events", "Frequency", "calibration_plots/no_events_histogram.png");
  make_hist(red_chi2s_hist, "histogram of chi squares", "Reduced chi square", "Frequency", "calibration_plots/reduced_chi_squares_histogram.png");
  make_hist(sigmas_hist, "histogram of standard deviations", "Standard deviation (ADC)", "Frequency", "calibration_plots/standard_deviations_histogram.png");
  
  std::cout << "End of macro.C" << std::endl;
  return 0;
}
