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
  TCanvas* canv = new TCanvas(title, title, 900, 600);
  gStyle->SetOptStat(10);
  graph->SetXTitle(xtitle);
  graph->SetYTitle(ytitle);
  graph->Draw("hist l");
  canv->SaveAs(save_title);
}

void make_tgraph(int n, double* x, double* y, int upper_limit, const char* title, const char* xtitle, const char* ytitle, const char* save_title){
  /*
  function to make and save a canvas with TGraph plot
  */
  TGraph *graph = new TGraph(n, x, y);
  
  TCanvas* canv = new TCanvas(title, title, 900, 600);
  graph->SetLineColor(9);
  graph->Draw();
  graph->GetXaxis()->SetLimits(0, upper_limit);
  graph->SetTitle(title);
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

  TCanvas* canv = new TCanvas(title, title, 900, 600);

  gStyle->SetOptFit(100);

  TGraphErrors *g1 = new TGraphErrors(vecy1.size(),x1,y1,xerr1,yerr1);
  TGraphErrors *g2 = new TGraphErrors(vecy2.size(),x2,y2,xerr2,yerr2);
   
  auto graph = new TMultiGraph();
  graph->Add(g1);
  graph->Add(g2);

  graph->GetXaxis()->SetLimits(0, 1024);
  graph->Fit(fitfunc);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  graph->Draw("AP");
  canv->SaveAs(save_title);

  std::cout << "vecy1 size:" << vecy1.size() << " vecy2 size:" << vecy2.size();
}

void make_cal_tgraph_single_wave(Int_t n, Double_t *x, Double_t *y, Double_t *xerr, Double_t *yerr, const char* title, const char* xtitle, const char* ytitle, const char *fitfunc, const char* save_title, Double_t *cal){
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

  TCanvas* canv = new TCanvas(title, title, 900, 600);

  gStyle->SetOptFit(100);

  TGraphErrors *g1 = new TGraphErrors(vecy1.size(),x1,y1,xerr1,yerr1);
  TGraphErrors *g2 = new TGraphErrors(vecy2.size(),x2,y2,xerr2,yerr2);
   
  auto graph = new TMultiGraph();
  graph->Add(g1);
  graph->Add(g2);

  graph->GetXaxis()->SetLimits(0, 1024);
  graph->Fit(fitfunc);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  graph->Draw("AP");
  canv->SaveAs(save_title);

  std::cout << "vecy1 size:" << vecy1.size() << " vecy2 size:" << vecy2.size();
}

void calc_standard_dev (std::vector<std::vector<Double_t>> vec, Double_t *standard_dev_array, Double_t *calibration_const, const char *cal_check){

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

  TTree* tree = (TTree*)file->Get("SimData");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("header.eventnumber",1);
  tree->SetBranchAddress("header.eventnumber", &eventnumber);
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
  
  int chosen_om = 124; // choose optical module
  int single_waveform_event_no = 0; // number of event to plot single waveform
  int no_cells = 1024;
  std::vector<Double_t> tracking_event;
  
  std::vector<std::vector<double> > all_cells(no_cells); // memory cell * event vector, so we have the value for each event at each memory cell
  int max_entries = tree->GetEntries();
  
  // ------------------------------------------------------------------------------------------------------------------------
  // build waveforms
  
  TH1D *waveform = new TH1D("waveform", "waveform", no_cells, 0, no_cells);
  TH1D *swaveform = new TH1D("swaveform", "singular waveform", no_cells, 0, no_cells);
  TH1D *avwaveform = new TH1D("avwaveform", "average waveform", no_cells, 0, no_cells);
  
  //for calibrating single waveform
  Double_t single_waveT[1024];
  std::vector<std::vector<Double_t>> cell_ordered_waveforms;

  for (int event = 0; event < max_entries; event++) {  // i: loop on event number
    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; k++) {  // k: loop on calo hit number
      int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k); // get OM number

      if (om_num == chosen_om){ // select only events for chosen OM
        tracking_event.push_back(event);
        std::vector<double> temp_ADC_vec(no_cells);
        
        std::vector<short> wave_k = wave->at(k); // get the waveform corresponding to calo hit k
        int fcr_k = fcr->at(k); // get the fcr corresponding to k
        
        std::cout << "event no: " << event << "\nOM no: " << om_num << "\nFCR: " << fcr_k << "\n------------" << std::endl;
        
        for (int bin = 0; bin < no_cells; bin++){ // bin: loop on waveform bin
        
          double wave_at_bin = wave_k.at(bin); // wave value corresponding to calo hit number and waveform bin
          int cell = get_cell(bin, fcr_k); // convert waveform bin to memory cell
            
          waveform->Fill(bin, wave_at_bin); // fill sum waveform without fcr shift
          
          if (event == single_waveform_event_no){
              swaveform->Fill(cell, wave_at_bin);
              single_waveT[cell] = wave_at_bin;
          }
          
          if (bin >= 960){
            temp_ADC_vec.at(cell) = 0;
          
            if (event == single_waveform_event_no){
              //single_waveT[cell] = 0;
            }
            
          } else{
            temp_ADC_vec.at(cell) = wave_at_bin;
            
            all_cells[cell].push_back(wave_at_bin); // input the wave value into xth vector of all_events (xth memory cell)
            
            if (event == single_waveform_event_no){
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
  make_hist(swaveform, "singular waveform", "Cell", "ADC", "waveform_histograms/single_waveform.png");
  make_hist(waveform, "full waveform", "Bin", "ADC", "waveform_histograms/waveform.png");
  
  // ------------------------------------------------------------------------------------------------------------------------
  // find average wavefom and calculate calibration constant for each memory cell

  std::vector<double> averages; // vector to store the average value over all events for each memory cell
  
  // define variables to plot against memory cell
  std::vector<double> means;
  std::vector<double> sigmas;
  std::vector<double> red_chi2s;
  std::vector<double> cell_nums;
  std::vector<double> no_events;
  
  Double_t empty_x_errT[1024];
  Double_t meansT[1024];
  Double_t sigmasT[1024];
  Double_t red_chi2sT[1024];
  Double_t cell_numsT[1024];
  Double_t no_eventsT[1024];
  
    for (int cell = 0; cell < no_cells; cell++){ // cell: loop on memory cell
    std::vector<double> ADC_vals = all_cells.at(cell); // select the memory cell from all_cells
    int no_vals = ADC_vals.size(); // number of ADC readings/entries
    
    double average = std::accumulate(ADC_vals.begin(), ADC_vals.end(), 0.0) / no_vals; // calculate average ADC across all events for this cell
    averages.push_back(average); // add average to averages vector
    avwaveform -> Fill(cell, average);

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
    TCanvas* canv_mc = new TCanvas(title, title, 900, 600);
    //gStyle->SetOptStat("eMR");
    gStyle->SetOptFit();
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus");
    hist->Draw("E1");
    canv_mc->SaveAs(save_name);

    // get fit parameters
    TF1 *fit = hist->GetFunction("gaus");
    double mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double chi2 = fit->GetChisquare();
    double ndof = fit->GetNDF();
    
    // fill variables to plot against memory cell
    cell_nums.push_back(cell+1);
    means.push_back(mean);
    sigmas.push_back(sigma);
    red_chi2s.push_back(chi2/ndof);
    no_events.push_back(no_vals);
    
    empty_x_errT[cell] = 0;
    meansT[cell] = mean;
    sigmasT[cell] = sigma;
    red_chi2sT[cell] = chi2 / ndof;
    cell_numsT[cell] = cell + 1;
    no_eventsT[cell] = no_vals;
    /*
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "MEAN: " << mean;
    std::cout << "-----------------------------------------------------------------\n";
    */
  }
  
  // plot graphs against memory cell
  make_tgraph(no_cells, cell_numsT, meansT, 1024, "Mean ADC readings", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC_tgraph.png");
  make_tgraph(no_cells, cell_numsT, sigmasT, 1024, "Standard deviation", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas_tgraph.png");
  make_tgraph(no_cells, cell_numsT, red_chi2sT, 1024, "Chi squares", "Memory cell", "Reduced chi square", "mem_cell_plots/chi2_tgraph.png");
  make_tgraph(no_cells, cell_numsT, no_eventsT, 1024, "Number of events", "Memory cell", "Number of events", "mem_cell_plots/no_events_tgraph.png");
  
  // make average waveform histogram
  make_hist(avwaveform, "average waveform", "Cell", "ADC", "waveform_histograms/average_waveform.png");
  
  // plot tgraphs for calibrated and uncalibrated single wave
  make_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT, "ADC by cell - single waveform", "Cell", "ADC", "pol1", "calibration_plots/singular_waveform_pre_cal.png");
  make_cal_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal.png", meansT);
  
  Double_t pre_cal_std_dev [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev, meansT, "no");

  Double_t post_cal_std_dev [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev, meansT, "yes");
  
  plot_standard_dev(pre_cal_std_dev, post_cal_std_dev, sizeof(post_cal_std_dev)/sizeof(post_cal_std_dev[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation.png", "calibration_plots/post_cal_standard_deviation.png", 1024);
  //for_mean_tgraph(all_cells, means_arr, chosen_om, tracking_event);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // calculate calibration constant for every 64 memory cells
  
  std::vector<std::vector<double> > all_cells64(64);

  Double_t meansT_64[64];
  Double_t sigmasT_64[64];
  Double_t red_chi2sT_64[64];
  Double_t cell_numsT_64[64];
  Double_t no_eventsT_64[64];
  
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
    TH1D *hist = new TH1D(title, title, no_vals, min_val, max_val);
    for (int val = 0; val < no_vals; val++){
      hist -> Fill(ADC_vals[val]);
    }

    // save in a canvas
    TCanvas* canv_mc = new TCanvas(title, title, 900, 600);
    gStyle->SetOptFit();
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus");
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
    
    /*
    // fill variables to plot against memory cell
    cell_nums64.push_back(cell+1);
    means64.push_back(mean);
    sigmas64.push_back(sigma);
    red_chi2s64.push_back(chi2/ndof);
    no_events64.push_back(no_vals);
    
    std::cout << "cell num: " << cell+1 << "\n";
    std::cout << "mean: " << mean << "\n";
    std::cout << "sigma: " << sigma << "\n";
    std::cout << "chi2: " << chi2/ndof << "\n";
    std::cout << "no events: " << no_vals << "\n";
    */
  }
  
  // plot graphs against memory cell
  make_tgraph(64, cell_numsT_64, meansT_64, 64, "Mean ADC readings overlayed", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC64_tgraph.png");
  make_tgraph(64, cell_numsT_64, sigmasT_64, 64, "Standard deviation overlayed", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas64_tgraph.png");
  make_tgraph(64, cell_numsT_64, red_chi2sT_64, 64, "Chi squares overlayed", "Memory cell", "Reduced chi square", "mem_cell_plots/chi264_tgraph.png");
  make_tgraph(64, cell_numsT_64, no_eventsT_64, 64, "Number of events overlayed", "Memory cell", "Number of events", "mem_cell_plots/no_events64_tgraph.png");
  
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
  make_cal_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT_64_full, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal64.png", meansT_64_full);
  
  Double_t pre_cal_std_dev64 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev64, meansT_64_full, "no");

  Double_t post_cal_std_dev64 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev64, meansT_64_full, "yes");

  plot_standard_dev(pre_cal_std_dev64, post_cal_std_dev64, sizeof(post_cal_std_dev64)/sizeof(post_cal_std_dev64[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation64.png", "calibration_plots/post_cal_standard_deviation64.png", 64);
  //for_mean_tgraph(all_cells, means_arr64_full, chosen_om, tracking_event);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // calculate calibration constant for every 16 memory cells
  
  std::vector<std::vector<double> > all_cells16(16);

  Double_t meansT_16[16];
  Double_t sigmasT_16[16];
  Double_t red_chi2sT_16[16];
  Double_t cell_numsT_16[16];
  Double_t no_eventsT_16[16];
  
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
    gStyle->SetOptFit();
    hist->SetXTitle("ADC");
    hist->SetYTitle("Frequency");
    hist->Fit("gaus");
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
    
    /*
    // fill variables to plot against memory cell
    cell_nums16.push_back(cell+1);
    means16.push_back(mean);
    sigmas16.push_back(sigma);
    red_chi2s16.push_back(chi2/ndof);
    no_events16.push_back(no_vals);
    
    std::cout << "cell num: " << cell+1 << "\n";
    std::cout << "mean: " << mean << "\n";
    std::cout << "sigma: " << sigma << "\n";
    std::cout << "chi2: " << chi2/ndof << "\n";
    std::cout << "no events: " << no_vals << "\n";
    */
  }

  // plot graphs against memory cell
  make_tgraph(16, cell_numsT_16, meansT_16, 16, "Mean ADC readings overlayed", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC16_tgraph.png");
  make_tgraph(16, cell_numsT_16, sigmasT_16, 16, "Standard deviation overlayed", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas16_tgraph.png");
  make_tgraph(16, cell_numsT_16, red_chi2sT_16, 16, "Chi squares overlayed", "Memory cell", "Reduced chi square", "mem_cell_plots/chi216_tgraph.png");
  make_tgraph(16, cell_numsT_16, no_eventsT_16, 16, "Number of events overlayed", "Memory cell", "Number of events", "mem_cell_plots/no_events16_tgraph.png");
  
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
  make_cal_tgraph_single_wave(1024, cell_numsT, single_waveT, empty_x_errT, sigmasT_16_full, "Calibrated ADC by cell - single waveform", "Cell", "ADC - Calibration constant", "pol1", "calibration_plots/singular_waveform_post_cal16.png", meansT_16_full);
  
  Double_t pre_cal_std_dev16 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, pre_cal_std_dev16, meansT_16_full, "no");

  Double_t post_cal_std_dev16 [cell_ordered_waveforms.size()];
  calc_standard_dev(cell_ordered_waveforms, post_cal_std_dev16, meansT_16_full, "yes");

  plot_standard_dev(pre_cal_std_dev16, post_cal_std_dev16, sizeof(post_cal_std_dev16)/sizeof(post_cal_std_dev16[0]), chosen_om, "calibration_plots/pre_cal_standard_deviation16.png", "calibration_plots/post_cal_standard_deviation16.png", 16);
  //for_mean_tgraph(all_cells, means_arr16_full, chosen_om, tracking_event);
  
  // ------------------------------------------------------------------------------------------------------------------------
  // overlayed plots
  
  TCanvas* canv_overlayed_64 = new TCanvas("Every 64 means overlyed", "Every 64 means overlyed", 900, 600);
  
  double* means_64 = slice(means, 1, 64);
  double* cell_nums_64 = slice(cell_nums, 1, 64);
  
  TGraph *gr64 = new TGraph(64, cell_nums_64, means_64);
  gr64->Draw("AL");
  
  for (int i = 1; i < 16; i++){
    double* means_slice = slice(means, 1+64*i, 64*(i+1));
    TGraph *gr = new TGraph(64, cell_nums_64, means_slice);
    gr->SetLineColor(i+1);
    gr->Draw("LSame");
  }
  
  canv_overlayed_64->SaveAs("mem_cell_plots/overlayed_means_every_64.png");
  
  TCanvas* canv_overlayed_16 = new TCanvas("Every 16 means overlyed", "Every 16 means overlyed", 900, 600);
  
  double* means_16 = slice(means, 1, 16);
  double* cell_nums_16 = slice(cell_nums, 1, 16);
  
  TGraph *gr16 = new TGraph(16, cell_nums_16, means_16);
  gr16->Draw("AL");
  
  for (int i = 1; i < 64; i++){
    double* means_slice = slice(means, 1+16*i, 16*(i+1));
    TGraph *gr = new TGraph(16, cell_nums_16, means_slice);
    gr->SetLineColor(i+1);
    gr->Draw("LSame");
  }
  
  canv_overlayed_16->SaveAs("mem_cell_plots/overlayed_means_every_16.png");
  
  // ------------------------------------------------------------------------------------------------------------------------
  /*
  std::cout << "CELL_ORDERED SIZE: " << cell_ordered_waveforms.size() << "\n";
  for (int i = 0; i < cell_ordered_waveforms.size(); i++){
    std::cout << "CELL_ORDERED[I] SIZE: " << cell_ordered_waveforms.at(i).size() << "\n";
  }
  for (int i = 0; i < cell_ordered_waveforms.size(); i++){
    for (int j = 0; j < cell_ordered_waveforms.at(i).size(); j++){
      std::cout << "---------------------------------------------------------------------------------------\n";
      std::cout << cell_ordered_waveforms.at(i).at(j) << "\n";
      std::cout << "---------------------------------------------------------------------------------------\n";
    }
  }
  */
  std::cout << "Number of events: " << counter << std::endl;
  TFile *newfile = new TFile("output_from_macro.root", "RECREATE");
  newfile->cd();
  waveform->Write();
  swaveform->Write();
  newfile->Close();
  
  std::cout << "End of macro.C" << std::endl;
  return 0;
}
