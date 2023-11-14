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

void make_tgraph(Int_t n, Double_t* x, Double_t* y, const char* title, const char* xtitle, const char* ytitle, const char* save_title){
  /*
  function to make and save a canvas with TGraph plot
  */
  TGraph *graph = new TGraph(n, x, y);
  
  TCanvas* canv = new TCanvas(title, title, 900, 600);
  graph->SetLineColor(9);
  graph->Draw();
  graph->GetXaxis()->SetLimits(0, 1024);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(xtitle);
  graph->GetYaxis()->SetTitle(ytitle);
  canv->SaveAs(save_title);
}

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
  
  std::vector<std::vector<double> > all_cells(no_cells); // memory cell * event vector, so we have the value for each event at each memory cell
  int max_entries = tree->GetEntries();
  
  TH1D *waveform = new TH1D("waveform", "waveform", no_cells, 0, no_cells);
  TH1D *swaveform = new TH1D("swaveform", "singular waveform", no_cells, 0, no_cells);
  TH1D *avwaveform = new TH1D("avwaveform", "average waveform", no_cells, 0, no_cells);

  for (int event = 0; event < max_entries; event++) {  // i: loop on event number
    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; k++) {  // k: loop on calo hit number
      int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k); // get OM number

      if (om_num == chosen_om){ // select only events for chosen OM
        
        std::vector<short> wave_k = wave->at(k); // get the waveform corresponding to calo hit k
        int fcr_k = fcr->at(k); // get the fcr corresponding to k
        
        std::cout << "event no: " << event << "\nOM no: " << om_num << "\nFCR: " << fcr_k << "\n------------" << std::endl;
        
        for (int bin = 0; bin < no_cells; bin++){ // bin: loop on waveform bin
        
          double wave_at_bin = wave_k.at(bin); // wave value corresponding to calo hit number and waveform bin
          int cell = get_cell(bin, fcr_k); // convert waveform bin to memory cell
            
          waveform->Fill(bin, wave_at_bin); // fill sum waveform without fcr shift
          
          if (event == single_waveform_event_no){
              swaveform->Fill(cell, wave_at_bin);
          }
          
          if (bin >= 960){
            continue;
          } else{
            all_cells[cell].push_back(wave_at_bin); // input the wave value into xth vector of all_events (xth memory cell)
          }
        }
      }
    counter++;
    }
  }

  std::vector<double> averages; // vector to store the average value over all events for each memory cell
  
  // define variables to plot against memory cell
  Double_t means[no_cells];
  std::vector<double> means_vec;
  Double_t mean_errs[no_cells];
  Double_t sigmas[no_cells];
  Double_t sigma_errs[no_cells];
  Double_t red_chi2s[no_cells];
  Double_t cell_nums[no_cells];
  Double_t no_events[no_cells];
  
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
    double mean_err = fit->GetParError(1);
    double sigma = fit->GetParameter(2);
    double sigma_err = fit->GetParError(2);
    double chi2 = fit->GetChisquare();
    double ndof = fit->GetNDF();
    
    // fill variables to plot against memory cell
    std::cout << "chi square: " << chi2 << "\n";
    cell_nums[cell] = cell + 1;
    means[cell] = mean;
    means_vec.push_back(mean);
    mean_errs[cell] = mean_err;
    sigmas[cell] = sigma;
    sigma_errs[cell] = sigma_err;
    red_chi2s[cell] = chi2 / ndof;
    no_events[cell] = no_vals;
  }
  
  Double_t means_1to64(means_vec.begin(), means_vec.begin()+63);
  Double_t means_64to128(means_vec.begin()+64, means_vec.begin()+127);
  
  std::cout << "Number of events: " << counter << std::endl;
  TFile *newfile = new TFile("output_from_macro.root", "RECREATE");
  newfile->cd();
  waveform->Write();
  swaveform->Write();
  newfile->Close();
  
  // plot graphs against memory cell
  make_tgraph(no_cells, cell_nums, means, "Mean ADC readings", "Memory cell", "Mean ADC", "mem_cell_plots/meanADC_tgraph.png");
  make_tgraph(no_cells, cell_nums, sigmas, "Standard deviation", "Memory cell", "Standard deviation", "mem_cell_plots/sigmas_tgraph.png");
  make_tgraph(no_cells, cell_nums, red_chi2s, "Chi squares", "Memory cell", "Reduced chi square", "mem_cell_plots/chi2_tgraph.png");
  make_tgraph(no_cells, cell_nums, no_events, "Number of events", "Memory cell", "Number of events", "mem_cell_plots/no_events_tgraph.png");
  
  make_hist(swaveform, "singular waveform", "Cell", "ADC", "waveform_histograms/single_waveform.png");
  make_hist(waveform, "full waveform", "Bin", "ADC", "waveform_histograms/waveform.png");
  make_hist(avwaveform, "average waveform", "Cell", "ADC", "waveform_histograms/average_waveform.png");
  
  std::cout << "End of macro.C" << std::endl;
  return 0;
}
