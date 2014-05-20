#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
using namespace std;

map<TString,TString> plot_names;
vector<TString> string_labels;
map<TString,int> string_map;
bool debug;

//----------------------------------------------------------------------------------------
void set_bins(vector<double> values, int n_bins, bool is_log_x, TString var_type, double *xbins) {
  if (is_log_x) {
    double log_xmin = log(powf(10.0f, floorf(log10f(values[0]))));  // round down to the nearest power of 10
    double log_xmax = log(powf(10.0f, ceilf(log10f(values[values.size()-1])))); // round up to the nearest power of 10
    double log_dx = (log_xmax - log_xmin) / n_bins;
    log_xmin -= 0.1*log_dx;  // expand a little to avoid overflows
    log_xmax += 0.1*log_dx;
    log_dx = (log_xmax - log_xmin) / n_bins;
    for (int ib=0; ib<=n_bins; ib++) {
      double low_edge = exp(log_xmin + ib * log_dx);
      xbins[ib] = low_edge;
    }
  } else {
    double dx = (values[values.size()-1] - values[0]) / n_bins;
    double xmin = values[0] - 0.1*dx;  // expand a little to avoid overflows
    double xmax = values[values.size()-1] + 0.1*dx;  // expand a little to avoid overflows
    if (var_type=="string") {
      xmin = 0.5;
      xmax = string_map.size() + 0.5;
    }
    dx = (xmax - xmin) / n_bins;  // then recalculate dx
    for (int ib=0; ib<=n_bins; ib++) {
      xbins[ib] = xmin + ib * dx;
    }
  }
}
//----------------------------------------------------------------------------------------
TH1F make_hist(TString infname, TString data_type, TString var_type, TString log) {
  int n_bins = 30;
  if (var_type=="string")
    n_bins = string_map.size();
  double xbins[n_bins+1];  // NOTE has to be n_bins + 1

  ifstream ifs(infname);
  if(!ifs.is_open()) {
    cout << "    " << infname << " d.n.e." << endl;
    return TH1F();
  }
  string line;
  vector<double> values;
  while(getline(ifs,line)) {
    stringstream ss(line);
    if (var_type=="double") {
      float value;
      ss >> value;
      values.push_back(value);
    } else if (var_type=="string") {
      string value;
      ss >> value;
      if (string_map.find(value) == string_map.end()) {
	cout << value << " not found!" << endl;
	assert(0);
      }
      values.push_back(double(string_map[value] + 1));  // string_map is a crappy name, but I can't come up with anything better
                                                        // NOTE add one to get to 1-base indexing for hist bins
   } else {
      assert(0);
    }
  }
  if (values.size() == 0)
    return TH1F();
  sort(values.begin(), values.end());
  cout << "  " << values.size() << " values in " << infname << endl;

  set_bins(values, n_bins, log.Contains("x"), var_type, xbins);
  TH1F hist("h"+data_type, "", n_bins, xbins);
  hist.Sumw2();
  for (unsigned iv=0; iv<values.size(); iv++) {
    hist.Fill(values[iv]);
  }
  // make sure there's no overflows
  if(hist.GetBinContent(0) != 0 || hist.GetBinContent(hist.GetNbinsX()+1) != 0) {
    cout << infname << endl;
    for (unsigned iv=0; iv<values.size(); iv++) {
      cout
	<< setw(12) << values[iv];
    }
    cout << endl;
    for (int ib=0; ib<hist.GetNbinsX()+2; ib++) {
      cout
      	<< setw(12) << ib
      	<< setw(12) << hist.GetBinLowEdge(ib)
      	<< setw(12) << hist.GetBinContent(ib)
      	<< endl;
    }
    assert(0);
  }

  hist.Scale(1./hist.Integral());
  
  return hist;
}
//----------------------------------------------------------------------------------------
void draw(TString human, TString var, TString region, TString imatch, TString var_type, TString log) {
  TCanvas c1("c1","",700,600);
  TString naivety("M");
  TString basedir("data/human-beings/" + human + "/" + naivety);
  TH1F hdata = make_hist(basedir + "/data/" + var + "/" + region + "-" + imatch + ".txt", "data", var_type, log);
  TH1F hsimu = make_hist(basedir + "/simu/" + var + "/" + region + "-" + imatch + ".txt", "simu", var_type, log);

  if (hdata.GetEntries() == 0 || hsimu.GetEntries() == 0) {
    if (debug)
      cout << "  no entries for " << var << " " << region << " " << imatch << endl;
    return;
  }

  double xmin(min(hdata.GetBinLowEdge(1), hsimu.GetBinLowEdge(1)));
  double xmax(max(hdata.GetXaxis()->GetBinUpEdge(hdata.GetNbinsX()), hsimu.GetXaxis()->GetBinUpEdge(hsimu.GetNbinsX())));
  TH1F hframe("hframe", "", hdata.GetNbinsX(), xmin, xmax);
  if (var_type=="string") {
    assert(hframe.GetNbinsX() == int(string_labels.size()));
    for (int ib=1; ib<hframe.GetNbinsX()+1; ib++) {
      hframe.GetXaxis()->SetBinLabel(ib, string_labels[ib-1]);
    }
  }
  hframe.SetMaximum(1.35*(max(hdata.GetMaximum(), hsimu.GetMaximum())));
  hframe.SetTitle(region + " match " + imatch + ";" + plot_names[var] + ";entries");
  hframe.Draw("txt");
  TLegend leg(0.65, 0.75, 0.94, 0.9);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(&hdata, "data", "lp");
  leg.AddEntry(&hsimu, "simulation", "l");
  leg.Draw();

  // for (int ib=0; ib<hdata.GetNbinsX()+2; ib++)
  //   cout
  //     << setw(12) << ib
  //     << setw(12) << hdata.GetBinLowEdge(ib)
  //     << setw(12) << hdata.GetBinContent(ib)
  //     << endl;

  // cout << "mins "
  //      << setw(12) << hdata.GetBinLowEdge(1)
  //      << setw(12) << hsimu.GetBinLowEdge(1)
  //      << setw(12) << xmin
  //      << endl;
  // cout << "maxs "
  //      << setw(12) << hdata.GetXaxis()->GetBinUpEdge(hdata.GetNbinsX())
  //      << setw(12) << hsimu.GetXaxis()->GetBinUpEdge(hsimu.GetNbinsX())
  //      << setw(12) << xmax
  //      << endl;

  hsimu.SetLineColor(kRed+1);
  hsimu.SetMarkerSize(0);
  hsimu.SetLineWidth(4);
  hsimu.Draw("ehist same");
  hdata.Draw("e same");
  c1.SetLogx(log.Contains("x"));
  c1.SetLogy(log.Contains("y"));
  TString plotdir("/var/www/sharing/dralph/work/plotting/human-beings/" + human + "/M/" + var + "/plots");
  gSystem->mkdir(plotdir, true);
  c1.SaveAs(plotdir + "/" + region + "-" + imatch + ".png");
}
//----------------------------------------------------------------------------------------
void plot_validation() {
  debug = false;
  gStyle->SetOptStat(0);

  // correspondence between variable name and variable title
  plot_names["found_strings"] = "find order";
  plot_names["ali_from"] = "alignment start";
  plot_names["ali_to"] = "alignment end";
  plot_names["evalue"] = "E value";
  plot_names["vd_insertion_length"] = "VD insertion";
  plot_names["dj_insertion_length"] = "DJ insertion";
  // correspondence between string value and integer
  string_labels.push_back("vjd");
  string_labels.push_back("vj");
  string_labels.push_back("vd");
  string_labels.push_back("jvd");
  string_labels.push_back("jd");
  string_labels.push_back("v");
  string_labels.push_back("None");
  string_labels.push_back("vdj");
  string_labels.push_back("jv");
  string_labels.push_back("jdv");
  string_labels.push_back("d");
  string_labels.push_back("j");
  string_labels.push_back("dvj");
  for (unsigned is=0; is<string_labels.size(); is++)
    string_map[string_labels[is]] = is;

  vector<TString> humans;
  humans.push_back("A");
  humans.push_back("B");
  humans.push_back("C");
  for (unsigned ih=0; ih<humans.size(); ih++) {
    vector<TString> regions,imatches;
    regions.push_back("v");
    regions.push_back("d");
    regions.push_back("j");
    imatches.push_back("0");
    imatches.push_back("1");
    imatches.push_back("2");
    for (unsigned ir=0; ir<regions.size(); ir++) {
      for (unsigned im=0; im<imatches.size(); im++) {
	draw(humans[ih], "found_strings"      , regions[ir], imatches[im], "string", "");
	draw(humans[ih], "ali_from"           , regions[ir], imatches[im], "double", "");
	draw(humans[ih], "ali_to"             , regions[ir], imatches[im], "double", "");
	draw(humans[ih], "evalue"             , regions[ir], imatches[im], "double", "x");
	draw(humans[ih], "vd_insertion_length", regions[ir], imatches[im], "double", "");
	draw(humans[ih], "dj_insertion_length", regions[ir], imatches[im], "double", "");
      }
    }
  }
}
