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

#include "HistUtils.h"
using namespace std;

map<TString,TString> plot_names;
vector<TString> string_labels;
bool debug;

//----------------------------------------------------------------------------------------
void draw(TString human, TString var, TString region, TString imatch, TString var_type, TString log, map<TString,int> string_map) {
  TCanvas c1("c1","",700,600);
  TString naivety("M");
  TString basedir("data/human-beings/" + human + "/" + naivety);
  TH1F hdata = make_hist(basedir + "/data/" + var + "/" + region + "-" + imatch + ".txt", "data", var_type, log, string_map);
  TH1F hsimu = make_hist(basedir + "/simu/" + var + "/" + region + "-" + imatch + ".txt", "simu", var_type, log, string_map);

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

  hsimu.SetLineColor(kRed+1);
  hsimu.SetMarkerSize(0);
  hsimu.SetLineWidth(4);
  hsimu.Draw("ehist same");
  hdata.Draw("e same");
  c1.SetLogx(log.Contains("x"));
  c1.SetLogy(log.Contains("y"));
  TString plotdir("/var/www/sharing/dralph/work/plotting/human-beings/" + human + "/M/" + var + "/plots");
  gSystem->mkdir(plotdir, true);
  // c1.SaveAs(plotdir + "/" + region + "-" + imatch + ".png");
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
  map<TString,int> string_map;
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
	draw(humans[ih], "found_strings"      , regions[ir], imatches[im], "string", "", string_map);
	draw(humans[ih], "ali_from"           , regions[ir], imatches[im], "double", "", string_map);
	draw(humans[ih], "ali_to"             , regions[ir], imatches[im], "double", "", string_map);
	draw(humans[ih], "evalue"             , regions[ir], imatches[im], "double", "x",string_map);
	draw(humans[ih], "vd_insertion_length", regions[ir], imatches[im], "double", "", string_map);
	draw(humans[ih], "dj_insertion_length", regions[ir], imatches[im], "double", "", string_map);
      }
    }
  }
}
