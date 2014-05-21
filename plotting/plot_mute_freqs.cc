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

//----------------------------------------------------------------------------------------
void plot_mute_freqs() {
  gStyle->SetOptStat(0);
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
