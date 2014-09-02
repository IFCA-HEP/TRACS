#include "utilities.h"


// function to export to a root histogram from a dolfin 1d function
TH2D utilities::export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max)
{
  TH2D hist = TH2D(hist_name,hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;
  for (int i=1; i<=n_bins_x; i++) {
    for (int j=1; j<=n_bins_y; j++) {
    double x_value = (i - 0.5)*step_x;
    double y_value = (j - 0.5)*step_y;
    hist.Fill(x_max-x_value,y_max-y_value, func(x_value, y_value));
    }
  }
  return hist;
}

// function to pain a TH2D root histogram in a color map
void  utilities::paint_TH2D_qcp(TH2D hist, QCPColorMap * color_map)
{
  // get some required variables
  int n_bins_x = hist.GetNbinsX();
  int n_bins_y = hist.GetNbinsY();
  double x_min = hist.GetXaxis()->GetXmin();
  double x_max = hist.GetXaxis()->GetXmax();
  double y_min = hist.GetYaxis()->GetXmin();
  double y_max = hist.GetYaxis()->GetXmax();
  color_map->data()->setSize(n_bins_x,n_bins_y);
  color_map->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double content = hist.GetBinContent(x,y);
      color_map->data()->setCell(x, y, content);
    }
  }
  color_map->setGradient(QCPColorGradient::gpPolar);
  color_map->rescaleDataRange(true);
}
