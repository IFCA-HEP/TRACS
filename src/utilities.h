#ifndef UTILITIES_H
#define UTILITIES_H


#include <dolfin.h>

#include <TH2D.h>
#include <TString.h>

using namespace dolfin;

namespace utilities
{
  TH2D export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max)
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
}

#endif // UTILITIES_H
