#ifndef UTILITIES_H
#define UTILITIES_H


#include <dolfin.h>

#include <TH2D.h>
#include <TString.h>

#include "qcustomplot.h"

using namespace dolfin;

namespace utilities
{
  TH2D export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max);
  void paint_TH2D_qcp(TH2D hist, QCPColorMap * color_map);
}

#endif // UTILITIES_H
