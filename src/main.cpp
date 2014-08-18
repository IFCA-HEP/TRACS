#include <SMSDetector.h>

int main()
{


  double pitch = 80.;
  double width = 80.;
  double depth = 200.;
  int nns = 0;

  int n_cells_x = 150;
  int n_cells_y = 150;

  char bulk_type = 'p';
  char implant_type = 'n';

  SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y);

  detector.solve_w_u();
  Function * w_u = detector.get_w_u();

  // Plot solution
  plot(*w_u,"Weighting Potential","auto");

  double v_bias = 150.0;
  double v_depletion = 61.0;
  detector.set_voltages(v_bias, v_depletion);

  detector.solve_d_u();
  Function * d_u = detector.get_d_u();
  // Plot solution
  plot(*d_u,"Drift Potential","auto");

  detector.solve_w_f_grad();
  Function * w_f_grad = detector.get_w_f_grad();
  // Plot solution
  plot((*w_f_grad)[1],"Weighting Field (Y)","auto");

  detector.solve_d_f_grad();
  Function * d_f_grad = detector.get_d_f_grad();
  // Plot solution
  plot((*d_f_grad)[1],"Drifting Field (Y)","auto");



  interactive();


  return 0;
}
