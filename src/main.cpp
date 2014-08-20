#include <SMSDetector.h>

#include <utilities.h>

#include <Carrier.h>

#include <TFile.h>

#include <fstream>
#include <iterator>

int main()
{

  parameters["allow_extrapolation"] = true;

  double pitch = 80.;
  double width = 30.;
  double depth = 200.;
  int nns = 3;

  double x_min = 0.;
  double x_max = pitch * (2*nns+1);
  double y_min = 0.;
  double y_max = depth;

  int n_cells_x = 150;
  int n_cells_y = 150;

  char bulk_type = 'p';
  char implant_type = 'n';

  SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y);

  detector.solve_w_u();
  Function * w_u = detector.get_w_u();

  // Plot solution
  //plot(*w_u,"Weighting Potential","auto");

  TH2D w_u_hist = utilities::export_to_histogram(*w_u, "hist", "hist", n_cells_x ,x_min,x_max, n_cells_y,y_min,y_max);

  TFile file("hist.root", "recreate");
  w_u_hist.Write();
  file.Close();

  // Save solution in VTK format
  File file_u("periodic.pvd");
  file_u << *w_u;

  double v_bias = 150.0;
  double v_depletion = 61.0;
  detector.set_voltages(v_bias, v_depletion);

  detector.solve_d_u();
  Function * d_u = detector.get_d_u();
  // Plot solution
  //plot(*d_u,"Drift Potential","auto");

  detector.solve_w_f_grad();
  Function * w_f_grad = detector.get_w_f_grad();
  // Plot solution
  //plot((*w_f_grad)[1],"Weighting Field (Y)","auto");

  detector.solve_d_f_grad();
  Function * d_f_grad = detector.get_d_f_grad();
  // Plot solution
  //plot((*d_f_grad)[1],"Drifting Field (Y)","auto");

  // Create carrier and observe movement
  SMSDetector * dec_pointer = &detector;



  std::vector<double> result(10000);
  std::vector<double> temp(10000);

  int times = 100000;

  for (int i = 0 ; i < times; i++)
  {
    Carrier electron = Carrier( 'h', 1. , 300., 0.3 , dec_pointer, 0.0);
    temp = electron.simulate_drift( 1e-12 , 10e-9);
    std::cout << "Number of electron: " << i << std::endl;
  }

  std::ofstream output_file("./example.txt");
  output_file << "Results vector[]={";
  std::copy(result.begin(), result.end(), std::ostream_iterator<double>(output_file , ", "));
  output_file << "}" <<  std::endl;

  //interactive();


  return 0;
}
