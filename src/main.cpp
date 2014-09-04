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

  //TH2D w_u_hist = utilities::export_to_histogram(*w_u, "hist", "hist", n_cells_x ,x_min,x_max, n_cells_y,y_min,y_max);

  //TFile file("hist.root", "recreate");
  //w_u_hist.Write();
  //file.Close();

  // Save solution in VTK format
  File file_u("periodic.pvd");
  file_u << *w_u;

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

  // Create carrier and observe movement
  SMSDetector * dec_pointer = &detector;



  double dt = 1.e-11;
  double max_time = 10.e-9;
  // get number of steps from time
  const int max_steps = (int) std::floor(max_time / dt);


  std::valarray<double> curr_holes((size_t) max_steps);
  std::valarray<double> curr_elec((size_t) max_steps);
  int times = 100;

  curr_holes = 0;
  curr_elec = 0;

  Carrier hole('h', 1. , 300., 190. , dec_pointer, 0.0);
  Carrier electron('e', 1. , 300., 190. , dec_pointer, 0.0);

  for (int i = 0 ; i < times; i++)
  {
    curr_holes += hole.simulate_drift( dt , max_time, 300., 190.);
    curr_elec += electron.simulate_drift( dt , max_time, 300., 190.);
    std::cout << "Number of carrier: " << 2*i << std::endl;
  }

  std::ofstream output_file("./example.txt");
  output_file << "curr_holes[]={";
  for ( const double &value : curr_holes )
  {
      output_file << value << ", ";
  }
  output_file << "}" <<  std::endl;

  output_file << "curr_elec[]={";
  for ( const double &value : curr_elec )
  {
    output_file << value << ", ";
  }
  output_file << "}" <<  std::endl;

  output_file.close();

  interactive();


  return 0;
}
