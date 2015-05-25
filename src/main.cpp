#include "SMSDetector.h"
#include "utilities.h"
#include "Carrier.h"
#include "CarrierCollection.h"

#include <TFile.h>
#include <TH2D.h>

#include <fstream>
#include <iterator>

/*
 ************************TRACS-MAIN********************************
 *
 * Main function for TRACS program it performs all the required
 * tasks 
 *
 * Required data: - Pitch
 * 		  - Width
 * 		  - Depth
 * 		  - nns
 * 		  - Starting point (X,Y)
 * 		  - End point (X,Y)
 * 		  - cell x
 * 		  - cell y
 * 		  - Bulk type
 * 		  - Implant Type
 * 		  - Vias Voltaje
 * 		  - Depletion Voltaje
 * 		  - Time Step
 * 		  - Total time
 * 		  - Neighbour strips
 * 		  - 
 *
 *                            TODO
 *
 *
 *                     - Wrapping for loop depth -> DONE
 *                     - Wrapping for loop voltaje -> DONE
 *                     - Point to a different file to read each iteration -> WIP
 *                     - Read from file -> TODO
 *                     - Try/debug -> TODO
 *                     - Call it a final producto -> Long till done
 *
 *
 *
 *
 */

int main()
{

  double pitch = 100.;
  double width = 20.;
  double depth = 300.;
  int nns = 3; // Neighbour Strips

  double x_min = 0.;
  double x_max = pitch * (2*nns+1);
  double y_min = 0.;
  double y_max = depth;

  int n_cells_x = 150;
  int n_cells_y = 150;

  char bulk_type = 'p';
  char implant_type = 'n';

  // for loop that scans the detector in all desired depths
// 	for (int iteration_depth = 0; iteration_depth <= depth; iteration_depth += 5)
//	{
  	
  	// for loop that simulates currents for voltajes [0-500]V 
  	// in steps of 50V. The position and angle at which the laser
  	// illuminates the detector is controlled via the read file
  	// which contains the CC position at generation time.
//		for (int voltaje_iteration = 0; voltaje_iteration < 11; ++voltaje_iteration)
//		{

		      parameters["allow_extrapolation"] = true;



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
		      // File file_u("periodic.pvd");
		      // file_u << *w_u;

		      double v_bias = 500.;//*voltaje_iteration; // Vias Voltage
		      double v_depletion = 180.0; // Depletion Voltaje
		      detector.set_voltages(v_bias, v_depletion);

		      detector.solve_d_u();
		      Function * d_u = detector.get_d_u();
		      // Plot solution
		      // plot(*d_u,"Drift Potential","auto");

		      detector.solve_w_f_grad();
		      Function * w_f_grad = detector.get_w_f_grad();
		      // Plot solution
		      // plot((*w_f_grad)[1],"Weighting Field (Y)","auto");

		      detector.solve_d_f_grad();
		      Function * d_f_grad = detector.get_d_f_grad();
		      // Plot solution
		      // plot((*d_f_grad)[1],"Drifting Field (Y)","auto");

		      // Create carrier and observe movement
		      SMSDetector * dec_pointer = &detector;

		      double dt = 6.25e-12;
		      double max_time = 10.e-9;

		      // get number of steps from time
		      const int max_steps = (int) std::floor(max_time / dt);

		      // filename to get carrier distribution
		      
		      QString filename = "etct.carriers";//+QString::number(iteration_depth)+".carriers";
		      CarrierCollection * carrier_collection = new CarrierCollection(dec_pointer);
		      carrier_collection->add_carriers_from_file(filename);

		      // init arrays
		      std::valarray<double> curr_elec((size_t) max_steps);
		      std::valarray<double> curr_hole((size_t) max_steps);
		      std::valarray<double> curr_total((size_t) max_steps);

		      // array of shifts
		      const int n_steps_y = 60; // Simulation Steps
		      double border = 7.15284952e+01; // WUT??
		      double c_value = 150.; // WUT??
		      double step_size_y = 2*3.51; // UNITS??
		      std::vector<double>  shift_y_array(n_steps_y+1);

		      for (int i = 0; i < n_steps_y + 1; i++ ) 
		      {
			shift_y_array[i] = i*step_size_y - c_value - border;
		      }

		      TString hist_name = "curr_map_total";
		      TString hist_title = "curr_map_total";
		      TH2D curr_map_total = TH2D(hist_name, hist_title,
						 n_steps_y + 1, y_min, y_max + 2*border,
						 max_steps, 0.0, max_time );

		      TH2D curr_map_shaped = TH2D("curr_map_shaped", "curr_map_shaped",
						n_steps_y + 1, y_min, y_max + 2*border,
						max_steps, 0.0, max_time - dt);

		      // For loop that does all the steps
		      for (int i = 0; i < n_steps_y + 1; i++) 
		      {
			  std::cout << "### Step number ### " << i << " of " << n_steps_y + 1 << std::endl;
			  curr_hole = 0;
			  curr_elec = 0;
			  curr_total = 0;

			  // Simulate the drift of both electrons and holes
			  carrier_collection -> simulate_drift(dt, max_time, 0.0, shift_y_array [i],  curr_elec, curr_hole);

			  // Add up total charge
			  curr_total = curr_elec + curr_hole;

			  QVector<double> x_elec(max_steps), y_elec(max_steps);
			  QVector<double> x_hole(max_steps), y_hole(max_steps);
			  QVector<double> x_total(max_steps), y_total(max_steps);

			  // Compute time + format vectors for writting to file
			  for (int j=0; j < max_steps; j++)
			  {
			     x_elec[j] = j*dt;
			     x_hole[j] = j*dt;
			     x_total[j] = j*dt;
			     y_elec[j] = curr_elec[j];
			     y_hole[j] = curr_hole[j];
			     y_total[j] = curr_total[j];
			     curr_map_total.SetBinContent(i+1,j+1, y_total[j] );
			  }

			  // RC shaping
			  double C = 15e-12; // 15 pF
			  double RC = 50.*C; // 50 ohms*C
			  double dt = 6.25e-12;
			  double alfa = dt/(RC+dt);
			  QVector<double> y_shaped(max_steps);
			  y_shaped[0]=y_total[0];
			  curr_map_shaped.SetBinContent(i+1, 1, y_shaped[0]);
			  for (int j = 1; j <max_steps; j++) 
			  {
			    y_shaped[j]=alfa*y_total[j]+(1-alfa)*y_shaped[j-1];
			    curr_map_shaped.SetBinContent(i+1,j+1,y_shaped[j]);
			  }


			  // save results
			  // REDUNDANT STEP??
			  //
			  // TODO: Check
			  //
			  QVector<QVector<double>> raw_results;
			  raw_results.resize(4);
			  raw_results[0] = x_total;
			  raw_results[1] = y_total;
			  raw_results[2] = y_elec;
			  raw_results[3] = y_hole;

			  QString out_filename = "result_step_"+QString::number(i)+".plot";

			  utilities::write_results_to_file(out_filename, raw_results);

		      }

		      // Open a ROOT file to save result
//		      TFile *tfile = new TFile("results.root","RECREATE" );
//		      curr_map_total.Write();
//		      curr_map_shaped.Write();
//		      tfile->Close();

//		} // End for loop over voltajes
    
//	} // End for loop over depths

  // No plot now
  // interactive(); 
  // Needed so that plots are not closed as program finishes executing


  return 0;
}
