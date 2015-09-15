#include "SMSDetector.h"
#include "utilities.h"
#include "Carrier.h"
#include "CarrierCollection.h"
#include <functional>

#include <TFile.h>
#include <TH2D.h> // 2 Dimensional ROOT histogram
#include <TH1D.h> // 1 Dimesional ROOT histogram

#include <fstream> 
#include <iterator>
#include <thread>
#include <limits>  // std::numeric_limits
#include <functional>

extern TH1D *H1DConvolution( TH1D *htct , Double_t Cend=0. ) ; 
void call_from_thread(CarrierCollection & cCollection, double dt, double max_time, double shift_x, double shift_y, std::vector<double> curr_elec, std::vector<double> curr_hole, int thr_id);

void call_from_thread2(CarrierCollection &cCollection, double dt, double max_time, double shift_x, double shift_y, double curr_elec, int id);

//------------

void call_from_thread(CarrierCollection & cCollection, double dt, double max_time, double shift_x, double shift_y, std::vector<double> curr_elec, std::vector<double> curr_hole, int id)
{
	int nTimeSteps = curr_elec.size();
	std::valarray<double> i_elec((size_t) nTimeSteps);	
	std::valarray<double> i_hole((size_t) nTimeSteps);
	cCollection.simulate_drift( dt, max_time, shift_x, shift_y, i_elec, i_hole, id);
	curr_elec.assign(std::begin(i_elec), std::end(i_elec));
	curr_hole.assign(std::begin(i_hole), std::end(i_hole));
}

void call_from_thread2(CarrierCollection &cCollection, double dt, double max_time, double shift_x, double shift_y, double curr_elec, int thr_id)
{
	std::cout << "Hello World" << std::endl;
}

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
 * 		  - Bias Voltage
 * 		  - Depletion Voltage
 * 		  - Time Step
 * 		  - Total time
 * 		  - Neighbour strips
 * 		  - 
 */

int main()
{
	int nThreads = 0, nns = 0, n_cells_y = 0, n_cells_x = 0, waveLength = 0, n_vSteps = 0, n_zSteps;
	double pitch = 0, width = 0, depth = 0, temp = 0, trapping = 0, fluence = 0, C = 0, dt = 0, max_time = 0, v_bias = 0, v_init = 0, deltaV = 0, v_max = 0, v_depletion = 0, margin = 0, deltaZ = 0;
	std::string scanType = "defaultString";
	char bulk_type = '\0', implant_type = '\0';


	utilities::parse_config_file("Config.TRACS", depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y, bulk_type, implant_type, waveLength, scanType, C, dt, max_time, v_bias, v_init, deltaV, v_max, v_depletion, margin, deltaZ);
//	// Number of threads
//	int nThreads = 1;
	std::thread t[nThreads-1];
//
//	// Physical properties of the detector
//	double pitch = 80.00; // in microns
//	double width = 25.00; // in microns
//	double depth = 300.; // in microns
//	int nns = 2; // Neighbour Strips
//
//	// "Enviroment" properties
//	double temp = 300; // Temperature of the detector in Kelvin
//	double trapping = 3.e-9; // Trapping time in seconds
//	double fluence = 1.e14; // Fluence in neutron eq
	std::string trap;
	if (fluence <= 0) // if no fluence -> no trapping
	{
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
	}
	else
	{
		trap = std::to_string((int) std::floor(1.e9*trapping));
	}

//	// Mesh properties
//	int n_cells_x = 150;
//	int n_cells_y = 150;
//
//	// Doping
//	char bulk_type = 'p';
//	char implant_type = 'n';
//
//	// Laser
//	int waveLength = 1064; // in nm
//	std::string scanType = "edge"; // edge/top/bottom
//
//
//	// RC shaping
//	double C = 5.e-12; // in Farad
//  //double RC = 50.*C; // Ohms*Farad
//	double dt = 5.e-11; // in seconds this is also the simulation timestep
//	double max_time = 1.5e-8;
//
//	// Voltages (might be modified inside the for loops)
//	double v_bias = 350.;// Still needed, should get rid of it later
//	double v_init = 350;
//	double deltaV = 350;
//	double v_max = 350.; 
//	double v_depletion = 250.0; // Depletion Voltaje
//
//	// Shifting CC from laser beam 
//	double margin = 10;
//	double deltaZ = 3.0000; // microns 
	TH1D *hnoconv , *hconv ;


	n_zSteps  = (int) std::floor((depth+margin*2)/deltaZ); // Simulation Steps
	n_vSteps = (int) std::floor((v_init-v_max)/deltaV);


	std::cout << depth << std::endl;

	// Convert relevant simulation numbers to string for fileNaming	
	std::string dtime = std::to_string((int) std::floor(dt*1.e12));
	std::string neigh = std::to_string(nns);
	std::string stepV = std::to_string((int) std::floor(deltaV)); 
	std::string stepZ = std::to_string((int) std::floor(deltaZ));
	std::string cap = std::to_string((int) std::floor(C*1.e12));
	std::string z_step  = std::to_string((int) std::floor(deltaZ));
	std::string voltage = std::to_string((int) std::floor(v_bias));

	parameters["allow_extrapolation"] = true;

	SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence);

	detector.set_voltages(v_init, v_depletion);

	// Create carrier and observe movement
	SMSDetector * dec_pointer = &detector;

	// get number of steps from time
	int n_time_steps = (int) std::floor(max_time / dt);

	std::string file_carriers = "etct.carriers";
	QString filename = "etct.carriers";
	CarrierCollection * carrier_collection = new CarrierCollection(dec_pointer);

	// carrier_collection is now a #thr-dimensional vector
	carrier_collection->add_carriers_from_file(filename, nThreads); // input #threads

	// init arrays
	// These should are now #thr-dimensioal vectors

//	std::valarray<std::valarray<double>> i_elec;
//	std::valarray<std::valarray<double>> i_hole;
//	std::valarray<double> i_total((size_t) n_time_steps);

	
//	std::valarray<double> i_elec((size_t) n_time_steps);
//	std::valarray<double> i_hole((size_t) n_time_steps);
	std::valarray<double> i_total((size_t) n_time_steps);
	std::vector< std::vector<double> > vva_elec(nThreads, std::vector<double>(n_time_steps)); // vector of size N-threads where every element is a vector
	std::vector< std::vector<double> > vva_hole(nThreads, std::vector<double>(n_time_steps)); // same as above for holes current
	
	// Shifting  of Charge Carriers
	// the laser gets shifted in x and/or z direction depending on the arrays 
	// fed to the simulate_drift method inside the main for loop
	std::vector<double>  z_shifts(n_zSteps+1);
	double x_shift = 0.0; // laser shift in X axis to center laser focus over read-out strip
	std::vector<double>  voltages(n_vSteps+1);

	for (int i = 0; i < n_vSteps + 1; i++ ) 
	{
		voltages[i] = (i*deltaV)+v_init;
	}

	for (int i = 0; i < n_zSteps + 1; i++ ) 
	{
		z_shifts[i] = (i*deltaZ)-margin;
	}

	TString hist_name = "i_ramo";
	TString hist_title = "i_ramo";
	TH2D i_ramo = TH2D(hist_name, hist_title,
			n_time_steps, 0.0, max_time, 
			n_zSteps + 1, -margin, depth+margin);


	hnoconv = new TH1D("hnoconv","Ramo current",n_time_steps, 0.0, max_time);
	hconv   = new TH1D("hconv","Amplifier convoluted",n_time_steps, 0.0, max_time);

	// Convert Z to milimeters
	std::vector<double> z_chifs(n_zSteps+1);
	z_chifs = z_shifts;
	std::transform(z_chifs.begin(), z_chifs.end(), z_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));

	// filename for data analysis
	std::string hetct_filename = "NOirrad_"+dtime+"ps_"+cap+"pF_t"+trap+"ns_"+stepZ+"um_"+stepV+"V_"+neigh+"nns_"+scanType+".hetct";
	
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_filename, detector, C, dt, z_chifs, waveLength, scanType, file_carriers, voltages);

	//Loop on voltages
	
	for (int k = 0; k < n_vSteps + 1; k++) 
	{
		detector.solve_w_u();
		detector.set_voltages(voltages[k], v_depletion);
		detector.solve_d_u();
		detector.solve_w_f_grad();
		detector.solve_d_f_grad();
		detector.get_mesh()->bounding_box_tree();

		
		
		TH2D *i_rc ;
		       
	        // Loop on depth
		for (int i = 0; i < n_zSteps + 1; i++) 
		{
			std::cout << "Height " << z_shifts[i] << " of " << z_shifts.back() << " || Voltage " << voltages[k] << " of " << voltages.back() << std::endl;
			
			i_total= 0;
	

			// Launch all threads but one
			for (int thrID = 0; thrID < nThreads-1; thrID++) 
			{
				vva_hole[thrID] = std::vector<double>(n_time_steps, 0);
				vva_elec[thrID] = std::vector<double>(n_time_steps, 0);

				// Simulate the drift of both electrons and holes
				// for loop over threads, each thread simulates the corresponding carrier_list (thr_id)

				//TODO OOOOOOO
				//t[thrID]= std::thread(call_from_thread, std::ref(carrier_collection), dt, max_time, x_shift, z_shifts[i], i_elec[thrID], i_hole[thrID], thrID); 
				t[thrID]= std::thread(call_from_thread, std::ref(*carrier_collection), dt, max_time, x_shift, z_shifts[i], vva_elec[thrID], vva_hole[thrID], thrID); 

				// input should now include thr_id and only a 1-D member of i_whatevercarrier
			}

			vva_hole[nThreads-1] = std::vector<double>(n_time_steps, 0);
			vva_elec[nThreads-1] = std::vector<double>(n_time_steps, 0);
			//Launch the last thread AKA main thread
			
			std::valarray<double> i_elec((size_t) n_time_steps);	
			std::valarray<double> i_hole((size_t) n_time_steps);
			carrier_collection->simulate_drift( dt, max_time, x_shift, z_shifts[i], i_elec, i_hole, nThreads-1);
			vva_elec[nThreads-1].assign(std::begin(i_elec), std::end(i_elec));
			vva_hole[nThreads-1].assign(std::begin(i_hole), std::end(i_hole));

			

			// Join all threads
			for (int id = 0; id < nThreads-1; id++)
			{
			t[id].join();	//join all threads
			}

			// calculate totalcurrent
			for (int id = 0; id < nThreads; id++)
			{
				for (int tPos = 0; tPos < n_time_steps; tPos++) 
				{
					i_total[tPos] = i_total[tPos] + (vva_hole[id][tPos] + vva_elec[id][tPos]);
				}
			}

//									QVector<double> x_elec(n_time_steps), y_elec(n_time_steps);
//									QVector<double> x_hole(n_time_steps), y_hole(n_time_steps);
//									QVector<double> x_total(n_time_steps), y_total(n_time_steps);


									// Compute time + format vectors for writting to file
									for (int j=0; j < n_time_steps; j++)
									{
//										x_elec[j] = j*dt;
//										x_hole[j] = j*dt;
//										x_total[j] = j*dt;
//										y_elec[j] = i_elec[j];
//										y_hole[j] = i_hole[j];
//										y_total[j] = i_total[j];
			i_ramo.SetBinContent(j+1,i+1, i_total[j] );
			hnoconv->SetBinContent( j+1 , i_total[j] );
											}

			hconv = H1DConvolution( hnoconv , C*1.e12 );
		        if (i==0)
				{
					i_rc = new TH2D("i_rc", "i_rc", hconv->GetNbinsX(), hconv->GetXaxis()->GetXmin(),hconv->GetXaxis()->GetXmax() , n_zSteps + 1, -margin, depth+margin);
				}
			

			

#ifdef RCONLY
			QVector<double> y_shaped(n_time_steps);
			y_shaped[0]=y_total[0];
			i_rc->SetBinContent(1, i+1, y_shaped[0]);
			double alfa = dt/(RC+dt);

			for (int j = 1; j <n_time_steps; j++) 
			{
				y_shaped[j]=y_shaped[j-1]+alfa*(y_total[j]-y_shaped[j-1]);
				i_rc->SetBinContent(j+1,i+1,y_shaped[j]);
			}
#endif

			for (int j = 1; j <=hconv->GetNbinsX(); j++)
			{
				i_rc->SetBinContent(j, i+1 , hconv->GetBinContent(j) );
			}



											// save results
											// REDUNDANT STEP??
											//
											// TODO: Check
											//
//											QVector<QVector<double>> raw_results;
//											raw_results.resize(4);
//											raw_results[0] = x_total;
//											raw_results[1] = y_total;
//											raw_results[2] = y_elec;
//											raw_results[3] = y_hole;

			// Write file from TH1D
			utilities::write_to_file_row(hetct_filename, hnoconv, detector.get_temperature(), z_shifts[i], voltages[k]);
		}

		std::string root_filename = "NOirrad_"+dtime+"ps_"+cap+"pF_t"+trap+"ns_"+voltage+"V_"+neigh+"nns_"+scanType+".root";

		// Open a ROOT file to save result
		TFile *tfile = new TFile(root_filename.c_str(), "RECREATE" );
		i_ramo.Write();
		i_rc->Write();
		tfile->Close();
		
		delete i_rc;

	} // End for loop over voltages
        	
	delete carrier_collection;
	return 0;
}
