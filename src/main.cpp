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

// Declaring external convolution function and threaded function
extern TH1D *H1DConvolution( TH1D *htct , Double_t Cend=0. ) ; 
void call_from_thread(CarrierCollection & cCollection, double dt, double max_time, double shift_x, double y_shifts, std::vector<double> curr_elec, std::vector<double> curr_hole, int thr_id);

//------------

// Threaded function to call carrier_collection
void call_from_thread(CarrierCollection & cCollection, double dt, double max_time, double shift_x, double y_shifts, std::vector<double> curr_elec, std::vector<double> curr_hole, int id)
{
	int nTimeSteps = curr_elec.size();
	std::valarray<double> i_elec((size_t) nTimeSteps);	
	std::valarray<double> i_hole((size_t) nTimeSteps);
	cCollection.simulate_drift( dt, max_time, shift_x, y_shifts, i_elec, i_hole, id);
	curr_elec.assign(std::begin(i_elec), std::end(i_elec));
	curr_hole.assign(std::begin(i_hole), std::end(i_hole));
}
/*
 ************** MAIN FUNCTION OF TRACS ***************
 *****************************************************
 ********** AKA: where the magic happens *************
 */

int main()
{
	// Declare variables with default values
	double pitch = 0,
		   width = 0,
		   depth = 0,
		   temp = 0,
		   trapping = 0,
		   fluence = 0,
		   C = 0,
		   dt = 0,
		   max_time = 0,
		   v_bias = 0,
		   vInit = 0,
		   deltaV = 0,
		   vMax = 0,
		   v_depletion = 0,
		   deltaZ = 0,
		   zInit = 0.,
		   zMax = depth,
		   yInit = 0.,
		   yMax = 0, //((2*nns)+1)*pitch,
		   deltaY = 5.;

	int nThreads = 0,
		nns = 0,
		n_cells_y = 0,
		n_cells_x = 0,
		waveLength = 0,
		n_vSteps = 0,
		n_zSteps = 0,
		n_ySteps = 0;

	char bulk_type = '\0', 
		 implant_type = '\0';

	std::string scanType = "defaultString";

	utilities::parse_config_file("Config.TRACS", depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y, bulk_type, implant_type, waveLength, scanType, C, dt, max_time, v_bias, vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY);
	
	// Create vector of (n-1) threads as the nth thread is the main thread
	std::thread t[nThreads-1];

	// Decide if there is trapping and make corresponding string
	std::string trap, start;
	if (fluence <= 0) // if no fluence -> no trapping
	{
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
		start = "NOirrad";
	}
	else
	{
		trap = std::to_string((int) std::floor(1.e9*trapping));
		start = "irrad";
	}

	TH1D *hnoconv , *hconv ;

	n_zSteps = (int) std::floor((zMax-zInit)/deltaZ); // Simulation Steps
	n_vSteps = (int) std::floor((vMax-vInit)/deltaV);
	n_ySteps = (int) std::floor((yMax-yInit)/deltaY);

	// Convert relevant simulation numbers to string for fileNaming	
	std::string dtime = std::to_string((int) std::floor(dt*1.e12));
	std::string neigh = std::to_string(nns);
	std::string stepV = std::to_string((int) std::floor(deltaV)); 
	std::string stepZ = std::to_string((int) std::floor(deltaZ));
	std::string stepY = std::to_string((int) std::floor(deltaY));
	std::string cap = std::to_string((int) std::floor(C*1.e12));
	//std::string z_step  = std::to_string((int) std::floor(deltaZ));
	std::string voltage = std::to_string((int) std::floor(v_bias));

	parameters["allow_extrapolation"] = true;

	SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence);

	detector.set_voltages(vInit, v_depletion);

	// Create carrier and observe movement
	SMSDetector * dec_pointer = &detector;

	// get number of steps from time
	int n_tSteps = (int) std::floor(max_time / dt);

	std::string file_carriers = "etct.carriers";
	QString filename = "etct.carriers";
	CarrierCollection * carrier_collection = new CarrierCollection(dec_pointer);

	// carrier_collection is now a #thr-dimensional vector
	carrier_collection->add_carriers_from_file(filename, nThreads); // input #threads

	// init arrays
	// These should are now #thr-dimensioal vectors
	std::valarray<double> i_total((size_t) n_tSteps);
	std::vector< std::vector<double> > vva_elec(nThreads, std::vector<double>(n_tSteps)); // vector of size N-threads where every element is a vector
	std::vector< std::vector<double> > vva_hole(nThreads, std::vector<double>(n_tSteps)); // same as above for holes current
	
	// Shifting  of Charge Carriers
	// the laser gets shifted in x and/or z direction depending on the arrays 
	// fed to the simulate_drift method inside the main for loop
	std::vector<double>  z_shifts(n_zSteps+1);
	std::vector<double>  y_shifts(n_ySteps+1) ; // laser shift in X axis to center laser focus over read-out strip
	std::vector<double>  voltages(n_vSteps+1);

	// Creat voltages
	for (int i = 0; i < n_vSteps + 1; i++ ) 
	{
		voltages[i] = (i*deltaV)+vInit;
	}
	// Creat shifts in Z
	for (int i = 0; i < n_zSteps + 1; i++ ) 
	{
		z_shifts[i] = (i*deltaZ)+zInit;
	}
	// Create shifts in Y
	for (int i = 0; i < n_ySteps + 1; i++ ) 
	{
		y_shifts[i] = (i*deltaY)+yInit;
	}
	TString hist_name = "i_ramo";
	TString hist_title = "i_ramo";
	TH2D i_ramo = TH2D(hist_name, hist_title,
			n_tSteps, 0.0, max_time, 
			n_zSteps + 1, zInit, zMax);

	hnoconv = new TH1D("hnoconv","Ramo current",n_tSteps, 0.0, max_time);
	hconv   = new TH1D("hconv","Amplifier convoluted",n_tSteps, 0.0, max_time);

	// Convert Z to milimeters
	std::vector<double> z_chifs(n_zSteps+1);
	z_chifs = z_shifts;
	std::transform(z_chifs.begin(), z_chifs.end(), z_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	
	// Convert Z to milimeters
	std::vector<double> y_chifs(n_ySteps+1);
	y_chifs = y_shifts;
	std::transform(y_chifs.begin(), y_chifs.end(), y_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));

	// filename for data analysis
	std::string hetct_conv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_conv.hetct";
	std::string hetct_noconv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_noconv.hetct";
	
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_conv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, file_carriers, voltages);
	utilities::write_to_hetct_header(hetct_noconv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, file_carriers, voltages);

	//Loop on voltages
	
	for (int k = 0; k < n_vSteps + 1; k++) 
	{
		detector.solve_w_u();
		detector.set_voltages(voltages[k], v_depletion);
		detector.solve_d_u();
		detector.solve_w_f_grad();
		detector.solve_d_f_grad();
		detector.get_mesh()->bounding_box_tree();
		
		// Loop on Y-axis
		for (int l = 0; l < n_ySteps + 1; l++) 
		{
			TH2D *i_rc ;

			// Loop on depth
			for (int i = 0; i < n_zSteps + 1; i++) 
			{
				std::cout << "Height " << z_shifts[i] << " of " << z_shifts.back()  <<  " || Y Position " << y_shifts[l] << " of " << y_shifts.back() << " || Voltage " << voltages[k] << " of " << voltages.back() << std::endl;

				i_total= 0;

				// Launch all threads but one
				for (int thrID = 0; thrID < nThreads-1; thrID++) 
				{
					vva_hole[thrID] = std::vector<double>(n_tSteps, 0);
					vva_elec[thrID] = std::vector<double>(n_tSteps, 0);

					// Simulate the drift of both electrons and holes
					// for loop over threads, each thread simulates the corresponding carrier_list (thr_id)
					t[thrID]= std::thread(call_from_thread, std::ref(*carrier_collection), dt, max_time, y_shifts[l], z_shifts[i], vva_elec[thrID], vva_hole[thrID], thrID); 

					// input should now include thr_id and only a 1-D member of i_whatevercarrier
				}
				vva_hole[nThreads-1] = std::vector<double>(n_tSteps, 0);
				vva_elec[nThreads-1] = std::vector<double>(n_tSteps, 0);

				//Launch the last thread AKA main thread
				std::valarray<double> i_elec((size_t) n_tSteps);	
				std::valarray<double> i_hole((size_t) n_tSteps);
				carrier_collection->simulate_drift( dt, max_time, y_shifts[l], z_shifts[i], i_elec, i_hole, nThreads-1);
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
					for (int tPos = 0; tPos < n_tSteps; tPos++) 
					{
						i_total[tPos] = i_total[tPos] + (vva_hole[id][tPos] + vva_elec[id][tPos]);
					}
				}
				// Compute time + format vectors for writting to file
				for (int j=0; j < n_tSteps; j++)
				{
					i_ramo.SetBinContent(j+1,i+1, i_total[j] );
					hnoconv->SetBinContent( j+1 , i_total[j] );
				}
				hconv = H1DConvolution( hnoconv , C*1.e12 );
				if (i==0)
				{
					i_rc = new TH2D("i_rc", "i_rc", hconv->GetNbinsX(), hconv->GetXaxis()->GetXmin(),hconv->GetXaxis()->GetXmax() , n_zSteps + 1, zInit, zMax);
				}
				for (int j = 1; j <=hconv->GetNbinsX(); j++)
				{
					i_rc->SetBinContent(j, i+1 , hconv->GetBinContent(j) );
				}
				// Write file from TH1D
				utilities::write_to_file_row(hetct_conv_filename, hconv, detector.get_temperature(), y_shifts[l], z_shifts[i], voltages[k]);
				utilities::write_to_file_row(hetct_noconv_filename, hnoconv, detector.get_temperature(), y_shifts[l], z_shifts[i], voltages[k]);
			}
			 std::string root_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_"+voltage+"V_"+neigh+"nns_"+scanType+".root";
			 // std::string hetct_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+".hetct";

			// Open a ROOT file to save result
			TFile *tfile = new TFile(root_filename.c_str(), "RECREATE" );
			i_ramo.Write();
			i_rc->Write();
			tfile->Close();

			delete i_rc;
		}
	}
	delete carrier_collection;
	return 0;
}
