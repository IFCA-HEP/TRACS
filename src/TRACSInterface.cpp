#include "TRACSInterface.h"
#include <mutex>          // std::mutex
/*
 * Constructor of class TRACSInterface
 * it takes the configuration file path as the only input and gets the 
 * rest of the data from said file
 *
 * This class provides a modular interface to TRACS simulator via 
 * different methods that allows the user to leverage all TRACS 
 * functionalities with their own code and use TRACS as a library for 
 * silicon detectors simulation.
 *
 */
 std::mutex mtx2;           // mutex for critical section
//extern const int num_threads;
TRACSInterface::TRACSInterface(std::string filename)
{
	neff_param = std::vector<double>(8,0);
	//utilities::parse_config_file(filename, carrierFile, depth, width, pitch, nns, temp, trapping, fluence, n_cells_x, n_cells_y, bulk_type, implant_type, C, dt, max_time, vBias, vDepletion, zPos, yPos, neff_param, neffType);
	utilities::parse_config_file(filename, carrierFile, depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y, bulk_type, implant_type, waveLength, scanType, C, dt, max_time, vInit, deltaV, vMax, vDepletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, neff_param, neffType);

	// Initialize vectors / n_Steps / detector / set default zPos, yPos, vBias / carrier_collection 
	if (fluence <= 0) // if no fluence -> no trapping
	{
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
		start = "NOirrad";
	} else 
	{
		start = "irrad";
	}

	n_zSteps = (int) std::floor((zMax-zInit)/deltaZ); // Simulation Steps
	n_zSteps1 = n_zSteps / 2;
	n_zSteps2 = (int) std::floor (n_zSteps - n_zSteps1);
	// if more threads than points
	if(num_threads>n_zSteps+1)
		{
			std::cout << "No. of threads > No. of z points! Program will terminate." << std::endl;	
			exit(EXIT_FAILURE);
		}
	n_zSteps_array = (int) std::floor ((n_zSteps+1) / num_threads);

	n_vSteps = (int) std::floor((vMax-vInit)/deltaV);
	n_ySteps = (int) std::floor((yMax-yInit)/deltaY);
	//WRITING TO FILES!
	// Convert relevant simulation numbers to string for fileNaming	
	dtime = std::to_string((int) std::floor(dt*1.e12));
	neigh = std::to_string(nns);
	stepV = std::to_string((int) std::floor(deltaV)); 
	stepZ = std::to_string((int) std::floor(deltaZ));
	stepY = std::to_string((int) std::floor(deltaY));
	cap = std::to_string((int) std::floor(C*1.e12));
	//std::string z_step  = std::to_string((int) std::floor(deltaZ));
	voltage = std::to_string((int) std::floor(vInit));

	parameters["allow_extrapolation"] = true;

	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	//SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	//detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	//detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType):
	//pDetector = &detector;

	n_tSteps = (int) std::floor(max_time / dt);

	carrierCollection = new CarrierCollection(detector);
	QString carrierFileName = QString::fromUtf8(carrierFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);

	//currents
	i_elec.resize((size_t) n_tSteps);
	i_hole.resize ((size_t) n_tSteps);
	i_total.resize((size_t) n_tSteps);

	i_elec = 0;
	i_hole = 0;
	i_total = 0;

	z_shifts.resize((size_t) n_zSteps+1,0.);
	z_shifts1.resize((size_t) n_zSteps1+1,0.);
	z_shifts2.resize((size_t) n_zSteps2+1,0.);
	z_shifts_array.resize(num_threads);
  	for (int i = 0; i < num_threads; i++)
    {	
    	if(i<(num_threads-1))
    		z_shifts_array[i].resize(n_zSteps_array, 0.);
		else
			z_shifts_array[i].resize((int)(n_zSteps+1-(n_zSteps_array)*(num_threads-1)), 0.);

	}
	
	y_shifts.resize ((size_t) n_ySteps+1,0.);
	voltages.resize((size_t) n_vSteps+1,0.);

	// Create voltages
	for (int i = 0; i < n_vSteps + 1; i++ ) 
	{
		voltages[i] = (i*deltaV)+vInit;
	}
	// CreatE shifts in Z
	for (int i = 0; i < n_zSteps + 1; i++ ) 
	{
		z_shifts[i] = (i*deltaZ)+zInit;	
	}
	/*
	int j = 0, k = 0;
	for (int i = 0; i < n_zSteps + 1; i++ ) 
	{
		z_shifts[i] = (i*deltaZ)+zInit;		
    	z_shifts_array[j][k] = z_shifts[i];
    	if(k<n_zSteps_array)
    	{
    		k++;
    	}
    	else
    	{
    		k = 0;
    		j++;
    	}

	}
	*/

	//"sampling" z array
	int l = 0;
	
		for (int i = 0; i < num_threads; i++)
		{
			for (int j = 0; j < z_shifts_array[i].size(); j++ ) 
			{
				z_shifts_array[i][j] = z_shifts[l];
				l++;
			}
		}
		
	

	// Create shifts in Y
	for (int i = 0; i < n_ySteps + 1; i++ ) 
	{
		y_shifts[i] = (i*deltaY)+yInit;
	}

	vBias = vInit;
	set_tcount(0);
	/*detector->set_voltages(vBias, vDepletion);
	detector->solve_w_u();
	detector->solve_d_u();
	detector->solve_w_f_grad();
	detector->solve_d_f_grad();
	detector->get_mesh()->bounding_box_tree();
	*/
	//i_ramo  = new TH1D("ramo","Ramo current",n_tSteps, 0.0, max_time);
	//i_conv   = new TH1D("conv","Amplifier convoluted",n_tSteps, 0.0, max_time);

	// Convert Z to milimeters
/*	std::vector<double> z_chifs(n_zSteps+1);
	z_chifs = z_shifts;
	std::transform(z_chifs.begin(), z_chifs.end(), z_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	
	// Convert Z to milimeters
	std::vector<double> y_chifs(n_ySteps+1);
	y_chifs = y_shifts;
	std::transform(y_chifs.begin(), y_chifs.end(), y_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
*/

	// filename for data analysis
/*	hetct_conv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_conv.hetct";
	hetct_noconv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_noconv.hetct";
	
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_conv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_noconv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
*/
	
	i_ramo  = NULL;
	i_rc    = NULL;
	i_conv  = NULL;
}
// Reads values, initializes detector

// Destructor
TRACSInterface::~TRACSInterface()
{
//	delete i_ramo, i_rc, i_conv;
}

/*
 * Convert i_total to TH1D
 */
TH1D * TRACSInterface::GetItRamo()
{
	if (i_ramo != NULL) 
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
   		//float r = f1->GetRandom();
		TString htit, hname;
		htit.Form("ramo_%d_%d", tcount, count1);
		hname.Form("Ramo_current_%d_%d", tcount, count1);
		i_ramo  = new TH1D(htit,hname,n_tSteps, 0.0, max_time);
		//std::cout << htit << std::endl;

		// Compute time + format vectors for writting to file
		for (int j=0; j < n_tSteps; j++)
		{
			i_ramo->SetBinContent(j+1, i_total[j] );
		}
		count1++;

	}
		return i_ramo;
}

/*
 * Convert i_total to TH1D after simulating simple RC circuit
 */
TH1D * TRACSInterface::GetItRc()
{

	if (i_rc != NULL) 
	{
	}
	else
	{	
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
   		//float r = f1->GetRandom();
		TString htit, hname;
		htit.Form("ramo_rc%d%d", tcount, count2);
		hname.Form("Ramo_current_%d_%d", tcount, count2);
		i_rc    = new TH1D(htit,hname,n_tSteps, 0.0, max_time);
		std::valarray<double> i_shaped ((size_t) n_tSteps);	
		double RC = 50.*C; // Ohms*Farad
		double alfa = dt/(RC+dt);

		for (int j = 1; j <n_tSteps; j++) 
		{
			i_shaped[j]=i_shaped[j-1]+alfa*(i_total[j]-i_shaped[j-1]);
			i_rc->SetBinContent(j+1, i_shaped[j]);
		}
		count2++;

	}
		return i_rc;
}

/*
 * Convert i_total to TH1D after convolution with the amplifier TransferFunction
 */
TH1D * TRACSInterface::GetItConv()
{
	if (i_conv != NULL) 
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
   		//float r = f1->GetRandom();
		TString htit, hname;
		htit.Form("ramo_conv_%d_%d", tcount, count3);
		hname.Form("Ramo current_%d_%d", tcount, count3);
		i_conv  = new TH1D(htit,hname,n_tSteps, 0.0, max_time);
		i_ramo = GetItRamo();
		//mtx2.lock();
		i_conv = H1DConvolution( i_ramo , C*1.e12, tcount );
		//mtx2.unlock();
		count3++;
	}
		return i_conv;
}

/*
 * Performs the simulation for all given carriers and stores the current in a valarray.
 * No variable is returned. To get the current one must choose the apropiate Getter for 
 * one's needs
 */
void TRACSInterface::simulate_ramo_current()
{
	i_rc = NULL;
	i_ramo = NULL;
	i_conv = NULL;
	i_hole = 0;
	i_elec = 0;
	i_total = 0;

	carrierCollection->simulate_drift( dt, max_time, yPos, zPos, i_elec, i_hole);

	//for (int tPos = 0; tPos < n_tSteps; tPos++) 
	//				{
	//					i_total[tPos] = i_elec[tPos] + i_hole[tPos];
	//				}
	//strange behaviour!?
	i_total = i_elec + i_hole;
}

/*
 * Sets the desired Neff parameters in the detector. Fields should be calculated 
 * again before simulating any current. Note that different neff parametrizations 
 * use different parameters so not all may be used at once.
 */
void TRACSInterface::set_NeffParam(std::vector<double> newParam)
{
	if ( newParam.size() == 8)
	{
		neff_param.assign(std::begin(newParam), std::end(newParam));
	}
	else
	{
		std::cout << "Error setting up new Neff, incorrect number of parameters" << std::endl;
	}

	detector->set_neff_param(neff_param);
}

/*
 * Sets the trapping time in the detector to the input value. 
 * Remember that the trapping time must be a positive number.
 * The smaller the trapping time the bigger the signal loss.
 */
void TRACSInterface::set_trappingTime(double newTrapTime)
{
	trapping = newTrapTime;
	detector->set_trapping_time(trapping);
}

/*
 * Sets how much the carriers will be displaced in the Z axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is relevant mainly for 
 * edge-TCT simulations.
 */
void TRACSInterface::set_zPos(double newZPos)
{
	zPos = newZPos;
}

/*
 * Sets how much the carriers will be displaced in the Y axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is used to center the red 
 * pulse in redTCT and to center the focus in TPA and edgeTCT
 */
void TRACSInterface::set_yPos(double newYPos)
{
	if (std::abs(newYPos) > (2*nns+1)*pitch)
	{
		std::cout << "Watch out! You probably set the laser out of the detector" << std::endl;
	}
	yPos = newYPos;
}

/*
 * Sets bias voltages in the detector, fields should be recalculated again 
 * before simulating any transients
 */
void TRACSInterface::set_vBias(double newVBias)
{
	vBias = newVBias;
	detector->set_voltages(vBias, vDepletion);
}

/*
 * Sets a number (current thread).
 * Used to index the different output files 
 *
 */
 void TRACSInterface::set_tcount(int tid)
{
	tcount = tid;
}

/*
 * Calculates the electric field and potential inside the detector. It is 
 * required after any modification of the Neff or the bias voltage applied. 
 * Weighting field and potential need not be calculated again since they 
 * are independent on those parameters.
 */
void TRACSInterface::calculate_fields()
{
	// Get detector ready
	//SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	//pDetector = &detector;
	detector->solve_w_u();
	detector->solve_d_u();
	detector->solve_w_f_grad();
	detector->solve_d_f_grad();
	detector->get_mesh()->bounding_box_tree();
	//detector->solve_d_f_grad();
	//detector->solve_d_u();
	
	
	
}

/*
 * Change parametrization of the Neff. Possibilities are: Trilinear (default)
 * Linear, Triconstant. More information on this three parametrizarions can
 * be found in the documentation (Config.TRACS and README.md)
 */
void TRACSInterface::set_neffType(std::string newParametrization)
{
	neffType = newParametrization;
	detector->set_neff_type(neffType);

}


/*
 * Allows the user to select a new carrier distribution from a file that
 * complies with the TRACS format. 
 */
void TRACSInterface::set_carrierFile(std::string newCarrFile)
{
	QString carrierFileName = QString::fromUtf8(newCarrFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);
}


/*
 * A loop through one parameter (depends on argument)
 *	"v": voltage, "z": z-axis, "y": y-axis
 *	example: TRACSsim->loop_on("v");
 * OVERLOADED!
 */
 void TRACSInterface::loop_on(std::string par)
 {
 	params[0] = 0; //zPos 
 	params[1] = 0; //yPos;
 	params[2] = 0; //vPos;
 	int p;
 	if (!par.compare("z"))
 	{
 		n_par0 = n_zSteps;
 		p = 0;
 	}
	 	else if (!par.compare("y"))
	 		{
	 		n_par0 = n_ySteps;
	 		p = 1;
	 		}
		 		else if (!par.compare("v"))
		 		{
		 			n_par0 = n_vSteps;
		 			p = 2;
		 		}
		 			else 
		 			{
		 				std::cout<<"Error in loop_on(string), invalid input argument! Use \"z\", \"y\", or \"v\"!" << std::endl;
		 				return;
		 			}
	//loop through the values and calculate
	for (params[p] = 0; params[p] < n_par0; params[p]++)
	{
		set_zPos(z_shifts[params[0]]);
		set_yPos(z_shifts[params[1]]);
		detector->set_voltages(voltages[params[2]], vDepletion);
		calculate_fields();
		simulate_ramo_current();
		GetItRamo();
		GetItRc();
		GetItConv();

	}
 	n_par0 = 0;
 }
/*
 * A loop through one parameter (depends on argument)
 *	"v": voltage, "z": z-axis, "y": y-axis
 *	example: TRACSsim->loop_on("v", "Z");
 */
void TRACSInterface::loop_on(std::string par1, std::string par2)
 {
 	params[0] = 0; //zPos 
 	params[1] = 0; //yPos;
 	params[2] = 0; //vPos;
 	char e1 = 0, e2 = 0, e3 = 0;
 	//v
 	if ((!par1.compare("v"))||((!par2.compare("v"))))
 	{
 		n_par2 = n_vSteps;
 		e1 = 1;
 	}
 	else
 		n_par2 = 1;
 	//y
 	if ((!par1.compare("y"))||((!par2.compare("y"))))
 	{
 		n_par1 = n_ySteps;
 		e2 = 1;
  	}
 	else
 		n_par1 = 1;
 	//z
 	if ((!par1.compare("z"))||((!par2.compare("z"))))
 	{
		n_par0 = n_zSteps;
  		e3 = 1;
 	}
 	else
 		n_par0 = 1;
 		 		
	if ((e1&e2)||(e1&e3)||(e2&e3))
		 	 	{
		 	 			//loop
		 		 	for (params[2] = 0; params[2] < n_par2; params[2]++)
		 			{
		 	 			detector->set_voltages(voltages[params[2]], vDepletion);
						calculate_fields();

						for (params[1] = 0; params[1] < n_par1; params[1]++)
						{
							set_yPos(y_shifts[params[1]]);

							for (params[0] = 0; params[0] < n_par0; params[0]++)
							{
								set_zPos(z_shifts[params[0]]);
								simulate_ramo_current();
								GetItRamo();
								GetItRc();
								GetItConv();
								std::cout<<"Success!" << std::endl;

							}
						}
	 		 		}
		 	 	}	
	else
		{
 			std::cout<<"Error in loop_on(str1, str2), invalid input arguments! Use \"z\", \"y\" and \"v\"!" << std::endl;
 		} 	 		
	 		 		
	 	
 		
	 	
 	n_par0 = 0;
 	n_par1 = 0;
 	n_par2 = 0;

 }

/*
 * A loop through all three parameters
 *	"v": voltage, "z": z-axis, "y": y-axis
 *	example: TRACSsim->loop_on("x","v","y");
 * 
 */
  void TRACSInterface::loop_on(std::string par1, std::string par2, std::string par3)
 {
 	params[0] = 0; //zPos 
 	params[1] = 0; //yPos;
 	params[2] = 0; //vPos;
 	char error = 1;
 	if ((!par1.compare("v"))||((!par2.compare("v")))||((!par3.compare("v"))))
 	{
 		 	if ((!par1.compare("y"))||((!par2.compare("y")))||((!par3.compare("y"))))
 		 	{
 		 		if ((!par1.compare("z"))||((!par2.compare("z")))||((!par3.compare("z"))))
		 		{
		 	 		n_par0 = n_zSteps;
		 	 		n_par1 = n_ySteps;
			 		n_par2 = n_vSteps;
			 		error = 0;
	 		 		//loop
		 		 	for (params[2] = 0; params[2] < n_par2 + 1; params[2]++)
		 			{
		 	 			detector->set_voltages(voltages[params[2]], vDepletion);
						calculate_fields();

						for (params[1] = 0; params[1] < n_par1 + 1; params[1]++)
						{
							set_yPos(y_shifts[params[1]]);
							for (params[0] = 0; params[0] < n_par0 + 1; params[0]++)
							{
								std::cout << "Height " << z_shifts[params[0]] << " of " << z_shifts.back()  <<  " || Y Position " << y_shifts[params[1]] << " of " << y_shifts.back() << " || Voltage " << voltages[params[2]] << " of " << voltages.back() << std::endl;								
								set_zPos(z_shifts[params[0]]);
								simulate_ramo_current();
								i_ramo = GetItRamo();
								//GetItRc();
								i_conv = GetItConv();
								//write to file
								//utilities::write_to_file_row(hetct_conv_filename, i_conv, detector->get_temperature(), y_shifts[params[1]], z_shifts[params[0]], voltages[params[2]]);
								//utilities::write_to_file_row(hetct_noconv_filename, i_ramo, detector->get_temperature(), y_shifts[params[1]], z_shifts[params[0]], voltages[params[2]]);
							}
						}
						std::string root_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_"+voltage+"V_"+neigh+"nns_"+scanType+".root";
						 // std::string hetct_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+".hetct";

						// Open a ROOT file to save result
						//TFile *tfile = new TFile(root_filename.c_str(), "RECREATE" );
						//i_ramo->Write();
						//i_rc->Write();
						//tfile->Close();
	 		 		}
	 		 	}
		 	}
 	}
 	if(error)
 		{
 			std::cout<<"Error in loop_on(str1, str2, str3), invalid input arguments! Use (\"z\", \"y\", \"v\")!" << std::endl;
 		}
	 	
 	n_par0 = 0;
 	n_par1 = 0;
 	n_par2 = 0;

 }


 /*
 * MULTITHREADING TEST!!! 
 *A loop through all three parameters
 *	"v": voltage, "z": z-axis, "y": y-axis
 *	example: TRACSsim->loop_on("x","v","y");
 * 
 */
  void TRACSInterface::loop_on(int tid)
 {
 	params[0] = 0; //zPos 
 	params[1] = 0; //yPos;
 	params[2] = 0; //vPos;
 	char error = 1;
 	//int n_par0, n_par1, n_par2;
 /*	std::vector<double> z_shifts_loc;
 	switch(tid)
 	{
 		case 0: n_par0 = n_zSteps1;
 		z_shifts_loc.resize((size_t) n_zSteps1+1,0.);
 		z_shifts_loc = z_shifts1;
 		std::cout <<  "z_shifts_1" << std::endl;
 		break;								

 		case 1: n_par0 = n_zSteps2;
 		z_shifts_loc.resize((size_t) n_zSteps2+1,0.);
 		z_shifts_loc = z_shifts2;
 		std::cout <<  "z_shifts_2" << std::endl;								
		break;

 		default: n_par0 = n_zSteps;
 		z_shifts_loc.resize((size_t) n_zSteps+1,0.);
 		z_shifts_loc = z_shifts;
 		std::cout <<  "z_shifts" << std::endl;								
 		break;


 	}
 	*/				n_par0 = (int) z_shifts_array[tid].size()-1;
		 	 		//n_par0 = n_zSteps_array;
		 	 		n_par1 = n_ySteps;
			 		n_par2 = n_vSteps;
			 		error = 0;
	 		 		//loop
		 		 	for (params[2] = 0; params[2] < n_par2 + 1; params[2]++)
		 			{
		 	 			detector->set_voltages(voltages[params[2]], vDepletion);
						calculate_fields();

						for (params[1] = 0; params[1] < n_par1 + 1; params[1]++)
						{
							set_yPos(y_shifts[params[1]]);
							for (params[0] = 0; params[0] < n_par0 + 1; params[0]++)
							{
								std::cout << "Height " << z_shifts_array[tid][params[0]] << " of " << z_shifts.back()  <<  " || Y Position " << y_shifts[params[1]] << " of " << y_shifts.back() << " || Voltage " << voltages[params[2]] << " of " << voltages.back() << std::endl;								
								set_zPos(z_shifts_array[tid][params[0]]);
								simulate_ramo_current();
								i_ramo = GetItRamo();
								i_ramo = NULL;
								i_rc = GetItRc();
								mtx2.lock();
								i_conv = GetItConv();
								mtx2.unlock();
								//write to file
								//mtx2.lock();
								utilities::write_to_file_row(hetct_conv_filename, i_conv, detector->get_temperature(), y_shifts[params[1]], z_shifts_array[tid][params[0]], voltages[params[2]]);
								utilities::write_to_file_row(hetct_noconv_filename, i_ramo, detector->get_temperature(), y_shifts[params[1]], z_shifts_array[tid][params[0]], voltages[params[2]]);
								//mtx2.unlock();
							}
						}
						//std::string root_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_"+voltage+"V_"+neigh+"nns_"+scanType+".root";
						 // std::string hetct_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+".hetct";

						// Open a ROOT file to save result
						//TFile *tfile = new TFile(root_filename.c_str(), "RECREATE" );
						//i_ramo->Write();
						//i_rc->Write();
						//tfile->Close();
	 		 		}
 	
	 		 	
 	/*if(error)
 		{
 			std::cout<<"Error in loop_on(str1, str2, str3), invalid input arguments! Use (\"z\", \"y\", \"v\")!" << std::endl;
 		}
	 */	
 	n_par0 = 0;
 	n_par1 = 0;
 	n_par2 = 0;

 }

 /*
  *
  * Write to file header. The input int is used to label files (multithreading)! 
  *
  *
  *
  */
  void TRACSInterface::write_header(int tid)
  {
  	// Convert Z to milimeters
	std::vector<double> z_chifs(n_zSteps+1);
	z_chifs = z_shifts;
	std::transform(z_chifs.begin(), z_chifs.end(), z_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	
	// Convert Z to milimeters
	std::vector<double> y_chifs(n_ySteps+1);
	y_chifs = y_shifts;
	std::transform(y_chifs.begin(), y_chifs.end(), y_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
  	// filename for data analysis
	hetct_conv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_"+std::to_string(tcount)+"_conv.hetct";
	hetct_noconv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_"+std::to_string(tcount)+"_noconv.hetct";
	
   
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_conv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_noconv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);

  }