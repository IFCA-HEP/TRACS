#include "TRACSInterface.h"


TRACSInterface::TRACSInterface(std::string filename)
{
	neff_param = std::vector<double>(8,0);
	utilities::parse_config_file(filename, carrierFile, depth, width, pitch, nns, temp, trapping, fluence, n_cells_x, n_cells_y, bulk_type, implant_type, C, dt, max_time, vBias, vDepletion, zPos, yPos, neff_param, neffType);

	// Initialize vectors / n_Steps / detector / set default zPos, yPos, vBias / carrier_collection 
	if (fluence <= 0) // if no fluence -> no trapping
	{
		trapping = std::numeric_limits<double>::max();
	} else {}

	n_tSteps = (int) std::floor(max_time / dt);

	std::valarray<double> i_elec ((size_t) n_tSteps);	
	std::valarray<double> i_hole ((size_t) n_tSteps);
	std::valarray<double> i_total((size_t) n_tSteps);
	i_hole = 0;
	i_elec = 0;
	i_total = 0;

	parameters["allow_extrapolation"] = true;

	SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	pDetector = &detector;
	carrierCollection = new CarrierCollection(pDetector);
	QString carrierFileName = QString::fromUtf8(carrierFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);


	pDetector->set_voltages(vBias, vDepletion);
	pDetector->solve_w_u();
	pDetector->solve_d_u();
	pDetector->solve_w_f_grad();
	pDetector->solve_d_f_grad();
	pDetector->get_mesh()->bounding_box_tree();

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
		i_ramo  = new TH1D("hnoconv","Ramo current",n_tSteps, 0.0, max_time);
		// Compute time + format vectors for writting to file
		for (int j=0; j < n_tSteps; j++)
		{
			i_ramo->SetBinContent(j+1, i_total[j] );
		}
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
		i_rc    = new TH1D("hnoconv","Ramo current",n_tSteps, 0.0, max_time);
		std::valarray<double> i_shaped ((size_t) n_tSteps);	
		double RC = 50.*C; // Ohms*Farad
		double alfa = dt/(RC+dt);

		for (int j = 1; j <n_tSteps; j++) 
		{
			i_shaped[j]=i_shaped[j-1]+alfa*(i_total[j]-i_shaped[j-1]);
			i_rc->SetBinContent(j+1, i_shaped[j]);
		}

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
		i_conv  = new TH1D("hnoconv","Ramo current",n_tSteps, 0.0, max_time);
		i_ramo = GetItRamo();
		i_conv = H1DConvolution( i_ramo , C*1.e12 );
	}
		return i_conv;
}

/*
 * Performs the simulation for all given carriers and stores the current in a valarray.
 * No variable is returned. To get the current une must choose the apropiate Getter for 
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

	i_total = i_elec + i_hole;
}

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

	pDetector->set_neff_param(neff_param);
}

void TRACSInterface::set_trappingTime(double newTrapTime)
{
	if (newTrapTime > 0) 
	{
		trapping = newTrapTime;
	}
	else
	{
		std::cout << "Error setting up new Trapping Time, it must be a positive value" << std::endl;
	}

	pDetector->set_trapping_time(trapping);

}

void TRACSInterface::set_zPos(double newZPos)
{
	if (newZPos <= 0 && newZPos >= depth) 
	{
		std::cout << "Watch out! You probably set the laser out of the detector" << std::endl;
	}
	zPos = newZPos;
}

void TRACSInterface::set_yPos(double newYPos)
{
	if (std::abs(newYPos) > (2*nns+1)*pitch)
	{
		std::cout << "Watch out! You probably set the laser out of the detector" << std::endl;
	}
	yPos = newYPos;
}

void TRACSInterface::set_vBias(double newVBias)
{
	if (newVBias == vBias) 
	{
		std::cout << "The new Bias voltage is the same as the previous" << std::endl;
	}
	else
	{
		vBias = newVBias;
	}
	pDetector->set_voltages(vBias, vDepletion);
}

void TRACSInterface::calculate_fields()
{
	// Get detector ready
	pDetector->solve_d_f_grad();
	pDetector->solve_d_u();
}

void TRACSInterface::set_neffType(std::string newParametrization)
{
	neffType = newParametrization;
	pDetector->set_neff_type(neffType);

}


void TRACSInterface::set_carrierFile(std::string newCarrFile)
{
	QString carrierFileName = QString::fromUtf8(newCarrFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);
}
