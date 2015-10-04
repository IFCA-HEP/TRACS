#include "CarrierCollection.h"
// _carrier_list should be a N_thr-dimensional vector

/*
 * Comments
 *
 * Simulate_drift is an overloaded function 
 * 		WATCHOUT
 */

CarrierCollection::CarrierCollection(SMSDetector * detector) :
	_detector(detector)
{

}

void CarrierCollection::add_carriers_from_file(QString filename, int nThreads) // should get N_thr=1 as input
{

	// get char representation and make ifstream
	char * char_fn = filename.toLocal8Bit().data();
	std::ifstream infile(char_fn);
	// define temporal 1-D vector
	std::vector<Carrier> allCarriers;

	if (!infile.good()){
		std::cout<<"Warning: File not found"<<std::endl;
	}
	else
	{
		// process line by line
		std::string line;
		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			char carrier_type;
			double q, x_init, y_init, gen_time;
			if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) 
			{
				std::cout << "Error while reading file" << std::endl; 
				break;
			}

			// read all carriers
			Carrier carrier(carrier_type, q, x_init, y_init , _detector, gen_time);
			allCarriers.push_back(carrier);
		}

		// get #carriers, #carriers/N_thr, initialize carrier_collection
		int nCarriers = allCarriers.size();
		int carrierPerThread = (int) std::ceil(nCarriers/nThreads);
		Carrier emptyCarrier('e', 0.0, -10.0, -10.0, _detector, 0);
		int count = 0;
		std::vector<Carrier> defaultVector (carrierPerThread, emptyCarrier);
		// 2 for-loops (N_thr(#carriers/N_thr)) to fill each dimension 
		for (int i = 0; i < nThreads; i++)
		{
			_carrier_list.push_back(defaultVector);
			for (int j = 0; j < carrierPerThread; j++) 
			{
				if (count < nCarriers) 
				{
					_carrier_list[i][j] = allCarriers[count];
					count++;
				}
			}
		}
	}
}

// to thread or not to thread, that is the question (better thread-it)
void CarrierCollection::simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thrId)
{
	// range for through the carriers
	for (auto carrier : _carrier_list[thrId])
	{
		char carrier_type = carrier.get_carrier_type();
		// simulate drift and add to proper valarray
		if (carrier_type == 'e')
		{
			curr_elec += carrier.simulate_drift( dt , max_time);
		}
		else if (carrier_type =='h')
		{ 
			curr_hole += carrier.simulate_drift( dt , max_time);
		}
	}
	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);
	}
}

// input should also include thrId to simulate and store only 1-D vectors for each thread
void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x, double shift_y, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thrId)
{
	// range for through the carriers
	for (auto carrier : _carrier_list[thrId])
	{
		char carrier_type = carrier.get_carrier_type();

		// get and shift carrier position
		std::array< double,2> x = carrier.get_x();
		double x_init = x[0]+shift_x;
		double y_init = x[1]+shift_y;

		// simulate drift and add to proper valarray
		if (carrier_type == 'e')
		{
			curr_elec += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}
		else if (carrier_type =='h')
		{
			curr_hole += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}
	}
	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);
	}
}

//TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y,  TString hist_name, TString hist_title, int thrId)
//{
//  // get detector limits
//  double x_min = _detector->get_x_min();
//  double x_max = _detector->get_x_max();
//  double y_min = _detector->get_y_min();
//  double y_max = _detector->get_y_max();
//
//  // create histogram object
//  TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);
//
//   // range for through the carriers and fill the histogram
//  for (auto carrier : _carrier_list[thrId])
//  {
//    char carrier_type = carrier.get_carrier_type();
//    if (carrier_type == 'e')
//    {
//      std::array< double,2> x = carrier.get_x();
//      double q = carrier.get_q();
//      e_dist.Fill(x[0], x[1], q);
//    }
//  }
//  return e_dist;
//}

/*
 ********************** OVERLOADED FUNCTIONS FOR GUI COMPATIBILITY **************************
 */


void CarrierCollection::add_carriers_from_file(QString filename)
{
	// get char representation and make ifstream
	char * char_fn = filename.toLocal8Bit().data();
	std::ifstream infile(char_fn);

	// process line by line
	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		char carrier_type;
		double q, x_init, y_init, gen_time;
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) { break;}  //TODO show error?

		Carrier carrier(carrier_type, q, x_init, y_init , _detector, 1e-9);
		_carrier_list_sngl.push_back(carrier);
	}
}

void CarrierCollection::simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole)
{
	// range for through the carriers
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		// simulate drift and add to proper valarray
		if (carrier_type == 'e')
		{
			curr_elec += carrier.simulate_drift( dt , max_time);
		}
		else if (carrier_type =='h')
		{
			curr_hole += carrier.simulate_drift( dt , max_time);
		}
	}
	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);
	}
}

void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x, double shift_y, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole)
{
	// range for through the carriers
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		// simulate drift and add to proper valarray
		if (carrier_type == 'e')
		{
			// get and shift carrier position
			std::array< double,2> x = carrier.get_x();
			double x_init = x[0]+shift_x;
			double y_init = x[1]+shift_y;
			curr_elec += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}
		else if (carrier_type =='h')
		{
			// get and shift carrier position
			std::array< double,2> x = carrier.get_x();
			double x_init = x[0]+shift_x;
			double y_init = x[1]+shift_y;
			curr_hole += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}
	}
	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);
	}
}

TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y,  TString hist_name, TString hist_title)
{
	// get detector limits
	double x_min = _detector->get_x_min();
	double x_max = _detector->get_x_max();
	double y_min = _detector->get_y_min();
	double y_max = _detector->get_y_max();

	// create histogram object
	TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

	// range for through the carriers and fill the histogram
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		if (carrier_type == 'e')
		{
			std::array< double,2> x = carrier.get_x();
			double q = carrier.get_q();
			e_dist.Fill(x[0], x[1], q);
		}
	}
	return e_dist;
}


/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	COLLECTION **************************
 */
CarrierCollection::~CarrierCollection()
{

}
