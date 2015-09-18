#include "SMSDetector.h"
#include "Source.h"
#include "utilities.h"
#include "Carrier.h"
#include "CarrierCollection.h"
#include <TH1D.h> // 1 Dimesional ROOT histogram 
#include <iterator>
#include <limits>  // std::numeric_limits
#include <functional>

class TRACSInterface
{
	private:

		// Declaring external convolution function and threaded function
		extern TH1D *H1DConvolution( TH1D *htct , Double_t Cend=0. ) ; 

		double pitch; 
		double width; 
		double depth; 
		double temp; 
		double trapping; 
		double fluence; 
		double C; 
		double dt; 
		double max_time; 
		double v_bias; 
		double vInit; 
		double deltaV; 
		double vMax; 
		double v_depletion; 
		double deltaZ; 
		double zInit; 
		double zMax = depth; 
		double yInit; 
		double yMax;

		unsigned int nns; 
		unsigned int n_cells_y; 
		unsigned int n_cells_x; 
		unsigned int waveLength; 
		unsigned int n_vSteps; 
		unsigned int n_zSteps; 
		unsigned int n_ySteps; 
		unsigned int n_tSteps;

		char bulk_type; 
		char implant_type;

		std::string scanType;
		std::vector<double> neff_param;

		TH1D *hnoconv; 
		TH1D *hconv; 
		TH1D *i_ramo; 
		TH1D *i_rc; 
		TH1D *i_conv;

		SMSDetector detector;

		// Create carrier and observe movement
		SMSDetector * dec_pointer;

		// get number of steps from time
		QString filename;
		CarrierCollection * carrier_collection;

	public:

		// Constructor
		TRACSInterface(std::string filename); // Reads values, initializes detector

		// Destructor
		~TRACSInterface();
		TH1D GetItRamo();
		TH1D GetItRc();
		TH1D GetItConv();
		void set_neffParam(std::vector<double> newParam);
		void set_trappingTime(double newTrapTime);
		void set_zPos(double newZPos);
		void set_yPos(double newYPos);
		void set_vBias(double newVBias);
		
		

}
