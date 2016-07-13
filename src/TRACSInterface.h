#ifndef TRACSINTERFACE_H
#define TRACSINTERFACE_H
#include "SMSDetector.h"
#include "Source.h"
#include "utilities.h"
#include "Carrier.h"
#include "CarrierCollection.h"
#include <TH1D.h> // 1 Dimesional ROOT histogram 
#include <iterator>
#include <limits>  // std::numeric_limits
#include <cmath>
#include <functional>

extern TH1D *H1DConvolution( TH1D *htct, Double_t Cend=0. ) ; 

class TRACSInterface
{
	private:

		// Declaring external convolution function 
		double pitch; 
		double width; 
		double depth; 
		double temp; 
		double trapping; 
		double fluence; 
		double C; 
		double dt; 
		double max_time; 
		double vBias; 
		double vDepletion; 
		double zPos; 
		double yPos; 

		int nns; 
		int n_cells_y; 
		int n_cells_x; 
		int n_tSteps;

		char bulk_type; 
		char implant_type;

		std::vector<double> neff_param = {0};
		std::valarray<double> i_total;
		std::valarray<double> i_elec;
		std::valarray<double> i_hole;

		std::string carrierFile;
		std::string neffType;

		TH1D *i_ramo; 
		TH1D *i_rc; 
		TH1D *i_conv;

		// Pointer to detector and carrier collection
		SMSDetector * detector;//URBAN
		//SMSDetector * pDetector;
		CarrierCollection * carrierCollection;

	public:

		// Constructor
		TRACSInterface(std::string filename); // Reads values, initializes detector

		// Destructor
		~TRACSInterface();

		// Getters
		TH1D *GetItRamo();
		TH1D *GetItRc();
		TH1D *GetItConv();

		// Simulations
		void simulate_ramo_current();
		void calculate_fields();

		// Setters
		void set_NeffParam(std::vector<double> newParam);
		void set_trappingTime(double newTrapTime);
		void set_zPos(double newZPos);
		void set_yPos(double newYPos);
		void set_vBias(double newVBias);
		void set_neffType(std::string newParametrization);
		void set_carrierFile(std::string newCarrFile);

		
};

#endif // TRACSINTERFACE_H
