#include <math.h>

/*
 * SOURCE TERM
 *
 * This class contains the source term for the poisson equation
 * Space charge distribution is to be parametrized here.
 *
 * At the moment only one space charge distribution is implemented 
 * but many more may come in the future together with a way for the 
 * user to select which one to use.
 */
 class Source : public Expression
 {
 public:

   /*
    * 3 ZONE space distribution
    *
    * This is the most general space charge distribution we expect
    * to use; it consists in 3 different straight lines corresponding 
    * to 3 different charge distributions
    */
	 
   void eval(Array<double>& values, const Array<double>& x) const
     {
	// Concentratio//n in transition and extremal points
	double y0 = -25.; // Neff(z0)
	double y1 = 0.02; // Neff(z1)
	double y2 = y1+y1*10; // Neff(z2)
	double y3 = 33; // Neff(z3)

	// Define diferent zones in depth
	double z0 = 0.;
	double z1 = 120.;
	double z2 = 220.;
	double z3 = 300.;

	// Charge concentration as a function of depthc (3zones)
	double neff_1 = ((y0-y1)/(z0-z1))*(x[1]-z0) + y0;
	double neff_2 = ((y1-y2)/(z1-z2))*(x[1]-z1) + y1;
	double neff_3 = ((y2-y3)/(z2-z3))*(x[1]-z2) + y2;

	// For continuity and smoothness purposes
	double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
	double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
	double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));

	double	neff = 0.5*(neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3);


/*
 * 1 ZONE approximatin
 *
 * First aproximation to the after-irradiation space charge distribution
 * Consists on a simple straight line with neff_max = 1e12
 *
 */


//	// Concentration in extremal points
//	double y0 = -12.; // Neff(z0)
//	double y3 = 13.; // Neff(z1)
//
//	// Define diferent zones in depth
//	double z0 = 0;
//	double z3 = 300;
//
//	// Function definition
//	double neff = ((y0-y3)/(z0-z3))*(x[1]-z0) + y0;

	// Fix problem with units from the PdC version
	values[0] = neff*0.00152132;

	
     }

 };
