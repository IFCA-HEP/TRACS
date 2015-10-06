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

	// Concentration in transition and extremal points
	double y0 = -25.; // Neff(z0)
	double y1 = 0.02; // Neff(z1)
	double y2 = y1+y1*10; // Neff(z2)
	double y3 = 33; // Neff(z3)

	// Define diferent zones in depth
	double z0 = 0.;
	double z1 = 120.;
	double z2 = 220.;
	double z3 = 300.;
	std::string NeffApproach = "Trilinear";
	 
	void eval(Array<double>& values, const Array<double>& x) const
	{
		if (NeffApproach == "Triconstant") 
		{
			/*
			 * 3 ZONE constant space distribution
			 *
			 * We define here a Neff distribution consisting in 3 different zones
			 * each zone is defined as a constant within the given region.
			 * It uses all but the last parameter (y3 = Neff(z3)). It takes zX 
			 * values as boundaries of the zones and the three first yX as the 
			 * value of Neff in each region
			 *
			 * Even though a function like this is generally not continuous, we 
			 * add the hyperbolic tangent bridges to ensure not only continuity 
			 * but also derivability.
			 *
			 */
			double neff_1 = y0;
			double neff_2 = y1;
			double neff_3 = y2;

			// For continuity and smoothness purposes
			double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
			double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
			double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));

			double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
			values[0] = neff*0.00152132;
		}
		else if (NeffApproach == "Linear") 
		{
			/*
			 * 1 ZONE approximatin
			 *
			 * First aproximation to the after-irradiation space charge distribution
			 * Consists on a simple straight line defined by the points (z0, y0) and 
			 * (z3, y3) and neglects the rest of the values.
			 *
			 */

			double neff = ((y0-y3)/(z0-z3))*(x[1]-z0) + y0;
			values[0] = neff*0.00152132;
		}
		else 
		{
			/*
			 * 3 ZONE space distribution
			 *
			 * It consists in 3 different straight lines corresponding to 3 different
			 * charge distributions. It uses all 8 parameters to compute the Neff.
			 *
			 * Continuity is assumed as straight lines have common points, continuity 
			 * is ensured by the hyperbolic tangent bridges
			 */
			double neff_1 = ((y0-y1)/(z0-z1))*(x[1]-z0) + y0;
			double neff_2 = ((y1-y2)/(z1-z2))*(x[1]-z1) + y1;
			double neff_3 = ((y2-y3)/(z2-z3))*(x[1]-z2) + y2;

			// For continuity and smoothness purposes
			double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
			double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
			double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));

			double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
			values[0] = neff*0.00152132;
		}
		// Fix units from the PdC version


	}

	
	void set_NeffApproach(std::string Neff_type)
	{
		NeffApproach = Neff_type;
	}

	void set_y0(double newValue)
	{
		y0 = newValue;
	}

	void set_y1(double newValue)
	{
		y1 = newValue;
	}

	void set_y2(double newValue)
	{
		y2 = newValue;
	}

	void set_y3(double newValue)
	{
		y3 = newValue;
	}

	void set_z0(double newValue)
	{
		z0 = newValue;
	}

	void set_z1(double newValue)
	{
		z1 = newValue;
	}

	void set_z2(double newValue)
	{
		z2 = newValue;
	}

	void set_z3(double newValue)
	{
		z3 = newValue;
	}

 };
