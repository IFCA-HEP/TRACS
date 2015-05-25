#include <CarrierMobility.h>

/*
 ******************************CARRIER MOBILITY*************************
 *
 * This function calculates the value of the mobility using data from 
 * desy experiment. The desired value is output as the return value of 
 * the "getter" method called: 'obtain_mobility' that requires the caller
 * to input the field at the desired point.
 *
 * Inputs: Carrier type ('e' for electrons / 'h' for holes)
 * 	   Temperature of the diode
 *
 * 	   - Field at desired point (only to obtain mobility values)
 *
 * "Outputs": Desired mobility.
 *
 *
 * Data Source:  http://www.desy.de/~beckerj/cc/files/Thesis.pdf
 *
 */

JacoboniMobility::JacoboniMobility( char carrier_type, double T)
{
  _T = T;
  if (carrier_type == 'e') // Electrons
  {
    _mu0 = 1440.e8*std::pow(_T/300., -2.260);
    _vsat = 1.054e11  * std::pow(_T/300., -0.602);
    _beta = 0.992 * std::pow(_T/300., 0.572); // <100> orientation
  }
  else if (carrier_type == 'h') // Holes
  {
    _mu0 = 474.e8 * std::pow(_T/300., -2.619);
    _vsat = 0.940e11  * std::pow(_T/300., -0.226);
    _beta = 1.181 * std::pow(_T/300., 0.633 ); // <100> orientation
  }
}

/*
 * OBTAIN MOBILITY - Getter method
 *
 * This method provides the value of the mobility using all the 
 * varibles initialized in the invokation.
 *
 */

double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return _mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta); // mum**2/ Vs
}

JacoboniMobility::~JacoboniMobility()
{

}
