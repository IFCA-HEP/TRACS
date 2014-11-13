#include <CarrierMobility.h>

JacoboniMobility::JacoboniMobility( char carrier_type, double T = 300)
{
  // parameters from http://www.desy.de/~beckerj/cc/files/Thesis.pdf
  _T = T;
  if (carrier_type == 'e')
  {
    _mu0 = 1440.e8*std::pow(_T/300., -2.260);
    _vsat = 1.054e11  * std::pow(_T/300., -0.602);
    _beta = 0.992 * std::pow(_T/300., 0.572); // <100> orientation
  }
  else if (carrier_type == 'h')
  {
    _mu0 = 474.e8 * std::pow(_T/300., -2.619);
    _vsat = 0.940e11  * std::pow(_T/300., -0.226);
    _beta = 1.181 * std::pow(_T/300., 0.633 ); // <100> orientation
  }
}

double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return _mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta); // mum**2/ Vs
}

JacoboniMobility::~JacoboniMobility()
{

}
