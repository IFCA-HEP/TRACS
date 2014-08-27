#include <CarrierMobility.h>

JacoboniMobility::JacoboniMobility( char carrier_type, double T = 300)
{
  _T = T;
  if (carrier_type == 'e')
  {
    _mu0 = 1430. * std::pow(_T/300., -2.5);
    _vsat = 1.04e7  * std::pow(_T/300., -0.85);
    _beta = 1.01 * std::pow(_T/300., 0.87);
  }
  else if (carrier_type == 'h')
  {
    _mu0 = 480. * std::pow(_T/300., -2.82);
    _vsat = 9.2e6  * std::pow(_T/300., -0.04);
    _beta = 0.924 * std::pow(_T/300., 1.0 );
  }
}

double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return 100000000.*_mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta); //mum/s
}

JacoboniMobility::~JacoboniMobility()
{

}
