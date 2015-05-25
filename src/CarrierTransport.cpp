#include <CarrierTransport.h>

DriftTransport::DriftTransport(char carrier_type, Function * d_f_grad, double givenT) :
  _mu(carrier_type, givenT)
 {
  _d_f_grad = d_f_grad;
  if (carrier_type == 'e') {
    _sign = -1;
  }
  else {
    _sign = 1;
  }
}

void DriftTransport::operator() ( const std::array<double,2>  &x , std::array<double,2>  &dxdt , const double /* t */ )
{
  // FIXME: avoid temporal creation overhead
  Array<double> e_field((std::size_t) 2); // temp wrap for e. field
  double e_field_mod;
  Point eval_point(x[0],x[1],0.0);
  (*_d_f_grad)(e_field, eval_point);
  e_field_mod = sqrt(e_field[0]*e_field[0] + e_field[1]*e_field[1]);
  dxdt[0] = _sign*_mu.obtain_mobility(e_field_mod) * e_field[0];
  dxdt[1] = _sign*_mu.obtain_mobility(e_field_mod) * e_field[1];
}

DriftTransport::~DriftTransport()
{

}
