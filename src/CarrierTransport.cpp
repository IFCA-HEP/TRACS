#include <CarrierTransport.h>

DriftTransport::DriftTransport(char carrier_type, Function * d_f_grad) :
  _mu(carrier_type, 300.)
{
  _d_f_grad = d_f_grad;
}

void DriftTransport::operator() ( const std::vector<double> &x , std::vector<double> &dxdt , const double /* t */ )
{
  std::vector< double > e_field(2);
  e_field[0] = ((*_d_f_grad)[0])(x[0],x[1]);
  e_field[1] = ((*_d_f_grad)[1])(x[0],x[1]);
  double e_field_mod = sqrt(e_field[0]*e_field[0] + e_field[1]*e_field[1]);
  dxdt[0] = _mu.obtain_mobility(e_field_mod) * e_field[0];
  dxdt[1] = _mu.obtain_mobility(e_field_mod) * e_field[1];
}

DriftTransport::~DriftTransport()
{

}

