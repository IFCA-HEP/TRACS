#ifndef CARRIERTRANSPORT_H
#define CARRIERTRANSPORT_H


#include <dolfin.h>

#include <CarrierMobility.h>

using namespace dolfin;

class DriftTransport
{
  private:
    JacoboniMobility _mu;
    Function * _d_f_grad;
    int _sign;


  public:
    DriftTransport(char carrier_type, Function * d_f_grad, double givenT = 253.);
		DriftTransport();
    ~DriftTransport();
    void operator() ( const std::array< double,2> &x , std::array< double,2> &dxdt , const double /* t */ );

};

#endif // CARRIERTRANSPORT_H
