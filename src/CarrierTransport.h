#include <dolfin.h>

#include <CarrierMobility.h>

using namespace dolfin;

class DriftTransport
{
  private:
    JacoboniMobility _mu;
    Function * _d_f_grad;


  public:
    DriftTransport(char carrier_type, Function * d_f_grad);
    ~DriftTransport();
    void operator() ( const std::array< double,2> &x , std::array< double,2> &dxdt , const double /* t */ );

};
