#include <dolfin.h>

#include <CarrierMobility.h>


class DriftTransport
{
  private:
    JacoboniMobility _mu;
    Function * _d_f_grad;


  public:
    DriftTransport(char carrier_type, Function * d_f_grad);
    ~DriftTransport();
    void operator() ( const std::vector<double> &x , std::vector<double> &dxdt , const double /* t */ );

};
