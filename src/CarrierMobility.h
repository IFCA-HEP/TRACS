#include <cmath>

class JacoboniMobility
{
  private:
    double _T;
    double _mu0;
    double _vsat;
    double _beta;

  public:
    JacoboniMobility(char carrier_type, double T);
    ~JacoboniMobility();
    double obtain_mobility(double e_field_mod);
};
