#ifndef CARRIERMOBILITY_H
#define CARRIERMOBILITY_H

#include <cmath>

/*
 ****************JACOBONI MOBILITY**************
 *
 * This class is able to calculate the mobility 
 * of a charge carrier given its type, the 
 * detector thickness and the electric field at
 * de desired point.
 *
 *
 */

class JacoboniMobility
{
  private:
    double _T; // Temperature of the detector
    double _mu0; // Mobility
    double _vsat;
    double _beta;

  public:
    JacoboniMobility(char carrier_type, double T);
    ~JacoboniMobility();
    double obtain_mobility(double e_field_mod);
};

#endif // CARRIERMOBILITY_H
