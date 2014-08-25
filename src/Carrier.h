#ifndef CARRIER_H
#define CARRIER_H

#include  <valarray>

#include <CarrierTransport.h>
#include <SMSDetector.h>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

using namespace boost::numeric::odeint;


class Carrier
{
  private:
    char _carrier_type;
    double _q; // charge
    double _gen_time; // instant of generation of the carrier
    std::array< double,2> _x; // carrier position array
    std::array< double,2> _e_field; // electric field at the carrier position
    std::array< double,2> _w_field; // weighting field at the carrier positions
    double _e_field_mod;

    SMSDetector * _detector;
    DriftTransport _drift;
    JacoboniMobility _mu;

  public:
    Carrier( char carrier_type, double q, double x_init, double y_init, SMSDetector * detector, double gen_time);
    ~Carrier();

    std::valarray<double> simulate_drift( double dt, double max_time);
    std::valarray<double> simulate_drift(double dt, double max_time, double x_init, double y_init );
};

#endif // CARRIER_H
