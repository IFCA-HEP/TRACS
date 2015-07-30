#ifndef CARRIER_H
#define CARRIER_H

#include  <valarray>

#include <CarrierTransport.h>
#include <SMSDetector.h>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#endif

using namespace boost::numeric::odeint;

/*
 **************************CARRIER************************
 *
 *
 *  Detailed Description
 *
 *
 *
 *
 */

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
    int _sign; // sign to describe if carrier moves in e field direction or opposite

    SMSDetector * _detector;
    double _myTemp; // Temperature of the detector
    DriftTransport _drift;
    JacoboniMobility _mu;
    double _trapping_time;

  public:
    Carrier( char carrier_type, double q, double x_init, double y_init, SMSDetector * detector, double gen_time);
    ~Carrier();

    char get_carrier_type();
    std::array< double,2> get_x();
    double get_q();

    std::valarray<double> simulate_drift( double dt, double max_time);
    std::valarray<double> simulate_drift(double dt, double max_time, double x_init, double y_init );
};

#endif // CARRIER_H
