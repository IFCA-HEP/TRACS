#include <vector>

#include <CarrierTransport.h>


class Carrier
{
  private:
    char _carrier_type;
    double _q; // charge
    double _gen_time; // instant of generation of the carrier
    std::vector<double> _x; // position state vector

    SMSDetector * _detector;
    DriftTransport _drift;
    JacoboniMobility _mu;

  public:
    Carrier( char carrier_type, , double q, std::vector<double> x, SMSDetector * detector, double gen_time);
    ~Carrier();

    std::vector<double> simulate_drift( double dt, int max_steps);
};
