

#include <Carrier.h>

Carrier::Carrier( char carrier_type, double q, std::vector<double> x_o, SMSDetector * detector, double gen_time = 0.0 ) :
  _carrier_type(carrier_type),
  _q = q,
  _x = x_o,
  _gen_time = gen_time,
  _detector(detector),
  _drift(_carrier_type, _detector->get_d_f_grad()),
  _mu(_carrier_type)
{
}

Carrier::simulate_drift(double dt, double max_time)
{
  int max_steps = (int) std::floor(max_time / dt);

  std::vector<double>  i_n(max_steps); // vector to save intensity

  runge_kutta< std::vector<double> > stepper;

  double t=0.0;

  for ( int i = 0 ; i < max_steps; i++)
  {

    if (t < _gen_time)
    {
      i_n[i] = 0;
    }
    else if (_detector->is_out(_x))
    {
      i_n[i] = 0;
      break;
    }
    else
    {
      std::vector< double > e_field;
      e_field[0] = ((*_detector->get_d_f_grad())[0])(x[0],x[1]);
      e_field[1] = ((*_detector->get_d_f_grad())[1])(x[0],x[1]);
      std::vector< double > w_field;
      w_field[0] = ((*_detector->get_w_f_grad())[0])(x[0],x[1]);
      w_field[1] = ((*_detector->get_w_f_grad())[1])(x[0],x[1]);
      double e_field_mod = sqrt(e_field[0]*e_field[0] + e_field[1]*e_field[1]);

      i_n[i] = _q * _mu.obtain_mobility(e_field_mod) * (e_field[0]*w_field[0] + e_field[1]*w_field[1])
      stepper.do_step(_drift, _x, t, dt);
    }
    t+=dt;
  }
  return i_n;
}
