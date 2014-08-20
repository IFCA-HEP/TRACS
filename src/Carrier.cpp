
#include <Carrier.h>

Carrier::Carrier( char carrier_type, double q,  double x_init, double y_init , SMSDetector * detector, double gen_time = 0.0 ) :
  _carrier_type(carrier_type),
  _q(q),
  _x(2),
  _gen_time(gen_time),
  _detector(detector),
  _drift(_carrier_type, _detector->get_d_f_grad()),
  _mu(_carrier_type, 300.)
{
  _x[0] = x_init;
  _x[1] = y_init;
}

std::vector<double> Carrier::simulate_drift(double dt, double max_time)
{
  int max_steps = (int) std::floor(max_time / dt);

  std::vector<double>  i_n(max_steps); // vector to save intensity

  runge_kutta4<std::vector<double>> stepper;

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
      std::cout << "Carrier Drift Time: " << t << " s" << std::endl;
      break;
    }
    else
    {
      std::vector< double > e_field(2);
      e_field[0] = ((*_detector->get_d_f_grad())[0])(_x[0],_x[1]);
      e_field[1] = ((*_detector->get_d_f_grad())[1])(_x[0],_x[1]);
      std::vector< double > w_field(2);
      w_field[0] = ((*_detector->get_w_f_grad())[0])(_x[0],_x[1]);
      w_field[1] = ((*_detector->get_w_f_grad())[1])(_x[0],_x[1]);
      double e_field_mod = sqrt(e_field[0]*e_field[0] + e_field[1]*e_field[1]);

      i_n[i] = _q * _mu.obtain_mobility(e_field_mod) * (e_field[0]*w_field[0] + e_field[1]*w_field[1]);
      stepper.do_step(_drift, _x, t, dt);
    }
    t+=dt;
  }
  return i_n;
}

Carrier::~Carrier()
{

}
