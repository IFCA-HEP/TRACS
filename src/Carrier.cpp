
#include <Carrier.h>

Carrier::Carrier( char carrier_type, double q,  double x_init, double y_init , SMSDetector * detector, double gen_time = 0.0 ) :
  _carrier_type(carrier_type),
  _q(q),
  _gen_time(gen_time),
  _detector(detector),
  _drift(_carrier_type, _detector->get_d_f_grad()),
  _mu(_carrier_type, 300.)
{
  _x[0] = x_init;
  _x[1] = y_init;
}

std::valarray<double> Carrier::simulate_drift(double dt, double max_time)
{
  // get number of steps from time
  int max_steps = (int) std::floor(max_time / dt);

  std::valarray<double>  i_n(max_steps); // valarray to save intensity

  runge_kutta4<std::array< double,2>> stepper;

  // wrapper for the arrays using dolphin array class
  Array<double> wrap_x(2, _x.data());
  Array<double> wrap_e_field(2, _e_field.data());
  Array<double> wrap_w_field(2, _w_field.data());


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
      _detector->get_d_f_grad()->eval(wrap_e_field, wrap_x);
      _detector->get_w_f_grad()->eval(wrap_w_field, wrap_x);
      _e_field_mod = sqrt(_e_field[0]*_e_field[0] + _e_field[1]*_e_field[1]);
      i_n[i] = _q * _mu.obtain_mobility(_e_field_mod) * (_e_field[0]*_w_field[0] + _e_field[1]*_w_field[1]);
      stepper.do_step(_drift, _x, t, dt);
    }
    t+=dt;
  }
  return i_n;
}

std::valarray<double> Carrier::simulate_drift(double dt, double max_time, double x_init, double y_init )
{
  _x[0] = x_init;
  _x[1] = y_init;

  // get number of steps from time
  int max_steps = (int) std::floor(max_time / dt);

  std::valarray<double>  i_n(max_steps); // valarray to save intensity

  runge_kutta4<std::array< double,2>> stepper;

  // wrapper for the arrays using dolphin array class
  Array<double> wrap_x(2, _x.data());
  Array<double> wrap_e_field(2, _e_field.data());
  Array<double> wrap_w_field(2, _w_field.data());

  if (_carrier_type == 'e') {
    _sign = -1;
  }
  else {
    _sign = 1;
  }


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
      _detector->get_d_f_grad()->eval(wrap_e_field, wrap_x);
      _detector->get_w_f_grad()->eval(wrap_w_field, wrap_x);
      _e_field_mod = sqrt(_e_field[0]*_e_field[0] + _e_field[1]*_e_field[1]);
      i_n[i] = _q * _sign * _mu.obtain_mobility(_e_field_mod) * (_e_field[0]*_w_field[0] + _e_field[1]*_w_field[1]);
      stepper.do_step(_drift, _x, t, dt);
    }
    t+=dt;
  }
  return i_n;
}

char Carrier::get_carrier_type()
{
  return _carrier_type;
}

std::array< double,2> Carrier::get_x()
{
  return _x;
}

double Carrier::get_q()
{
  return _q;
}


Carrier::~Carrier()
{

}
