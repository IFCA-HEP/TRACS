#include <SMSDetector.h>

SMSDetector::SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type, int n_cells_x = 100, int n_cells_y = 100, double tempK) :
    _pitch(pitch), //Distance between implants
    _width(width), //Size of the implant
    _depth(depth), //Vertical size of the pad (typically 300microns)
    _tempK(tempK),
    _nns(nns),
    _bulk_type(bulk_type), //Dopant type of the silicon (p/n)
    _implant_type(implant_type), //Dopant type of the implant, normally opposite of the bulk (n/p)
    _x_min(0.0), // Starting horizontal position for carrier generation (hereafter CG)
    _x_max(_pitch * (2*_nns+1)), // Endingvertical positio for CG
    _y_min(0.0), // Starting vertical position for CG (microns)
    _y_max(_depth), // Ending vertical position for CG (microns)
    _n_cells_x(n_cells_x),
    _n_cells_y(n_cells_y),
    _mesh(_x_min,_y_min,_x_max,_y_max, _n_cells_x, _n_cells_y),
    _periodic_boundary(_x_min, _x_max, _depth),
    _central_strip(_pitch, _width, _nns),
    _neighbour_strips(_pitch, _width, _nns),
    _backplane(_x_min, _x_max, _depth),
    _V_p(_mesh, _periodic_boundary),
    _a_p(_V_p, _V_p),
    _L_p(_V_p),
    _V_g(_mesh),
    _a_g(_V_g, _V_g),
    _L_g(_V_g),
    _w_u(_V_p),
    _d_u(_V_p),
    _w_f_grad(_V_g), // Weighting field
    _d_f_grad(_V_g)
{
}

void SMSDetector::set_voltages(double v_bias, double v_depletion)
{
  _v_strips = (_implant_type == 'n') ? v_bias : 0.0;
  _v_backplane = (_implant_type == 'p') ? v_bias : 0.0;
  _f_poisson = ((_bulk_type== 'p') ? +1.0 : -1.0)*(-2.0*v_depletion)/(_depth*_depth);
}

/*
 * Method for solving the weighting potential using Laplace equations 
 *
 *
 * 
 */
void SMSDetector::solve_w_u()
{

  // Solving Laplace equation f = 0
  Constant f(0);
  _L_p.f = f;

  // Set BC values
  Constant central_strip_V(1.0);
  Constant neighbour_strip_V(0.0);
  Constant backplane_V(0.0);
  // Set BC variables
  DirichletBC central_strip_BC(_V_p, central_strip_V, _central_strip);
  DirichletBC neighbour_strip_BC(_V_p, neighbour_strip_V, _neighbour_strips);
  DirichletBC backplane_BC(_V_p, backplane_V, _backplane);
  // Collect them
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&central_strip_BC);
  bcs.push_back(&neighbour_strip_BC);
  bcs.push_back(&backplane_BC);

  solve(_a_p == _L_p , _w_u, bcs);
}

/*
 * Getter for the weighting potential
 *
 *
 */
Function * SMSDetector::get_w_u()
{
  return &_w_u;
}

/*
 *
 * Method for solving the XXXXXXXXX d_u XXXXXXXXXXX using Poisson's equation
 * 
 *
 */

void SMSDetector::solve_d_u()
{
  // Solving Poisson equation for the d_u
  Constant f(_f_poisson);
  _L_p.f = f;

  // Set BC values
  Constant central_strip_V(_v_strips);
  Constant neighbour_strip_V(_v_strips);
  Constant backplane_V(_v_backplane);
  // Set BC variables
  DirichletBC central_strip_BC(_V_p, central_strip_V, _central_strip);
  DirichletBC neighbour_strip_BC(_V_p, neighbour_strip_V, _neighbour_strips);
  DirichletBC backplane_BC(_V_p, backplane_V, _backplane);
  // Collect them
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&central_strip_BC);
  bcs.push_back(&neighbour_strip_BC);
  bcs.push_back(&backplane_BC);

  solve(_a_p == _L_p , _d_u, bcs);
}

/*
 * Getter for the XXXXXXXXXX d_u XXXXXXXXXXXX
 *
 *
 */
Function * SMSDetector::get_d_u()
{
  return &_d_u;
}

/*
 * Method that calculates the weighting field inside the detector
 *
 *
 */
void SMSDetector::solve_w_f_grad()
{

  _L_g.u = _w_u;
  solve(_a_g == _L_g, _w_f_grad);
  // Change sign E = - grad(u)
  _w_f_grad = _w_f_grad * (-1.0);
}

/*
 * Getter for the weighting field
 *
 *
 */
Function * SMSDetector::get_w_f_grad()
{
  return &_w_f_grad;
}

/*
 * Method for calculating the d_f_grad XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  field
 *
 *
 */
void SMSDetector::solve_d_f_grad()
{
  _L_g.u = _d_u;
  solve(_a_g == _L_g, _d_f_grad);
  // Change sign E = - grad(u)
  _d_f_grad = _d_f_grad * (-1.0);

}

/*
 * Getter for the d_f_grad XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX field
 *
 */
Function * SMSDetector::get_d_f_grad()
{
  return &_d_f_grad;
}


/*
 * Method that checks if the carrier is inside or outside
 * of the detectore volume.
 *
 */
bool SMSDetector::is_out(const std::array< double,2> &x)
{
  bool out = true;
  if ( (x[0] > _x_min) && (x[0] < _x_max) && (x[1] > _y_min) && (x[1] < _y_max))
  {
    out = false;
  }
  return out;
}

/*
 * Getter for the minimum X value
 *
 *
 */
double  SMSDetector::get_x_min()
{
  return _x_min;
}

/*
 * Getter for the maximum X value
 *
 *
 */
double  SMSDetector::get_x_max()
{
  return _x_max;
}
/*
 *
 * Getter method for the temperature of the diode
 *
 */
double  SMSDetector::get_temperature()
{
  return _tempK;
}

/*
 * Getter for the minimum Y value.
 *
 */
double  SMSDetector::get_y_min()
{
  return _y_min;
}


/*
 * Getter for the maximum Y value.
 *
 *
 */
double  SMSDetector::get_y_max()
{
  return _y_max;
}

/*
 * Setter for the distance between implants.
 *
 *
 */
void SMSDetector::set_pitch(double pitch)
{
  _pitch = pitch;
}

/*
 * Setter for the detector width.
 *
 *
 */
void SMSDetector::set_width(double width)
{
  _width = width;
}

/*
 * Setter for the depth of the implants.
 *
 *
 */
void SMSDetector::set_depth(double depth)
{
  _depth = depth;
}

/*
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 *
 */
void SMSDetector::set_nns(int nns)
{
  _nns = nns;
}

/*
 *
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 */
void SMSDetector::set_bulk_type(char bulk_type)
{
  _bulk_type = bulk_type;
}

/*
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 *
 */
void SMSDetector::set_implant_type(char implant_type)
{
  _implant_type = implant_type;
}

/*
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 */
void SMSDetector::set_n_cells_x(int n_cells_x)
{
  _n_cells_x = n_cells_x;
}

/*
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 *
 */
void SMSDetector::set_n_cells_y(int n_cells_y)
{
  _n_cells_y  = n_cells_y;
}

// void SMSDetector::set_derived() {
//     _x_min = 0.0 ;
//     _x_max = _pitch * (2*_nns+1);
//     _y_min = 0;
//     _y_max = _depth;
//     _mesh = RectangleMesh(_x_min,_y_min,_x_max,_y_max, _n_cells_x, _n_cells_y);
//     _periodic_boundary = PeriodicLateralBoundary(_x_min, _x_max, _depth);
//     _central_strip = CentralStripBoundary(_pitch, _width, _nns);
//     _neighbour_strips = NeighbourStripBoundary(_pitch, _width, _nns);
//     _backplane = BackPlaneBoundary(_x_min, _x_max, _depth);
//     _V_p = Poisson::FunctionSpace(_mesh, _periodic_boundary);
//     _a_p = Poisson::BilinearForm(_V_p, _V_p);
//     _L_p = Poisson::LinearForm(_V_p);
//     _V_g = Gradient::FunctionSpace(_mesh);
//     _a_g = Gradient::BilinearForm(_V_g, _V_g);
//     _L_g = Gradient::LinearForm(_V_g);
//     _w_u = Function(_V_p);
//     _d_u = Function(_V_p);
//     _w_f_grad = Function(_V_g);
//     _d_f_grad = Function(_V_g);
// }


SMSDetector::~SMSDetector()
{

}




