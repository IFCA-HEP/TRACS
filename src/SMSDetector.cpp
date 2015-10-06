#include <SMSDetector.h>
#include <dolfin.h>
#include <Source.h>

SMSDetector::SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type, int n_cells_x, int n_cells_y, double tempK, double trapping, double fluence, std::vector<double> neff_param, std::string neff_type) :
    
    _pitch(pitch), //Distance between implants
    _width(width), //Size of the implant
    _depth(depth), //Vertical size of the pad (typically 300microns)
    _tempK(tempK), // Temperature of the detector
    _trapping_time(trapping), // Trapping constant simulates radiation-induced defects (traps)
    _fluence(fluence), // Irradiation fluence in neutron equivalent (neq)
    _nns(nns), // Number of neighbouring strips
    _bulk_type(bulk_type), //Dopant type of the silicon (p/n)
    _implant_type(implant_type), //Dopant type of the implant, normally opposite of the bulk (n/p)
	_neff_type(neff_type), // Select aproach to parametrize Neff (irradiation only)
    _neff_param(neff_param), // Parametrized description of the Space Charge distribution
    _x_min(0.0), // Starting horizontal position for carrier generation (hereafter CG)
    _x_max(_pitch * (2*_nns+1)), // Endingvertical positio for CG
    _y_min(0.0), // Starting vertical position for CG (microns)
    _y_max(_depth), // Ending vertical position for CG (microns)

    // Mesh properties
    _n_cells_x(n_cells_x),
    _n_cells_y(n_cells_y),
    _mesh(_x_min,_y_min,_x_max,_y_max, _n_cells_x, _n_cells_y),
    _periodic_boundary(_x_min, _x_max, _depth),

    // More detector properties/parts
    _central_strip(_pitch, _width, _nns),
    _neighbour_strips(_pitch, _width, _nns),
    _backplane(_x_min, _x_max, _depth), 

    // Functions & variables to solve the PDE
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

/*
 * Set right polarity for detector depending on the type (p-on-n/n-on-p)
 */
void SMSDetector::set_voltages(double v_bias, double v_depletion)
{
	_vdep = v_depletion; // Store depletion voltage
  _v_strips = (_implant_type == 'n') ? v_bias : 0.0;
  _v_backplane = (_implant_type == 'p') ? v_bias : 0.0;
  // neff defined in F/microns
  _f_poisson = ((_bulk_type== 'p') ? +1.0 : -1.0)*(-2.0*v_depletion)/(_depth*_depth);
}

/*
 * Method for solving the weighting potential using Laplace triangles 
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
 * Method for solving the XXXXXXXXX d_u XXXXXXXXXXX using Poisson's equation
 */

void SMSDetector::solve_d_u()
{
	Constant fpois(_f_poisson); 
	Source f;
	if (_fluence <= 0)
	{
		_L_p.f = fpois;
		// Idiot-proofing
		_trapping_time = std::numeric_limits<double>::max();
	}
	else 
	{
		f.set_NeffApproach(_neff_type);
		f.set_y0(_neff_param[0]);
		f.set_y1(_neff_param[1]);
		f.set_y2(_neff_param[2]);
		f.set_y3(_neff_param[3]);
		f.set_z0(_neff_param[4]);
		f.set_z1(_neff_param[5]);
		f.set_z2(_neff_param[6]);
		f.set_z3(_neff_param[7]);
		_L_p.f = f;
	}

	
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
 * Method that calculates the weighting field inside the detector
 */
void SMSDetector::solve_w_f_grad()
{

  _L_g.u = _w_u;
  solve(_a_g == _L_g, _w_f_grad);
  // Change sign E = - grad(u)
  _w_f_grad = _w_f_grad * (-1.0);
}

/*
 * Method for calculating the Electric field
 */
void SMSDetector::solve_d_f_grad()
{
  _L_g.u = _d_u;
  solve(_a_g == _L_g, _d_f_grad);
  // Change sign E = - grad(u)
  _d_f_grad = _d_f_grad * (-1.0);

}

/*
 * Method that checks if the carrier is inside or outside
 * of the detectore volume.
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

/************************************************************************
*************************************************************************
***                                                                   ***
***            BEHOLD!! BORING SETTERS AND GETTERS AHEAD              ***
***                                                                   ***
*************************************************************************
*************************************************************************/

/*
 * Getter for the weighting potential
 */
Function * SMSDetector::get_w_u()
{
  return &_w_u;
}

/*
 * Getter for the weighting field
 */
Function * SMSDetector::get_w_f_grad()
{
	return &_w_f_grad;
}

/*
 * Getter for the XXXXXXXXXX d_u XXXXXXXXXXXX
 */
Function * SMSDetector::get_d_u()
{

	return &_d_u;
}

/*
 * Getter for the drift field
 */
Function * SMSDetector::get_d_f_grad()
{
	return &_d_f_grad;
}

/*
 * Getter for the mesh of both fields
 */
RectangleMesh * SMSDetector::get_mesh()
{
	return &_mesh;
}

/*
 * Getter for the minimum X value
 */
double  SMSDetector::get_x_min()
{
  return _x_min;
}

/*
 * Getter for the maximum X value
 */
double  SMSDetector::get_x_max()
{
  return _x_max;
}
/*
 * Getter method for the temperature of the diode
 */
double  SMSDetector::get_temperature()
{
  return _tempK;
}

/*
 * Getter for the minimum Y value.
 */
double  SMSDetector::get_y_min()
{
  return _y_min;
}


/*
 * Getter for the maximum Y value.
 */
double SMSDetector::get_y_max()
{
  return _y_max;
}

/*
 * Getter for the trapping factor
 */
double SMSDetector::get_trapping_time()
{
  return _trapping_time;
}

/*
 * Getter for the fluence
 */
double SMSDetector::get_fluence()
{
  return _fluence;
}

/*
 * Getter for the pitch
 */
double SMSDetector::get_depth()
{
  return _depth;
}

/*
 * Getter for the pitch
 */
double SMSDetector::get_pitch()
{
  return _pitch;
}

/*
 * Getter for the width
 */
double SMSDetector::get_width()
{
  return _width;
}

/*
 * Getter for the number of neighbouring strips
 */
int SMSDetector::get_nns()
{
  return _nns;
}

/*
 * Getter for the bulk type
 */
char SMSDetector::get_bulk_type()
{
  return _bulk_type;
}

/*
 * Getter for the implant type
 */
char SMSDetector::get_implant_type()
{
  return _implant_type;
}

/*
 * Getter for the bias voltage
 */
double SMSDetector::get_vbias()
{
  return _v_strips-_v_backplane;
}


/*
 * Getter for the depletion voltage
 */
double SMSDetector::get_vdep()
{
  return _vdep;
}

///////////////////////SETTERS//////////////////////////////////

/*
 * Setter for the distance between implants.
 */
void SMSDetector::set_pitch(double pitch)
{
  _pitch = pitch;
}

/*
 * Setter for the detector width.
 */
void SMSDetector::set_width(double width)
{
  _width = width;
}

/*
 * Setter for the depth of the implants.
 */
void SMSDetector::set_depth(double depth)
{
  _depth = depth;
}

/*
 * Setter for the number of neighbouring strips
 */
void SMSDetector::set_nns(int nns)
{
  _nns = nns;
}

/*
 * Setter for changing the bulk type
 */
void SMSDetector::set_bulk_type(char bulk_type)
{
  _bulk_type = bulk_type;
}

/*
 * Setter for changing implant type
 */
void SMSDetector::set_implant_type(char implant_type)
{
  _implant_type = implant_type;
}

/*
 * Setter for changing the number of cells in x direction of the mesh
 */
void SMSDetector::set_n_cells_x(int n_cells_x)
{
  _n_cells_x = n_cells_x;
}

/*
 * Setter for changing the number of cells in y direction of the mesh
 */
void SMSDetector::set_n_cells_y(int n_cells_y)
{
  _n_cells_y  = n_cells_y;
}

/*
 * Setter for changing the temperature of the detector
 */
void SMSDetector::set_temperature(double temperature)
{
  _tempK = temperature;
}

/*
 * Setter for changing the number of cells in y direction of the mesh
 */
void SMSDetector::set_trapping_time(double trapping_tau)
{
  _trapping_time = trapping_tau;
}

/*
 * Setter for changing the number of cells in y direction of the mesh
 */
void SMSDetector::set_fluence(double fluencia)
{
  _fluence = fluencia;
}

/*
 * Setter for changing the space charge distribution
 */
void SMSDetector::set_neff_param(std::vector<double> neff_parameters)
{
  _neff_param[0] = neff_parameters[0]; // y0
  _neff_param[1] = neff_parameters[1]; // y1
  _neff_param[2] = neff_parameters[2]; // y2
  _neff_param[3] = neff_parameters[3]; // y3
  _neff_param[4] = 0.0; // z0 (Idiot proof)
  _neff_param[5] = neff_parameters[5]; // z1
  _neff_param[6] = neff_parameters[6]; // z2
  _neff_param[7] = _depth; // z3 (Idiot proof)
}

void SMSDetector::set_neff_type(std::string newApproach)
{
	_neff_type = newApproach;
}

SMSDetector::~SMSDetector()
{

}

