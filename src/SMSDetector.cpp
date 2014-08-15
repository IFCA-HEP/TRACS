#include <SMSDetector.h>

SMSDetector::SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type, int n_cells_x = 100, int n_cells_y = 100 )
{
  // set members
  _pitch = pitch;
  _width = width;
  _depth = depth;
  _nns = nns;
  _bulk_type = bulk_type;
  _implant_type = implant_type;

  // set derived values
  _x_min = 0.;
  _x_max = _pitch * (2*_nns+1);
  _y_min = 0.;
  _y_max = _depth;

  // default number of cells
  _n_cells_x = n_cells_x;
  _n_cells_y = n_cells_y;

  // set mesh and subdomains
  set_mesh();
  set_mesh_subdomains();
  set_pde_variables();



}


void SMSDetector::set_mesh()
{
  // clean previous object and create new one
  //delete _mesh;
  _mesh = new RectangleMesh(_x_min,_y_min,_x_max,_depth, _n_cells_x, _n_cells_y);
}

void SMSDetector::set_mesh_subdomains()
{
  // clean previous objects and create new ones
  //delete _periodic_boundary;
  //delete _central_strip;
  //delete _neighbour_strips;
  //delete _backplane;
  _periodic_boundary = new PeriodicLateralBoundary(_x_min, _x_max, _depth);
  _central_strip = new CentralStripBoundary(_pitch, _width, _nns);
  _neighbour_strips = new NeighbourStripBoundary(_pitch, _width, _nns);
  _backplane = new BackPlaneBoundary(_x_min, _x_max, _depth);
}

void SMSDetector::set_pde_variables()
{
  // clean previous objects
  //delete _V_p;
  //delete _a_p;
  //delete _L_p;
  //delete _V_g;
  //delete _a_g;
  //delete _L_g;
  //delete _w_u;
  //delete _d_u;
  //delete _w_f_grad;
  //delete _d_f_grad;
  // create new objects according to the new mesh
  _V_p = new Poisson::FunctionSpace(*_mesh, *_periodic_boundary);
  _a_p = new Poisson::BilinearForm(*_V_p, *_V_p);
  _L_p = new Poisson::LinearForm(*_V_p);
  _V_g = new Gradient::FunctionSpace(*_mesh);
  _a_g = new Gradient::BilinearForm(*_V_g, *_V_g);
  _L_g = new Gradient::LinearForm(*_V_g);
  _w_u = new Function(*_V_p);
  _d_u = new Function(*_V_p);
  _w_f_grad = new Function(*_V_g);
  _d_f_grad = new Function(*_V_g);
}

void SMSDetector::solve_w_u()
{

  // Solving Laplace equation f = 0
  Constant f(0);
  (*_L_p).f = f;

  // Set BC values
  Constant central_strip_V(1.0);
  Constant neighbour_strip_V(0.0);
  Constant backplane_V(0.0);
  // Set BC variables
  DirichletBC central_strip_BC(*_V_p, central_strip_V, *_central_strip);
  DirichletBC neighbour_strip_BC(*_V_p, neighbour_strip_V, *_neighbour_strips);
  DirichletBC backplane_BC(*_V_p, backplane_V, *_backplane);
  // Collect them
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&central_strip_BC);
  bcs.push_back(&neighbour_strip_BC);
  bcs.push_back(&backplane_BC);

  solve(*_a_p == *_L_p, *_w_u, bcs);

}

Function * SMSDetector::get_w_u()
{
  return _w_u;
}

SMSDetector::~SMSDetector()
{

}




