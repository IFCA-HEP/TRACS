
#include <dolfin.h>

#include "Poisson.h"
#include "Gradient.h"

#include <SMSDSubDomains.h>

using namespace dolfin;

#ifndef SMSDetectorClass
#define SMSDetectorClass

class SMSDetector
{
  private:
    // detector characteristics
    double _pitch; // in microns
    double _width; // in microns
    double _depth; // in microns
    int _nns;
    char _bulk_type;
    char _implant_type;

    // some useful derived variables
    double _x_min; // in microns
    double _x_max; // in microns
    double _y_min; // in microns
    double _y_max; // in microns

    // Meshing parameters
    int _n_cells_x;
    int _n_cells_y;

    // bias
    double _v_strips;
    double _v_backplane;

    // poisson term for solving electric field
    double _f_poisson;

    // meshes (one for each could be used)
    RectangleMesh _mesh; // mesh for both weighing and drifing potential

    // mesh subdomains
    PeriodicLateralBoundary _periodic_boundary;
    CentralStripBoundary _central_strip;
    NeighbourStripBoundary _neighbour_strips;
    BackPlaneBoundary _backplane;

    // Poisson PDE Function Space
    Poisson::FunctionSpace _V_p;
    Poisson::BilinearForm _a_p;
    Poisson::LinearForm _L_p;

    // Gradient PDE Function Space
    Gradient::FunctionSpace _V_g;
    Gradient::BilinearForm _a_g;
    Gradient::LinearForm _L_g;

    // potentials
    Function _w_u;  // function to store the weighting potential
    Function _d_u;  // function to store the drifting potential

    // fields
    Function _w_f_grad; // function to store the weighting field (vectorial)
    Function _d_f_grad; // function to store the drifting field (vectorial)

  public:
    // default constructor and destructor
    SMSDetector(double pitch, double width, double depth,
                int nns, char bulk_type, char implant_type, int n_cells_x, int n_cells_y);
    ~SMSDetector();
    // set methods
    void set_voltages(double v_bias, double v_depletion);

    // solve potentials
    void solve_w_u();
    void solve_d_u();
    void solve_w_f_grad();
    void solve_d_f_grad();

    // get methods
    Function * get_w_u();
    Function * get_d_u();
    Function * get_w_f_grad();
    Function * get_d_f_grad();

    // some other methods
    bool is_out(std::vector<double> x);

};

#endif
