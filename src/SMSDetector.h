#include <dolfin.h>

#include "Poisson.h"
#include "Gradient.h"

#include <SMSDSubDomains.h>

using namespace dolfin;

class SMSDetector {
  private:
    // detector characteristics
    double _pitch; // in microns
    double _width; // in microns
    double _depth: // in microns
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
    int _n_cells_y:

    // bias
    double _v_strips;
    double _v_backplane;

    // potentials
    Function *_w_u;  // function to store the weighting potential
    Function *_d_u;  // function to store the drifting potential

    // fields
    Function *_w_f_grad; // function to store the gradient weighting field
    Function *_d_f_grad; // function to store the gradient weighting field

    // meshes (one for each could be used)
    RectangleMesh *_mesh; // mesh for both weighing and drifing potential

    // mesh subdomains
    PeriodicLateralBoundary *_periodic_boundary;
    CentralStripBoundary *_central_strip;
    NeighbourStripBoundary *_neighbour_strips;
    BackPlaneBoundary *_backplane;

    // Poisson PDE Function Space
    Poisson::FunctionSpace *_V_p;
    Poisson::BilinearForm *_a_p;
    Poisson::LinearForm *_L_p;

    // Gradient PDE Function Space
    Gradient::FunctionSpace *_V_g;
    Gradient::BilinearForm *_a_g;
    Gradient::LinearForm *_L_g;





  public:
    // default constructor and destructor
    SMSDetector(double pitch, double width, double depth,
                double nns, char bulk_type, char implant_type);
    ~SMSDetector();
    // set methods
    void set_mesh();
    void set_mesh_subdomains();
    void set_pde_variables();
    void solve_w_u();
    void solve_d_u();
    // get methods
    Function * get_w_u():

};
