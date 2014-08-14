#include <dolfin.h>
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
    Function *_w_u;  // Function to store the weighting potential
    Function *_d_u;  // Function to store the drifting field









  public:
    SMSDetector();
    ~SMSDetector();


};
