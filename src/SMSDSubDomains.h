#include <dolfin.h>

using namespace dolfin;

class CentralStripBoundary: public SubDomain {
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
  public:
    CentralStripBoundary(double pitch, double width, int nns );
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class NeighbourStripBoundary: public SubDomain {
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
  public:
    NeighbourStripBoundary(double pitch, double width, int nns );
    bool inside(const Array<double>& x, bool on_boundary) const;
};


class BackPlaneBoundary: public SubDomain {
  private:
    double _x_min;      // min x value
    double _x_max;      // max x value
    double _depth;      // detector depth
  public:
    BackPlaneBoundary(double x_min, double x_max, double depth);
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class PeriodicLateralBoundary: public SubDomain {
  private:
    double _x_min;      // min x value
    double _x_max;      // max x value
    double _depth;      // detector depth
  public:
    PeriodicLateralBoundary(double x_min, double x_max, double depth);
    bool inside(const Array<double>& x, bool on_boundary) const;
    void map(const Array<double>& x, Array<double>& y) const;
};
