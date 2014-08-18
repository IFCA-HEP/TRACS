#include <SMSDSubDomains.h>

using namespace dolfin;


CentralStripBoundary::CentralStripBoundary(double pitch, double width, int nns )
{
  _pitch = pitch;
  _width = width;
  _nns = nns;
}

bool CentralStripBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] < DOLFIN_EPS ) && on_boundary) // y = 0 condition
  {
    double x_translated = x[0] - _pitch*_nns; // reference system change
    double l_lim = (_pitch - _width) / 2.0;
    double r_lim = l_lim + _width;
    if ((x_translated > l_lim*(1-DOLFIN_EPS)) && (x_translated < r_lim*(1+DOLFIN_EPS))) // check if strip
    {
      is_inside = true;
    }
  }
  return is_inside;
}

NeighbourStripBoundary::NeighbourStripBoundary(double pitch, double width, int nns )
{
  _pitch = pitch;
  _width = width;
  _nns = nns;
}

bool NeighbourStripBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] < DOLFIN_EPS ) && on_boundary) // y = 0 condition
  {
    int t_nns = 1 + 2*_nns; // total number of strips
    for ( int count = 0; count < t_nns; count++) // for loop for each strip
    {
      if (count != _nns) // not central strip
      {
        double x_translated = x[0] - _pitch*count; // reference system change
        double l_lim = (_pitch - _width) / 2.0;
        double r_lim = l_lim + _width;
        if ((x_translated > l_lim*(1-DOLFIN_EPS)) && (x_translated < r_lim*(1+DOLFIN_EPS))) // check if strip
        {
          is_inside = true;
        }
      }
    }
  }
  return is_inside;
}

BackPlaneBoundary::BackPlaneBoundary(double x_min, double x_max, double depth)
{
  _x_min = x_min;
  _x_max = x_max;
  _depth = depth;
}

bool BackPlaneBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] > (_depth - DOLFIN_EPS*_depth) ) && on_boundary) // y = depth condition
  {
    if ((x[0] > _x_min - DOLFIN_EPS) && (x[0] < _x_max + DOLFIN_EPS )) // within boundaries
    {
      is_inside = true;
    }
  }
  return is_inside;
}

PeriodicLateralBoundary::PeriodicLateralBoundary(double x_min, double x_max, double depth)
{
  _x_min = x_min;
  _x_max = x_max;
  _depth = depth;
}

// Left boundary is "target domain"
bool PeriodicLateralBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  return (std::abs(x[0]) < DOLFIN_EPS); // left lateral domain
}

// Map right boundary to left boundary
void PeriodicLateralBoundary::map(const Array<double>& x, Array<double>& y) const
{
  y[0] = x[0] - _x_max; //translate x coordinate
  y[1] = x[1];  // leave y equal
}








