// Copyright (C) 2007-2008 Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2007-07-11
// Last changed: 2012-11-12
//
// This demo program solves Poisson's equation,
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with homogeneous Dirichlet boundary conditions
// at y = 0, 1 and periodic boundary conditions at x = 0, 1.

#include <SMSDetector.h>

int main()
{


  double pitch = 80.;
  double width = 30.;
  double depth = 200.;
  int nns = 3;

  int n_cells_x = 300;
  int n_cells_y = 300;

  char bulk_type = 'p';
  char implant_type = 'n';

  SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y);

  detector.solve_w_u();
  Function * w_u = detector.get_w_u();

  // Plot solution
  plot(*w_u,"Weighting Potential","auto");
  interactive();

  return 0;
}
