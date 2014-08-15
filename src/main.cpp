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

#include <dolfin.h>
#include "Poisson.h"
#include "Gradient.h"
#include <SMSDSubDomains.h>

#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>

using namespace dolfin;

int main()
{

  double pitch = 80.;
  double width = 30.;
  double depth = 200.;
  int nns = 3;

  double x_min = 0.;
  double x_max = pitch * (2*nns+1);
  double y_min = 0.;
  double y_max = depth;


  // Create rectangle mesh
  RectangleMesh mesh(x_min,y_min,x_max,depth, 150, 150);

  // Create periodic boundary condition
  PeriodicLateralBoundary periodic_boundary(x_min, x_max, depth);

  // Create functions
  Constant f(0);

  // Define PDE
  Poisson::FunctionSpace V(mesh, periodic_boundary);
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  L.f = f;

  // Create central strip BC
  Constant central_strip_V(1.0);
  CentralStripBoundary central_strip(pitch, width, nns);
  DirichletBC central_strip_BC(V, central_strip_V, central_strip);

  // Create neighbour strips BC
  Constant neighbour_strip_V(0.0);
  NeighbourStripBoundary neighbour_strip(pitch, width, nns);
  DirichletBC neighbour_strip_BC(V, neighbour_strip_V, neighbour_strip);

  // Create backplane BC
  Constant backplane_V(0.0);
  BackPlaneBoundary backplane(x_min, x_max, depth);
  DirichletBC backplane_BC(V, backplane_V, backplane);

  // Collect boundary conditions
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&central_strip_BC);
  bcs.push_back(&neighbour_strip_BC);
  bcs.push_back(&backplane_BC);

  // Compute solution
  Function u(V);
  solve(a == L, u, bcs);

  int n_bins_x = (int) pitch*(2*nns+1);
  int n_bins_y = (int) depth;


  TH2D *hist = new TH2D("h2","h2", n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

  for (int i=1; i<=n_bins_x; i++) {
    for (int j=1; j<=n_bins_y; j++) {
      double x_value = i - 0.5;
      double y_value = j - 0.5;
      hist->Fill(n_bins_x-x_value,n_bins_y-y_value, u(x_value, y_value));
    }
  }

  TFile file("hist.root", "recreate");
  hist->Write();
  file.Close();

  // Obtain gradient
  Gradient::FunctionSpace V_g(mesh);
  Gradient::BilinearForm a_g(V_g, V_g);
  Gradient::LinearForm L_g(V_g);
  L_g.u = u;

  Function grad_u(V_g);
  solve(a_g == L_g, grad_u);

  // Save solution in VTK format
  //File file_u("periodic.pvd");
  //file_u << u;

  // Plot solution
  plot(u,"Weighting Potential","auto");
  plot(grad_u[1],"Weighting Field","auto");
  interactive();

  return 0;
}
