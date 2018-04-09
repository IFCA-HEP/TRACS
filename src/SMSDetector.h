/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef SMSDETECTOR_H
#define SMSDETECTOR_H

#include <dolfin.h>
#include <cmath> 
#include <limits>  // std::numeric_limits

#include "Poisson.h"
#include "Gradient.h"

#include <SMSDSubDomains.h>

using namespace dolfin;

class SMSDetector
{
  private:
    // detector characteristics
    double _pitch; // in microns
    double _width; // in microns
    double _depth; // in microns
    double _tempK; // temperature in Kelvin
    double _trapping_time; // radiation damage effect
    double _fluence; // irradiation fluence (in neq)
    int _nns; // numbre of neighbouring strips
    char _bulk_type; // p or n
    char _implant_type; // n or p
	std::string _neff_type;
	std::vector<double> _neff_param; // Neff parametrization
		double _vdep; // depletion voltage

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
    SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type, int n_cells_x = 100, int n_cells_y = 100, double tempK = 253., double trapping = 9e300, double fluence = 0.0, std::vector<double> neff_param = {0}, std::string neff_type = "Trilinear");
    ~SMSDetector();
    // set methods
    void set_voltages(double v_bias, double v_depletion);
    void set_pitch(double pitch);
    void set_width(double width);
    void set_depth(double depth);
    void set_nns(int nns);
    void set_bulk_type(char bulk_type);
    void set_implant_type(char implant_type);
    void set_n_cells_x(int n_cells_x);
    void set_n_cells_y(int n_cells_y);
    void set_derived(); // Properly sets all derived quantities
    void set_temperature(double temperature);
    void set_trapping_time(double trapping_tau);
    void set_fluence(double fluencia);
	void set_neff_param(std::vector<double> neff_parameters);
	void set_neff_type(std::string newApproach);
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
	RectangleMesh * get_mesh();
    double get_x_min();
    double get_x_max();
    double get_y_min();
    double get_y_max();
    double get_temperature();
    double get_trapping_time();
	double get_fluence();
	double get_depth();
	double get_pitch();
	double get_width();
	int get_nns();
	char get_bulk_type();
	char get_implant_type();
	double get_vbias();
	double get_vdep();

    // some other methods
    bool is_out(const std::array< double,2> &x);

};

#endif // SMSDETECTOR_H
