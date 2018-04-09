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

#ifndef CARRIERMOBILITY_H
#define CARRIERMOBILITY_H

#include <cmath>

/*
 ****************JACOBONI MOBILITY**************
 *
 * This class is able to calculate the mobility 
 * of a charge carrier given its type, the 
 * detector thickness and the electric field at
 * de desired point.
 *
 *
 */

class JacoboniMobility
{
  private:
    double _T; // Temperature of the detector
    double _mu0; // Mobility
    double _vsat;
    double _beta;

  public:
    JacoboniMobility(char carrier_type, double T);
		JacoboniMobility();
    ~JacoboniMobility();
    double obtain_mobility(double e_field_mod);
};

#endif // CARRIERMOBILITY_H
