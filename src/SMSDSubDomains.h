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

#ifndef SMSSUBDOMAINS_H
#define SMSSUBDOMAINS_H

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

#endif // SMSSUBDOMAINS_H
