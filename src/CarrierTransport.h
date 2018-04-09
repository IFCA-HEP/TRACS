
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

#ifndef CARRIERTRANSPORT_H
#define CARRIERTRANSPORT_H


#include <dolfin.h>

#include <CarrierMobility.h>

using namespace dolfin;

class DriftTransport
{
  private:
    JacoboniMobility _mu;
    Function * _d_f_grad;
    int _sign;


  public:
    DriftTransport(char carrier_type, Function * d_f_grad, double givenT = 253.);
		DriftTransport();
    ~DriftTransport();
    void operator() ( const std::array< double,2> &x , std::array< double,2> &dxdt , const double /* t */ );

};

#endif // CARRIERTRANSPORT_H
