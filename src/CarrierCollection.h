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

#ifndef CARRIER_COLLECTION_H
#define CARRIER_COLLECTION_H

#include <string>
#include <sstream>
#include <fstream>

#include <QString>

#include <TH2D.h>
#include <TString.h>

#include "Carrier.h"

/*
 ***********************************CARRIER COLLECTION***********************************
 *
 *
 *Detailed description
 *
 *
 *
 */

class CarrierCollection
{
  private:
    std::vector< std::vector<Carrier> > _carrier_list;
    std::vector<Carrier> _carrier_list_sngl;
    SMSDetector * _detector;

  public:
    CarrierCollection(SMSDetector * detector);
    ~CarrierCollection();

    void add_carriers_from_file(QString filename, int n_thr);
    void simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thr_id);
    void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thr_id);

//    TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist", int thr_id);

		/*
		 ********************* OVERLOADED FUNCTION FOR GUI COMPATIBILITY ****************************
		 */

    void add_carriers_from_file(QString filename);
    void simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole);
    void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec, std::valarray<double> &curr_hole);

    TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist");
    TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, double shift_x, double shift_y, TString hist_name = "e_dist", TString hist_title ="e_dist");


};




#endif // CARRIER_COLLECTION_H


