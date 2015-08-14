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
    SMSDetector * _detector;

  public:
    CarrierCollection(SMSDetector * detector);
    ~CarrierCollection();

    void add_carriers_from_file(QString filename, int n_thr = 1);
    void simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thr_id = 1);
    void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec, std::valarray<double> &curr_hole, int thr_id = 1);

    TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist", int thr_id = 1);

};




#endif // CARRIER_COLLECTION_H


