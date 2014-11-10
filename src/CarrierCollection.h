#ifndef CARRIER_COLLECTION_H
#define CARRIER_COLLECTION_H

#include <string>
#include <sstream>
#include <fstream>

#include <QString>

#include <TH2D.h>
#include <TString.h>

#include "Carrier.h"

class CarrierCollection
{
  private:
    std::vector<Carrier> _carrier_list;
    SMSDetector * _detector;

  public:
    CarrierCollection(SMSDetector * detector);
    ~CarrierCollection();

    void add_carriers_from_file(QString filename);
    void simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole);
    void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec, std::valarray<double> &curr_hole);

    TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist");

};




#endif // CARRIER_COLLECTION_H


