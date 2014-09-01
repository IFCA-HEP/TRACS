#ifndef CARRIER_COLLECTION_H
#define CARRIER_COLLECTION_H

#include <string>
#include <sstream>
#include <fstream>

#include <QString>

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
};




#endif // CARRIER_COLLECTION_H


