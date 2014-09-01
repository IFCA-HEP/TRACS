#include "CarrierCollection.h"


CarrierCollection::CarrierCollection(SMSDetector * detector) :
  _detector(detector)
{

}


void CarrierCollection::add_carriers_from_file(QString filename)
{
  // get char representation and make ifstream
  char * char_fn = filename.toLocal8Bit().data();
  std::ifstream infile(char_fn);

  // process line by line
  std::string line;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    char carrier_type;
    double q, x_init, y_init, gen_time;
    if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) { break;}  //show error?

    Carrier carrier(carrier_type, q, x_init, y_init , _detector, 0.0);
    _carrier_list.push_back(carrier);
  }
}

void CarrierCollection::simulate_drift( double dt, double max_time, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole)
{
  // range for through the carriers
  for (auto carrier : _carrier_list)
  {
    char carrier_type = carrier.get_carrier_type();
    // simulate drift and add to proper valarray
    if (carrier_type == 'e')
    {
      curr_elec += carrier.simulate_drift( dt , max_time);
    }
    else if (carrier_type =='h')
    {
      curr_hole += carrier.simulate_drift( dt , max_time);
    }
  }
}

