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

void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x, double shift_y, std::valarray<double> &curr_elec, std::valarray<double> &curr_hole)
{
  // range for through the carriers
  for (auto carrier : _carrier_list)
  {
    char carrier_type = carrier.get_carrier_type();
    // simulate drift and add to proper valarray
    if (carrier_type == 'e')
    {
      // get and shift carrier position
      std::array< double,2> x = carrier.get_x();
      double x_init = x[0]+shift_x;
      double y_init = x[1]+shift_y;
      curr_elec += carrier.simulate_drift( dt , max_time, x_init, y_init);
    }
    else if (carrier_type =='h')
    {
      // get and shift carrier position
      std::array< double,2> x = carrier.get_x();
      double x_init = x[0]+shift_x;
      double y_init = x[1]+shift_y;
      curr_hole += carrier.simulate_drift( dt , max_time, x_init, y_init);
    }
  }
}

TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y,  TString hist_name, TString hist_title)
{
  // get detector limits
  double x_min = _detector->get_x_min();
  double x_max = _detector->get_x_max();
  double y_min = _detector->get_y_min();
  double y_max = _detector->get_y_max();

  // create histogram object
  TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

   // range for through the carriers and fill the histogram
  for (auto carrier : _carrier_list)
  {
    char carrier_type = carrier.get_carrier_type();
    if (carrier_type == 'e')
    {
      std::array< double,2> x = carrier.get_x();
      double q = carrier.get_q();
      e_dist.Fill(x[0], x[1], q);
    }
  }
  return e_dist;
}

