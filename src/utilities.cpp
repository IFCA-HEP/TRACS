#include "utilities.h" 

// function to export to a root histogram from a dolfin 1d function
TH2D utilities::export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max)
{
  TH2D hist = TH2D(hist_name,hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;
  for (int i=1; i<=n_bins_x; i++) {
    for (int j=1; j<=n_bins_y; j++) {
    double x_value = (i - 0.5)*step_x;
    double y_value = (j - 0.5)*step_y;
    hist.Fill(x_max-x_value,y_max-y_value, func(x_value, y_value));
    }
  }
  return hist;
}

// function to paint a TH2D root histogram in a color map
void  utilities::paint_TH2D_qcp(TH2D hist, QCPColorMap * color_map)
{
  // get some required variables
  int n_bins_x = hist.GetNbinsX();
  int n_bins_y = hist.GetNbinsY();
  double x_min = hist.GetXaxis()->GetXmin();
  double x_max = hist.GetXaxis()->GetXmax();
  double y_min = hist.GetYaxis()->GetXmin();
  double y_max = hist.GetYaxis()->GetXmax();
  color_map->data()->setSize(n_bins_x,n_bins_y);
  color_map->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double content = hist.GetBinContent(x,y);
      color_map->data()->setCell(x, y, content);
    }
  }
  QCPColorGradient gpHot = QCPColorGradient::gpHot;
  color_map->setGradient(gpHot.inverted());
  color_map->rescaleDataRange(true);
}

// function to write results to file (in columns)
void utilities::write_results_to_file(QString filename, QVector<QVector<double>> results)
{
  // open file
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	  std::cout<<"File could not be open for writing"<<std::endl; 
      return;

  QTextStream out(&file);
  QString s_time = "#time";
  QString s_total = "#total";
  QString s_elec = "#elec";
  QString s_holes = "#holes";

  int s_width = 8;

  out << s_time.leftJustified(s_width) << "\t" << s_total.leftJustified(s_width) << "\t" << s_elec.leftJustified(s_width) << "\t" << s_holes.leftJustified(s_width);

  for (int i = 0; i < results[0].size(); i++ )
  {
    QString line = QString("\n%1\t%2\t%3\t%4").arg(results[0][i], s_width, 'E', 4)
                                              .arg(results[1][i], s_width, 'E', 4)
                                              .arg(results[2][i], s_width, 'E', 4)
                                              .arg(results[3][i], s_width, 'E', 4);
    out << line;
  }
}

// function to write results to file (in rows)
void utilities::write_to_file_row(std::string filename, QVector<QVector<double>> results, double dt)
{
  // open file
//  QFile file(filename);
//  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
//	  std::cout<<"File could not be open for writing"<<std::endl; 
//      return;
//
//  QTextStream out(&file);
  unsigned long int steps = results[0].size();

  std::ofstream out(filename); // open file
  std::ostringstream temp_oss; // create string stream for precision output
  temp_oss << std::fixed ; 
  temp_oss << std::setprecision(15) << dt; // set det with enough precision
  
  // Scan on z
  for (int j = 0; j < results.size(); j++) 
  {
	  std::string line = temp_oss.str()+" "+std::to_string(steps)+" ";

	  // Scan on times
	  for (unsigned int i = 0; i < steps; i++ )
	  {
		line +=  std::to_string(results[j][i]) + " ";
			
	  }
		out << line << "\n";
  }
  
		temp_oss << std::setprecision(0) << dt; // set det with enough precision
        out.close();
}


// function to write results to file (in rows)
// overloaded (now from TH1D)
void utilities::write_to_file_row(std::string filename, TH1D *hconv, double temp, double height, double voltage)
{
  unsigned long int steps = hconv->GetNbinsX();
 height = height/1000.;

  std::ofstream out; // open file
	out.open(filename, std::ios_base::app);
	if (out.is_open())
	{
		out << steps << " ";
		out << temp-273. << " ";
		out << voltage << " ";
		out << "0 0 " << height << " ";


	  // Scan on times
	  for (unsigned int i = 1; i <= steps; i++ )
	  {
			out << std::fixed << std::setprecision(9) << hconv->GetBinContent(i) << " "; 
			//std::cout <<  steps << " " << i << " " << hconv->GetBinContent(i) << std::endl;
	  }

		//out << std::resetiosflags;
		out << std::endl;
    out.close();
	}
	else // Error output
	{
	std::cout << "File could not be read/created"<<std::endl;
	}
}


// function to write results to file (in rows)
// overloaded (now from TH1D)
void utilities::write_to_hetct_header(std::string filename, SMSDetector detector, double C, double dt, std::vector<double> z_shifts, double landa, std::string type, std::string carriers_file, std::vector<double> voltages)
{
	// Initialize stream for outputting to file
	std::ofstream header;  

	// Derived quantities
	int nX = 0;
	int nY = 0;
	int nZ = z_shifts.size();
	int nV = voltages.size();
	int nScans = (nX + nY + nZ)*nV;
	double deltaZ = 0;
	double deltaV = 0;	
	if (nZ > 1)	deltaZ = std::abs((z_shifts.back()-z_shifts.front())/(nZ-1));
	if (nV > 1) deltaV = std::abs((voltages.back()-voltages.back())/(nV-1)); 	
	double temp = detector.get_temperature() - 273;
	
	// Convert relevant quantities to string for outputting
	std::string s_date = "2015-06-15_19-25";
	std::string s_zVector = utilities::vector_to_string(z_shifts);
	std::string s_vVector = utilities::vector_to_string(voltages);
	
	// Open file
	header.open (filename);  

	// Check the file was open
	if (header.is_open())
	{ 
		// Store header string 
		header << "================\n";
		header << "SSD simulation file\n";
		header << "version: 1.0\n";
		header << "================\n";
		header << "scanType: " << "TRACS\n";
		header << "startTime: " << s_date << "\n" ;
		header << "landaLaser: " << landa << "\n";
		header << "illumDirect: " << type << "\n";
		header << "ampGain: 0\n";
		header << "fluence: " << detector.get_fluence() << "\n";
		header << "annealing: 0\n";
		header << "numOfScans: " << nScans << "\n";
		header << "temperatureMin: " << temp  << "\n";
		header << "temperatureMax: " << temp  << "\n";
		header << "deltaTemperature: " << 0. << "\n";
		header << "nTemperature: 1\n";
		header << "temperatureVector: " << temp  << "\n";
		header << "voltageMin: " << voltages.front() << "\n";
		header << "voltageMax: " << voltages.back() << "\n";
		header << "deltaVoltage: " << deltaV << "\n";
		header << "nVoltage: " << nV << "\n";
		header << "voltageVector: " << s_vVector << "\n";
		header << "pulsePowerIntensityMin: 60.000\n";
		header << "pulsePowerIntensityMax: 60.000\n";
		header << "deltaPulsePowerIntensity: 0\n";
		header << "nPulsePowerIntensity: 1\n";
		header << "pulsePowerIntensityVector: 60.000 \n";
		header << "pulseWidthMin: 0.000\n";
		header << "pulseWidthMax: 0.000\n";
		header << "deltaPulseWidth: 0\n";
		header << "nPulseWidth: 1\n";
		header << "pulseWidthVector: 0.000 \n";
		header << "xMin: 0.000\n";
		header << "xMax: 0.000\n";
		header << "deltaX: 0\n";
		header << "nX: 1\n";
		header << "yMin: 0.000\n";
		header << "yMax: 0.000\n";
		header << "deltaY: 0\n";
		header << "nY: 1\n";
		header << "yVector: 0.000 \n";
		header << "zMin: " << z_shifts.front() << "\n";
		header << "zMax: " << z_shifts.back() << "\n";
		header << "deltaZ: " << deltaZ << "\n";
		header << "nZ: " << nZ << "\n";
		header << "zVector: " << s_zVector << "\n";
		header << "At: " << std::fixed << std::setprecision(15) << dt << "\n";
		header << "Capacitance[F]: " <<  std::fixed << std::setprecision(15) << C << "\n";
		header << "Bulk: " << detector.get_bulk_type() << "\n";
		header << "Implant: " << detector.get_implant_type() << "\n";
		header << "NStrips: " << detector.get_nns() << "\n";
		header << "Pitch: " << std::setprecision(0) << detector.get_pitch() << "\n";
		header << "Width: " << std::setprecision(0) << detector.get_width() << "\n";
		header << "Depth: " << std::setprecision(0) << detector.get_depth() << "\n";
		header << "Vdep: " << detector.get_vdep() << "\n";
		header << "Carriers File: " << carriers_file << "\n";
		header << "================\n";
		header <<	 "Nt T[C] Vset[V] x[mm] y[mm] z[mm] I(t)[A]\n";
		header <<	 "================\n";

		header.close();
	}
	else // Error output
	{
	std::cout << "File could not be created"<<std::endl;
	}
}


// Utility to convert large arrays/vecto into strings
std::string utilities::vector_to_string(std::vector<double> input_list)
{
	// Use stringstream to convert doubles to string
	unsigned long int points = input_list.size();
	std::stringstream converter; 

	for (unsigned int i = 0; i < points; i++) {
		converter << input_list[i] << " ";
	}

	return converter.str();
}
