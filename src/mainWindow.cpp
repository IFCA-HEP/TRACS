#include "mainWindow.h"
#include "ui_mainWindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    detector(80., 30., 200., 3, 'p', 'n', 150, 150)
{
  ui->setupUi(this);



  connect(ui->solve_fem_button, SIGNAL(clicked()), this, SLOT(solve_fem()));
  connect(ui->show_weighting_pot_button, SIGNAL(clicked()), this, SLOT(show_weighting_potential()));
  connect(ui->show_electric_pot_button, SIGNAL(clicked()), this, SLOT(show_electric_potential()));
  connect(ui->single_carrier_button, SIGNAL(clicked()),this, SLOT(drift_single_carrier()));

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::solve_fem()
{
  double v_bias = 150.0;
  double v_depletion = 61.0;
  detector.set_voltages(v_bias, v_depletion);

  detector.solve_w_u();
  detector.solve_d_u();
  detector.solve_w_f_grad();
  detector.solve_d_f_grad();
}



void MainWindow::show_weighting_potential()
{
  // Get weighting potential
  Function * w_u = detector.get_w_u();

  ui->potentials_qcustomplot_map->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
  ui->potentials_qcustomplot_map->axisRect()->setupFullAxesBox(true);

  QCPColorMap *color_map_w_u = new QCPColorMap(ui->potentials_qcustomplot_map->xAxis, ui->potentials_qcustomplot_map->yAxis);
  ui->potentials_qcustomplot_map->addPlottable(color_map_w_u);

  int n_bins_x = 150;
  int n_bins_y = 150;

  double x_min = detector.get_x_min();
  double x_max = detector.get_x_max();
  double y_min = detector.get_y_min();
  double y_max = detector.get_y_max();

  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  color_map_w_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));

  for (int x=0; x<n_bins_x; ++x)
    for (int y=0; y<n_bins_y; ++y)
      color_map_w_u->data()->setCell(x, y, (*w_u)(x*step_x, y*step_y));

  // add a color scale:
  QCPColorScale *colorScale = new QCPColorScale(ui->potentials_qcustomplot_map);
  ui->potentials_qcustomplot_map->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
  colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
  color_map_w_u->setColorScale(colorScale); // associate the color map with the color scale
  colorScale->axis()->setLabel("Weighting Potential");

  color_map_w_u->setGradient(QCPColorGradient::gpPolar);
  color_map_w_u->rescaleDataRange(true);

  // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
  QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->potentials_qcustomplot_map);
  ui->potentials_qcustomplot_map->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
  colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);


  ui->potentials_qcustomplot_map->rescaleAxes();
  ui->potentials_qcustomplot_map->replot();

  // Plot weighting potential in an external window
  plot(reference_to_no_delete_pointer(*w_u),"Weighting Potential","auto");
  interactive();
}

void MainWindow::show_electric_potential()
{
  // Get electric potential
  Function * d_u = detector.get_d_u();
  // Plot electric potential in an external window
  plot(reference_to_no_delete_pointer(*d_u),"Electric Potential","auto");
  interactive();
}

void MainWindow::drift_single_carrier()
{


  // Create carrier and observe movement
  SMSDetector * dec_pointer = &detector;

  double dt = 1.e-11;
  double max_time = 10.e-9;
  // get number of steps from time
  const int max_steps = (int) std::floor(max_time / dt);

  std::valarray<double> curr_elec((size_t) max_steps);

  curr_elec = 0;

  Carrier electron('h', 1. , 300., 190. , dec_pointer, 0.0);

  curr_elec += electron.simulate_drift( dt , max_time, 300., 100.);
  QVector<double> x(max_steps), y(max_steps);

  std::ofstream output_file("./example.txt");

  for (int i=0; i< max_steps; i++)
  {
     x[i] = i*dt;
     y[i] = curr_elec[i];
     output_file << x[i] << ", ";
     output_file << curr_elec[i] << ", ";
     output_file << y[i] << ", ";
  }

  output_file.close();

  ui->carrier_curr_qcp->addGraph();
  ui->carrier_curr_qcp->graph(0)->setData(x,y);
  ui->carrier_curr_qcp->xAxis->setRange(0, max_time);
  ui->carrier_curr_qcp->yAxis->setRange(curr_elec.min(), curr_elec.max());
  ui->carrier_curr_qcp->replot();

}

