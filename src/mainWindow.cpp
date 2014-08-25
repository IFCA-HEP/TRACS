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

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::solve_fem()
{

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

