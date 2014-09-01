#include "mainWindow.h"
#include "ui_mainWindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    detector( new SMSDetector(80., 30., 200., 3, 'p', 'n', 150, 150))
{
  ui->setupUi(this);

  // setup potential tab initial plot settings
  init_weighting_potential_plot();
  init_electric_potential_plot();

  // setup currents tab initial plot settings
  init_simple_currents_plot();
  init_carrier_map_qcp();

  // setup fields tab initial plot settings
  init_weighting_field_map_qcp();
  init_electric_field_map_qcp();
  init_weighting_field_cut_qcp();
  init_electric_field_cut_qcp();

  // potentials tab conector
  connect(ui->solve_fem_button, SIGNAL(clicked()), this, SLOT(solve_fem()));
  connect(ui->show_weighting_pot_2d_button, SIGNAL(clicked()), this, SLOT(show_weighting_potential_2d()));
  connect(ui->show_weighting_pot_3d_button, SIGNAL(clicked()), this, SLOT(show_weighting_potential_3d()));
  connect(ui->show_electric_pot_2d_button, SIGNAL(clicked()), this, SLOT(show_electric_potential_2d()));
  connect(ui->show_electric_pot_3d_button, SIGNAL(clicked()), this, SLOT(show_electric_potential_3d()));

  // currents tab connectors
  connect(ui->s_carrier_button, SIGNAL(clicked()),this, SLOT(drift_single_carrier()));
  connect(ui->l_carrier_button, SIGNAL(clicked()),this, SLOT(drift_line_carrier()));
  connect(ui->view_carrier_line_button, SIGNAL(clicked()),this, SLOT(show_carrier_map_line()));
  connect(ui->carrier_map_qcp, SIGNAL(mouseDoubleClick(QMouseEvent*)), this, SLOT(carrier_from_click(QMouseEvent*)));

  // fields tab connectors
  connect(ui->show_w_field_map_mod_button, SIGNAL(clicked()), this, SLOT(show_w_field_mod_map()));
  connect(ui->show_w_field_map_x_button, SIGNAL(clicked()), this, SLOT(show_w_field_x_map()));
  connect(ui->show_w_field_map_y_button, SIGNAL(clicked()), this, SLOT(show_w_field_y_map()));
  connect(ui->show_e_field_map_mod_button, SIGNAL(clicked()), this, SLOT(show_e_field_mod_map()));
  connect(ui->show_e_field_map_x_button, SIGNAL(clicked()), this, SLOT(show_e_field_x_map()));
  connect(ui->show_e_field_map_y_button, SIGNAL(clicked()), this, SLOT(show_e_field_y_map()));
  connect(ui->show_w_field_3d_mod_button, SIGNAL(clicked()), this, SLOT(show_w_field_mod_3d()));
  connect(ui->show_w_field_3d_x_button, SIGNAL(clicked()), this, SLOT(show_w_field_x_3d()));
  connect(ui->show_w_field_3d_y_button, SIGNAL(clicked()), this, SLOT(show_w_field_y_3d()));
  connect(ui->show_e_field_3d_mod_button, SIGNAL(clicked()), this, SLOT(show_e_field_mod_3d()));
  connect(ui->show_e_field_3d_x_button, SIGNAL(clicked()), this, SLOT(show_e_field_x_3d()));
  connect(ui->show_e_field_3d_y_button, SIGNAL(clicked()), this, SLOT(show_e_field_y_3d()));
  connect(ui->w_field_vert_button, SIGNAL(clicked()), this, SLOT(show_w_field_vert_cut()));
  connect(ui->e_field_vert_button, SIGNAL(clicked()), this, SLOT(show_e_field_vert_cut()));
  connect(ui->w_field_hor_button, SIGNAL(clicked()), this, SLOT(show_w_field_hor_cut()));
  connect(ui->e_field_hor_button, SIGNAL(clicked()), this, SLOT(show_e_field_hor_cut()));

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::solve_fem()
{
  ui->fem_progress_bar->setValue(5);
  // get detector and solver variables
  double pitch = ui->pitch_double_box->value();
  double width = ui->width_double_box->value();
  double depth = ui->depth_double_box->value();
  int nns = ui->nn_strips_int_box->value();
  char bulk_type = ui->bulk_type_combo_box->currentText().toStdString().c_str()[0];
  char implant_type = ui->implant_type_combo_box->currentText().toStdString().c_str()[0];
  int n_cellsx = ui->n_cellsx_int_box->value();
  int n_cellsy = ui->n_cellsy_int_box->value();

  ui->fem_progress_bar->setValue(10);
  detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cellsx, n_cellsy);

  double v_bias = ui->bias_voltage_double_box->value();
  double v_depletion = ui->depletion_voltage_double_box->value();
  detector->set_voltages(v_bias, v_depletion);

  // solve potentials and field using FEM methods
  ui->fem_progress_bar->setValue(20);
  detector->solve_w_u();
  ui->fem_progress_bar->setValue(40);
  detector->solve_d_u();
  ui->fem_progress_bar->setValue(60);
  detector->solve_w_f_grad();
  ui->fem_progress_bar->setValue(80);
  detector->solve_d_f_grad();
  ui->fem_progress_bar->setValue(95);

  // plot solutions in 2d
  show_weighting_potential_2d();
  show_electric_potential_2d();
  show_carrier_map_qcp();
  show_w_field_mod_map();
  show_e_field_mod_map();
  ui->fem_progress_bar->setValue(100);

}



void MainWindow::show_weighting_potential_2d()
{
  // get weighting potential from detector instance
  Function * w_u = detector->get_w_u();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_w_u = qobject_cast<QCPColorMap *>(ui->weighting_pot_qcp->plottable(0));
  color_map_w_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
    for (int y=0; y<n_bins_y; ++y)
      color_map_w_u->data()->setCell(x, y, (*w_u)(x*step_x, y*step_y));
  color_map_w_u->setGradient(QCPColorGradient::gpPolar);
  color_map_w_u->rescaleDataRange(true);
  ui->weighting_pot_qcp->rescaleAxes();
  ui->weighting_pot_qcp->replot();
}

void MainWindow::show_weighting_potential_3d()
{
  // get weighting potential from detector instance
  Function * w_u = detector->get_w_u();
  // Plot weighting potential in an external window
  plot(reference_to_no_delete_pointer(*w_u),"Weighting Potential","auto");
  interactive();
}

void MainWindow::show_electric_potential_2d()
{
  // get electric potential from detector instance
  Function * d_u = detector->get_d_u();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_d_u = qobject_cast<QCPColorMap *>(ui->electric_pot_qcp->plottable(0));
  color_map_d_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_d_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
    for (int y=0; y<n_bins_y; ++y)
      color_map_d_u->data()->setCell(x, y, (*d_u)(x*step_x, y*step_y));
  color_map_d_u->setGradient(QCPColorGradient::gpPolar);
  color_map_d_u->rescaleDataRange(true);
  ui->electric_pot_qcp->rescaleAxes();
  ui->electric_pot_qcp->replot();
}

void MainWindow::show_electric_potential_3d()
{
  // Get electric potential
  Function * d_u = detector->get_d_u();
  // Plot electric potential in an external window
  plot(reference_to_no_delete_pointer(*d_u),"Electric Potential","auto");
  interactive();
}

void MainWindow::show_carrier_map_qcp()
{
  // get electric potential from detector instance
  Function * d_u = detector->get_d_u();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_d_u = qobject_cast<QCPColorMap *>(ui->carrier_map_qcp->plottable(0));
  color_map_d_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_d_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
    for (int y=0; y<n_bins_y; ++y)
      color_map_d_u->data()->setCell(x, y, (*d_u)(x*step_x, y*step_y));
  color_map_d_u->setGradient(QCPColorGradient::gpGrayscale);
  color_map_d_u->rescaleDataRange(true);
  ui->carrier_map_qcp->rescaleAxes();
  ui->carrier_map_qcp->replot();
}


void MainWindow::show_w_field_mod_map()
{
  // get weighting field from detector instance
  Function * w_f_grad = detector->get_w_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_w_f = qobject_cast<QCPColorMap *>(ui->weighting_field_map_qcp->plottable(0));
  color_map_w_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double w_f_x = ((*w_f_grad)[0])(x*step_x, y*step_y);
      double w_f_y = ((*w_f_grad)[1])(x*step_x, y*step_y);
      double w_f_mod = sqrt(w_f_x*w_f_x +  w_f_y *w_f_y);
      color_map_w_f->data()->setCell(x, y, w_f_mod);
    }
  }
  color_map_w_f->setGradient(QCPColorGradient::gpPolar);
  color_map_w_f->rescaleDataRange(true);
  ui->weighting_field_map_qcp->rescaleAxes();
  ui->weighting_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->w_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->w_field_hor_double->setValue(y_min + (y_max -y_min) / 2);

}

void MainWindow::show_e_field_mod_map()
{
  // get drift electric field from detector instance
  Function * e_f_grad = detector->get_d_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_e_f = qobject_cast<QCPColorMap *>(ui->electric_field_map_qcp->plottable(0));
  color_map_e_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_e_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double e_f_x = ((*e_f_grad)[0])(x*step_x, y*step_y);
      double e_f_y = ((*e_f_grad)[1])(x*step_x, y*step_y);
      double e_f_mod = sqrt(e_f_x*e_f_x +  e_f_y *e_f_y);
      color_map_e_f->data()->setCell(x, y, e_f_mod);
    }
  }
  color_map_e_f->setGradient(QCPColorGradient::gpPolar);
  color_map_e_f->rescaleDataRange(true);
  ui->electric_field_map_qcp->rescaleAxes();
  ui->electric_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->e_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->e_field_hor_double->setValue(y_min + (y_max -y_min) / 2);

}


void MainWindow::show_w_field_x_map()
{
  // get weighting field from detector instance
  Function * w_f_grad = detector->get_w_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_w_f = qobject_cast<QCPColorMap *>(ui->weighting_field_map_qcp->plottable(0));
  color_map_w_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double w_f_x = ((*w_f_grad)[0])(x*step_x, y*step_y);
      color_map_w_f->data()->setCell(x, y, w_f_x);
    }
  }
  color_map_w_f->setGradient(QCPColorGradient::gpPolar);
  color_map_w_f->rescaleDataRange(true);
  ui->weighting_field_map_qcp->rescaleAxes();
  ui->weighting_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->w_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->w_field_hor_double->setValue(y_min + (y_max -y_min) / 2);
}

void MainWindow::show_e_field_x_map()
{
  // get drift electric field from detector instance
  Function * e_f_grad = detector->get_d_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_e_f = qobject_cast<QCPColorMap *>(ui->electric_field_map_qcp->plottable(0));
  color_map_e_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_e_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double e_f_x = ((*e_f_grad)[0])(x*step_x, y*step_y);
      color_map_e_f->data()->setCell(x, y, e_f_x);
    }
  }
  color_map_e_f->setGradient(QCPColorGradient::gpPolar);
  color_map_e_f->rescaleDataRange(true);
  ui->electric_field_map_qcp->rescaleAxes();
  ui->electric_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->e_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->e_field_hor_double->setValue(y_min + (y_max -y_min) / 2);
}

void MainWindow::show_w_field_y_map()
{
  // get weighting field from detector instance
  Function * w_f_grad = detector->get_w_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_w_f = qobject_cast<QCPColorMap *>(ui->weighting_field_map_qcp->plottable(0));
  color_map_w_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double w_f_y = ((*w_f_grad)[1])(x*step_x, y*step_y);
      color_map_w_f->data()->setCell(x, y, w_f_y);
    }
  }
  color_map_w_f->setGradient(QCPColorGradient::gpPolar);
  color_map_w_f->rescaleDataRange(true);
  ui->weighting_field_map_qcp->rescaleAxes();
  ui->weighting_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->w_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->w_field_hor_double->setValue(y_min + (y_max -y_min) / 2);
}

void MainWindow::show_e_field_y_map()
{
  // get drift electric field from detector instance
  Function * e_f_grad = detector->get_d_f_grad();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_x = (x_max -x_min)/n_bins_x;
  double step_y = (y_max -y_min)/n_bins_y;

  // get plot and set new data
  QCPColorMap *color_map_e_f = qobject_cast<QCPColorMap *>(ui->electric_field_map_qcp->plottable(0));
  color_map_e_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_e_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  for (int x=0; x<n_bins_x; ++x)
  {
    for (int y=0; y<n_bins_y; ++y)
    {
      double e_f_y = ((*e_f_grad)[1])(x*step_x, y*step_y);
      color_map_e_f->data()->setCell(x, y, e_f_y);
    }
  }
  color_map_e_f->setGradient(QCPColorGradient::gpPolar);
  color_map_e_f->rescaleDataRange(true);
  ui->electric_field_map_qcp->rescaleAxes();
  ui->electric_field_map_qcp->replot();

  // set vertical and horizontal cuts to the middle
  ui->e_field_vert_double->setValue(x_min + (x_max -x_min) / 2);
  ui->e_field_hor_double->setValue(y_min + (y_max -y_min) / 2);
}

void MainWindow::show_w_field_mod_3d()
{
  // Get weighting field grad
  Function * w_f_grad = detector->get_w_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer(*w_f_grad),"Weighting Field Modulus","auto");
  interactive();
}

void MainWindow::show_w_field_x_3d()
{
 // Get weighting field grad
  Function * w_f_grad = detector->get_w_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer((*w_f_grad)[0]),"Weighting Field X Component","auto");
  interactive();
}

void MainWindow::show_w_field_y_3d()
{
  // Get weighting field grad
  Function * w_f_grad = detector->get_w_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer((*w_f_grad)[1]),"Weighting Field Y Component","auto");
  interactive();
}

void MainWindow::show_e_field_mod_3d()
{
  // Get weighting field grad
  Function * e_f_grad = detector->get_d_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer(*e_f_grad),"Electric Field Modulus","auto");
  interactive();
}

void MainWindow::show_e_field_x_3d()
{
 // Get weighting field grad
  Function * e_f_grad = detector->get_d_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer((*e_f_grad)[0]),"Electric Field X Component","auto");
  interactive();
}

void MainWindow::show_e_field_y_3d()
{
  // Get weighting field grad
  Function * e_f_grad = detector->get_d_f_grad();
  // Plot weighting field in an external window
  plot(reference_to_no_delete_pointer((*e_f_grad)[1]),"Electric Field Y Component","auto");
  interactive();
}

void MainWindow::show_w_field_vert_cut()
{
  // Get weighting field grad
  Function * w_f_grad = detector->get_w_f_grad();

  int index = ui->w_field_vert_combo->currentIndex();
  double x_cut_value = ui->w_field_vert_double->value();

  // get some required variables
  int n_bins_y = ui->n_cellsy_int_box->value();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_y = (y_max -y_min)/n_bins_y;

  ui->weighting_field_cut_qcp->xAxis->setLabel("Y (mum)");

  // create and fill vectors to plot
  QVector<double> y_position(n_bins_y), w_field(n_bins_y);
  for (int i = 0; i < n_bins_y; i++) {
    double y_value = y_min + i*step_y;
    y_position[i] = y_value;
    double w_f_x = ((*w_f_grad)[0])(x_cut_value, y_value);
    double w_f_y = ((*w_f_grad)[1])(x_cut_value, y_value);
    if (index == 0) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Weighting Field Modulus (m^-1)");
      w_field[i] = sqrt(w_f_x*w_f_x + w_f_y*w_f_y);
    } else if ( index == 1) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Weighting Field X (m^-1)");
      w_field[i] = w_f_x;
    } else if ( index == 2 ) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Weighting Field Y (m^-1)");
      w_field[i] = w_f_y;
    }
  }

  // delete previous graph and create new graph
  ui->weighting_field_cut_qcp->removeGraph(ui->weighting_field_cut_qcp->graph(0));
  QCPGraph * graph = ui->weighting_field_cut_qcp->addGraph();
  graph->setData(y_position, w_field);

  // reescale and plot
  graph->rescaleAxes();
  ui->weighting_field_cut_qcp->replot();
}

void MainWindow::show_w_field_hor_cut()
{

  // Get weighting field grad
  Function * w_f_grad = detector->get_w_f_grad();

  int index = ui->w_field_hor_combo->currentIndex();
  double y_cut_value = ui->w_field_hor_double->value();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double step_x = (x_max -x_min)/n_bins_x;

  ui->weighting_field_cut_qcp->xAxis->setLabel("X (mum)");

  // create and fill vectors to plot
  QVector<double> x_position(n_bins_x), w_field(n_bins_x);
  for (int i = 0; i < n_bins_x; i++) {
    double x_value = x_min + i*step_x;
    x_position[i] = x_value;
    double w_f_x = ((*w_f_grad)[0])(x_value, y_cut_value);
    double w_f_y = ((*w_f_grad)[1])(x_value, y_cut_value);
    if (index == 0) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Electric Field Modulus (m^-1)");
      w_field[i] = sqrt(w_f_x*w_f_x + w_f_y*w_f_y);
    } else if ( index == 1) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Weighting Field X (m^-1)");
      w_field[i] = w_f_x;
    } else if ( index == 2 ) {
      ui->weighting_field_cut_qcp->yAxis->setLabel("Weighting Field Y (m^-1)");
      w_field[i] = w_f_y;
    }
  }

  // delete previous graph and create new graph
  ui->weighting_field_cut_qcp->removeGraph(ui->weighting_field_cut_qcp->graph(0));
  QCPGraph * graph = ui->weighting_field_cut_qcp->addGraph();
  graph->setData(x_position, w_field);

  // reescale and plot
  graph->rescaleAxes();
  ui->weighting_field_cut_qcp->replot();

}

void MainWindow::show_e_field_vert_cut()
{
  // Get electric field grad
  Function * e_f_grad = detector->get_d_f_grad();

  int index = ui->e_field_vert_combo->currentIndex();
  double x_cut_value = ui->e_field_vert_double->value();

  // get some required variables
  int n_bins_y = ui->n_cellsy_int_box->value();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();
  double step_y = (y_max -y_min)/n_bins_y;

  ui->electric_field_cut_qcp->xAxis->setLabel("Y (mum)");

  // create and fill vectors to plot
  QVector<double> y_position(n_bins_y), e_field(n_bins_y);
  for (int i = 0; i < n_bins_y; i++) {
    double y_value = y_min + i*step_y;
    y_position[i] = y_value;
    double e_f_x = ((*e_f_grad)[0])(x_cut_value, y_value);
    double e_f_y = ((*e_f_grad)[1])(x_cut_value, y_value);
    if (index == 0) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Electric Field Modulus (m^-1)");
      e_field[i] = sqrt(e_f_x*e_f_x + e_f_y*e_f_y);
    } else if ( index == 1) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Weighting Field X (m^-1)");
      e_field[i] = e_f_x;
    } else if ( index == 2 ) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Weighting Field Y (m^-1)");
      e_field[i] = e_f_y;
    }
  }

  // delete previous graph and create new graph
  ui->electric_field_cut_qcp->removeGraph(ui->electric_field_cut_qcp->graph(0));
  QCPGraph * graph = ui->electric_field_cut_qcp->addGraph();
  graph->setData(y_position, e_field);

  // reescale and plot
  graph->rescaleAxes();
  ui->electric_field_cut_qcp->replot();

}

void MainWindow::show_e_field_hor_cut()
{

  // Get weighting field grad
  Function * e_f_grad = detector->get_d_f_grad();

  int index = ui->e_field_hor_combo->currentIndex();
  double y_cut_value = ui->e_field_hor_double->value();

  // get some required variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double step_x = (x_max -x_min)/n_bins_x;

  ui->electric_field_cut_qcp->xAxis->setLabel("Y (mum)");

  // create and fill vectors to plot
  QVector<double> x_position(n_bins_x), e_field(n_bins_x);
  for (int i = 0; i < n_bins_x; i++) {
    double x_value = x_min + i*step_x;
    x_position[i] = x_value;
    double e_f_x = ((*e_f_grad)[0])(x_value, y_cut_value);
    double e_f_y = ((*e_f_grad)[1])(x_value, y_cut_value);
    if (index == 0) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Electric Field Modulus (m^-1)");
      e_field[i] = sqrt(e_f_x*e_f_x + e_f_y*e_f_y);
    } else if ( index == 1) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Weighting Field X (m^-1)");
      e_field[i] = e_f_x;
    } else if ( index == 2 ) {
      ui->electric_field_cut_qcp->yAxis->setLabel("Weighting Field Y (m^-1)");
      e_field[i] = e_f_y;
    }
  }

  // delete previous graph and create new graph
  ui->electric_field_cut_qcp->removeGraph(ui->electric_field_cut_qcp->graph(0));
  QCPGraph * graph = ui->electric_field_cut_qcp->addGraph();
  graph->setData(x_position, e_field);

  // reescale and plot
  graph->rescaleAxes();
  ui->electric_field_cut_qcp->replot();

}


void MainWindow::drift_single_carrier()
{

  double dt = 1.e-12 * ui->s_carrier_time_step_box->value();
  double max_time = 1.e-9 * ui->s_carrier_max_time_box->value();
  // get number of steps from time
  int max_steps = (int) std::floor(max_time / dt);
  // get carrier charge in coulombs from input
  double carrier_q = 1.6021765710e-19 * ui->s_carrier_q_box->value();

  std::valarray<double> curr_elec((size_t) max_steps);
  std::valarray<double> curr_hole((size_t) max_steps);
  std::valarray<double> curr_total((size_t) max_steps);

  curr_elec = 0;
  curr_hole = 0;

  double x_pos = ui->s_carrier_x_pos_box->value();
  double y_pos = ui->s_carrier_y_pos_box->value();

  Carrier electron('e', -1.*carrier_q , x_pos, y_pos , detector, 0.0);
  Carrier hole('h', carrier_q, x_pos, y_pos, detector, 0.0);

  int index = ui->s_carrier_type->currentIndex();

  if ( index == 0 || index == 1)
  {
  curr_elec += electron.simulate_drift( dt , max_time, x_pos, y_pos);
  }
  if ( index == 0 || index == 2)
  {
  curr_hole += hole.simulate_drift( dt , max_time, x_pos, y_pos);
  }
  curr_total = curr_elec + curr_hole;

  QVector<double> x_elec(max_steps), y_elec(max_steps);
  QVector<double> x_hole(max_steps), y_hole(max_steps);
  QVector<double> x_total(max_steps), y_total(max_steps);

  for (int i=0; i< max_steps; i++)
  {
     x_elec[i] = i*dt;
     x_hole[i] = i*dt;
     x_total[i] = i*dt;
     y_elec[i] = curr_elec[i];
     y_hole[i] = curr_hole[i];
     y_total[i] = curr_total[i];
  }

  // set new data
  ui->carrier_curr_qcp->graph(0)->setData(x_elec,y_elec);
  ui->carrier_curr_qcp->graph(1)->setData(x_hole,y_hole);
  ui->carrier_curr_qcp->graph(2)->setData(x_total, y_total);

  // reescale and plot
  if ( index == 0 || index == 1)
  {
  ui->carrier_curr_qcp->graph(0)->rescaleAxes();
  ui->carrier_curr_qcp->graph(1)->rescaleAxes(true);
  }
  if ( index == 0 || index == 2)
  {
  ui->carrier_curr_qcp->graph(1)->rescaleAxes();
  ui->carrier_curr_qcp->graph(0)->rescaleAxes(true);
  }
  ui->carrier_curr_qcp->graph(2)->rescaleAxes(true);
  ui->carrier_curr_qcp->replot();

}

void MainWindow::drift_line_carrier()
{
  double dt = 1.e-12 * ui->s_carrier_time_step_box->value();
  double max_time = 1.e-9 * ui->s_carrier_max_time_box->value();
  // get number of steps from time
  int max_steps = (int) std::floor(max_time / dt);
  // get carrier charge in coulombs from input
  double carrier_q = 1.6021765710e-19 * ui->s_carrier_q_box->value();

  std::valarray<double> curr_elec((size_t) max_steps);
  std::valarray<double> curr_hole((size_t) max_steps);
  std::valarray<double> curr_total((size_t) max_steps);

  curr_elec = 0;
  curr_hole = 0;
  curr_total = 0;

  // get data from ui
  double l_start_x = ui->l_carrier_startx_box->value();
  double l_start_y = ui->l_carrier_starty_box->value();
  double l_end_x = ui->l_carrier_endx_box->value();
  double l_end_y = ui->l_carrier_endy_box->value();
  double l_carrier_sep = ui->l_carrier_sep_box->value();

  // obtain displacement, length and normalized vector components
  double l_delta_x = l_end_x - l_start_x;
  double l_delta_y = l_end_y - l_start_y;
  double l_length = sqrt(l_delta_x*l_delta_x + l_delta_y*l_delta_y);
  double l_vecnorm_x = l_delta_x / l_length;
  double l_vecnorm_y = l_delta_y / l_length;

  int n_carriers = (int) std::floor(l_length / l_carrier_sep);

  // check carrier type
  int index = ui->s_carrier_type->currentIndex();

  // init default carrier position
  double x_pos = 0;
  double y_pos = 0;
  // create instance of carrier types
  Carrier electron('e', -1.*carrier_q , x_pos, y_pos , detector, 0.0);
  Carrier hole('h', carrier_q, x_pos, y_pos, detector, 0.0);

  for (int c=0; c< n_carriers; c++)
  {
    x_pos = l_start_x + l_vecnorm_x*c*l_carrier_sep;
    y_pos = l_start_y + l_vecnorm_y*c*l_carrier_sep;

    if ( index == 0 || index == 1)
    {
    curr_elec += electron.simulate_drift( dt , max_time, x_pos, y_pos);
    }
    if ( index == 0 || index == 2)
    {
    curr_hole += hole.simulate_drift( dt , max_time, x_pos, y_pos);
    }
  }

  curr_total = curr_elec + curr_hole;

  QVector<double> x_elec(max_steps), y_elec(max_steps);
  QVector<double> x_hole(max_steps), y_hole(max_steps);
  QVector<double> x_total(max_steps), y_total(max_steps);

  for (int i=0; i< max_steps; i++)
  {
     x_elec[i] = i*dt;
     x_hole[i] = i*dt;
     x_total[i] = i*dt;
     y_elec[i] = curr_elec[i];
     y_hole[i] = curr_hole[i];
     y_total[i] = curr_total[i];
  }

  // set new data
  ui->carrier_curr_qcp->graph(0)->setData(x_elec, y_elec);
  ui->carrier_curr_qcp->graph(1)->setData(x_hole, y_hole);
  ui->carrier_curr_qcp->graph(2)->setData(x_total, y_total);

  // reescale and plot
  if ( index == 0 || index == 1)
  {
  ui->carrier_curr_qcp->graph(0)->rescaleAxes();
  ui->carrier_curr_qcp->graph(1)->rescaleAxes(true);
  }
  if ( index == 0 || index == 2)
  {
  ui->carrier_curr_qcp->graph(1)->rescaleAxes();
  ui->carrier_curr_qcp->graph(0)->rescaleAxes(true);
  }
  ui->carrier_curr_qcp->graph(2)->rescaleAxes(true);
  ui->carrier_curr_qcp->replot();

}

void MainWindow::show_carrier_map_line()
{

  // get last line and exist if it exists
  ui->carrier_map_qcp->removeItem(ui->carrier_map_qcp->item(0));

  // create a line item object
  QCPItemLine * carrier_line = new QCPItemLine(ui->carrier_map_qcp);
  ui->carrier_map_qcp->addItem(carrier_line);

  // get data from ui
  double l_start_x = ui->l_carrier_startx_box->value();
  double l_start_y = ui->l_carrier_starty_box->value();
  double l_end_x = ui->l_carrier_endx_box->value();
  double l_end_y = ui->l_carrier_endy_box->value();

  // set coordinates, color and replot
  carrier_line->start->setType( QCPItemPosition::ptPlotCoords);
  carrier_line->start->setCoords(l_start_x, l_start_y);
  carrier_line->end->setType( QCPItemPosition::ptPlotCoords);
  carrier_line->end->setCoords(l_end_x, l_end_y);
  carrier_line->setPen(QPen(Qt::red));
  ui->carrier_map_qcp->replot();

}

void MainWindow::carrier_from_click(QMouseEvent * event)
{

  double x_pos = ui->carrier_map_qcp->xAxis->pixelToCoord(event->pos().x());
  double y_pos = ui->carrier_map_qcp->yAxis->pixelToCoord(event->pos().y());

  ui->s_carrier_x_pos_box->setValue(x_pos);
  ui->s_carrier_y_pos_box->setValue(y_pos);

  drift_single_carrier();

}

void MainWindow::init_weighting_potential_plot()
{
  // allow resizing and dragging
  ui->weighting_pot_qcp->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
  ui->weighting_pot_qcp->axisRect()->setupFullAxesBox(true);

  // create and add color map object
  QCPColorMap *color_map_w_u = new QCPColorMap(ui->weighting_pot_qcp->xAxis, ui->weighting_pot_qcp->yAxis);
  ui->weighting_pot_qcp->addPlottable(color_map_w_u);

  // add a color scale and set to the right of the main axis rect
  QCPColorScale *colorScale = new QCPColorScale(ui->weighting_pot_qcp);
  ui->weighting_pot_qcp->plotLayout()->addElement(0, 1, colorScale);
  colorScale->setType(QCPAxis::atRight);
  color_map_w_u->setColorScale(colorScale);
  colorScale->axis()->setLabel("Weighting Potential");

  // get inital variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();

  // set ranges and color gradients
  color_map_w_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  color_map_w_u->setGradient(QCPColorGradient::gpPolar);
  color_map_w_u->rescaleDataRange(true);

  // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
  QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->weighting_pot_qcp);
  ui->weighting_pot_qcp->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
  colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

  // reescale and replot
  ui->weighting_pot_qcp->rescaleAxes();
  ui->weighting_pot_qcp->replot();
}

void MainWindow::init_electric_potential_plot()
{
  // allow resizing and dragging
  ui->electric_pot_qcp->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
  ui->electric_pot_qcp->axisRect()->setupFullAxesBox(true);

  // create and add color map object
  QCPColorMap *color_map_d_u = new QCPColorMap(ui->electric_pot_qcp->xAxis, ui->electric_pot_qcp->yAxis);
  ui->electric_pot_qcp->addPlottable(color_map_d_u);

  // add a color scale and set to the right of the main axis rect
  QCPColorScale *colorScale = new QCPColorScale(ui->electric_pot_qcp);
  ui->electric_pot_qcp->plotLayout()->addElement(0, 1, colorScale);
  colorScale->setType(QCPAxis::atRight);
  color_map_d_u->setColorScale(colorScale);
  colorScale->axis()->setLabel("Electric Potential");

  // get inital variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();

  // set ranges and color gradients
  color_map_d_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_d_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  color_map_d_u->setGradient(QCPColorGradient::gpPolar);
  color_map_d_u->rescaleDataRange(true);

  // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
  QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->electric_pot_qcp);
  ui->electric_pot_qcp->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
  colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

  // reescale and replot
  ui->electric_pot_qcp->rescaleAxes();
  ui->electric_pot_qcp->replot();
}

void MainWindow::init_simple_currents_plot()
{
  ui->carrier_curr_qcp->legend->setVisible(true);
  ui->carrier_curr_qcp->yAxis->setLabel("Current (A)");
  ui->carrier_curr_qcp->xAxis->setLabel("Time (s)");
  ui->carrier_curr_qcp->addGraph();
  ui->carrier_curr_qcp->graph(0)->setPen(QPen(Qt::blue));
  ui->carrier_curr_qcp->graph(0)->setName("electron");
  ui->carrier_curr_qcp->addGraph();
  ui->carrier_curr_qcp->graph(1)->setPen(QPen(Qt::red));
  ui->carrier_curr_qcp->graph(1)->setName("hole");
  ui->carrier_curr_qcp->addGraph();
  ui->carrier_curr_qcp->graph(2)->setPen(QPen(Qt::black));
  ui->carrier_curr_qcp->graph(2)->setName("total");
  ui->carrier_curr_qcp->graph(0)->rescaleAxes();
  ui->carrier_curr_qcp->graph(1)->rescaleAxes(true);
  ui->carrier_curr_qcp->graph(2)->rescaleAxes(true);
  ui->carrier_curr_qcp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
  ui->carrier_curr_qcp->replot();
}

void MainWindow::init_carrier_map_qcp()
{
  // not allow resizing and dragging
  // ui->carrier_map_qcp->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
  ui->carrier_map_qcp->axisRect()->setupFullAxesBox(true);

  // create and add color map object
  QCPColorMap *color_map_d_u = new QCPColorMap(ui->carrier_map_qcp->xAxis, ui->carrier_map_qcp->yAxis);
  ui->carrier_map_qcp->addPlottable(color_map_d_u);

  // get inital variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();

  // set ranges and color gradients
  color_map_d_u->data()->setSize(n_bins_x,n_bins_y);
  color_map_d_u->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  color_map_d_u->setGradient(QCPColorGradient::gpGrayscale);
  color_map_d_u->rescaleDataRange(true);

  // reescale and replot
  ui->carrier_map_qcp->rescaleAxes();
  ui->carrier_map_qcp->replot();
}

void MainWindow::init_weighting_field_map_qcp()
{
  // allow resizing and dragging
  ui->weighting_field_map_qcp->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
  ui->weighting_field_map_qcp->axisRect()->setupFullAxesBox(true);

  // create and add color map object
  QCPColorMap *color_map_w_f = new QCPColorMap(ui->weighting_field_map_qcp->xAxis, ui->weighting_field_map_qcp->yAxis);
  ui->weighting_field_map_qcp->addPlottable(color_map_w_f);

  // add a color scale and set to the right of the main axis rect
  QCPColorScale *colorScale = new QCPColorScale(ui->weighting_field_map_qcp);
  ui->weighting_field_map_qcp->plotLayout()->addElement(0, 1, colorScale);
  colorScale->setType(QCPAxis::atRight);
  color_map_w_f->setColorScale(colorScale);
  colorScale->axis()->setLabel("Weighting Field");

  // get inital variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();

  // set ranges and color gradients
  color_map_w_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_w_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  color_map_w_f->setGradient(QCPColorGradient::gpPolar);
  color_map_w_f->rescaleDataRange(true);

  // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
  QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->weighting_field_map_qcp);
  ui->weighting_field_map_qcp->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
  colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

  // reescale and replot
  ui->weighting_field_map_qcp->rescaleAxes();
  ui->weighting_field_map_qcp->replot();
}


void MainWindow::init_electric_field_map_qcp()
{
  // allow resizing and dragging
  ui->electric_field_map_qcp->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
  ui->electric_field_map_qcp->axisRect()->setupFullAxesBox(true);

  // create and add color map object
  QCPColorMap *color_map_e_f = new QCPColorMap(ui->electric_field_map_qcp->xAxis, ui->electric_field_map_qcp->yAxis);
  ui->electric_field_map_qcp->addPlottable(color_map_e_f);

  // add a color scale and set to the right of the main axis rect
  QCPColorScale *colorScale = new QCPColorScale(ui->electric_field_map_qcp);
  ui->electric_field_map_qcp->plotLayout()->addElement(0, 1, colorScale);
  colorScale->setType(QCPAxis::atRight);
  color_map_e_f->setColorScale(colorScale);
  colorScale->axis()->setLabel("Electric Field");

  // get inital variables
  int n_bins_x = ui->n_cellsx_int_box->value();
  int n_bins_y = ui->n_cellsy_int_box->value();
  double x_min = detector->get_x_min();
  double x_max = detector->get_x_max();
  double y_min = detector->get_y_min();
  double y_max = detector->get_y_max();

  // set ranges and color gradients
  color_map_e_f->data()->setSize(n_bins_x,n_bins_y);
  color_map_e_f->data()->setRange(QCPRange(x_min, x_max), QCPRange(y_min, y_max));
  color_map_e_f->setGradient(QCPColorGradient::gpPolar);
  color_map_e_f->rescaleDataRange(true);

  // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
  QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->electric_field_map_qcp);
  ui->electric_field_map_qcp->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
  colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

  // reescale and replot
  ui->electric_field_map_qcp->rescaleAxes();
  ui->electric_field_map_qcp->replot();
}

void MainWindow::init_weighting_field_cut_qcp()
{
  ui->weighting_field_cut_qcp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
  ui->weighting_field_cut_qcp->axisRect()->setupFullAxesBox(true);
  ui->weighting_field_cut_qcp->replot();
}

void MainWindow::init_electric_field_cut_qcp()
{
  ui->electric_field_cut_qcp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
  ui->electric_field_cut_qcp->axisRect()->setupFullAxesBox(true);
  ui->electric_field_cut_qcp->replot();
}


