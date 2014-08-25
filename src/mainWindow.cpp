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
}



void MainWindow::show_weighting_potential()
{
  // Get weighting potential
  Function * w_u = detector.get_w_u();
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

