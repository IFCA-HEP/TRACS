#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>

#include "qcustomplot.h"

#include <SMSDetector.h>
#include <Carrier.h>

#include <fstream>
#include <iterator>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void solve_fem();
    void show_weighting_potential_2d();
    void show_weighting_potential_3d();
    void show_electric_potential_2d();
    void show_electric_potential_3d();
    void drift_single_carrier();
    void drift_line_carrier();

    void show_carrier_map_qcp();
    void show_carrier_map_line();

private:
    Ui::MainWindow *ui;
    SMSDetector * detector;
    void init_weighting_potential_plot();
    void init_electric_potential_plot();
    void init_simple_currents_plot();
    void init_carrier_map_qcp();
    void init_weighting_field_map_qcp();
    void init_electric_field_map_qcp();


};

#endif // MAINWINDOW_H
