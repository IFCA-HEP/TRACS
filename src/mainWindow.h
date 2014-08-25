#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>

#include <SMSDetector.h>
#include <Carrier.h>

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
    void show_weighting_potential();
    void show_electric_potential();

private:
    Ui::MainWindow *ui;
    SMSDetector detector;


};

#endif // MAINWINDOW_H
