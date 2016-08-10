#ifndef GLOBAL_H
#define GLOBAL_H
#include "TRACSInterface.h" 
#include <TH1D.h> // 1 Dimesional ROOT histogram 
#include <vector>
using std::vector;


extern vector<vector <TH1D*> >  i_ramo_array, i_conv_array, i_rc_array;
extern int num_threads;

#endif // GLOBAL_H
