/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

#include "TRACSInterface.h"
#include <thread>
#include <mutex>
#include "global.h"

//using namespace std;
//num_threads = 2;
//extern const int num_threads = 2;
std::mutex mtx;           // mutex for critical sections
std::string fnm="Config.TRACS";
std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
//int init_num_threads; // initial number of threads, which might change dynamically

void call_from_thread(int tid);

//int main(int nthreads = std::thread::hardware_concurrency())
int main(int argc, char* argv[])
{
	std::cout << "Maximum number of available cores = " << std::thread::hardware_concurrency() <<std::endl;
	if(argc<2)
	{
		num_threads = std::thread::hardware_concurrency(); // No. of threads = No. of cores
	}
	else
		num_threads = atoi(argv[1]);
	if(num_threads == 0){
		num_threads = 1;
	}
	std::cout << "The execution will use " << num_threads << " Thread(s) "  <<std::endl;

	TRACSsim.resize(num_threads);
	t.resize(num_threads);

	//Launching threads
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i);
	}

	//Joining threads
	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	//write output to single file!
	TRACSsim[0]->write_to_file(0);
	//Finalizing the execution
	//getter test
	std::vector<double> neff_test = TRACSsim[0]->get_NeffParam();
	std::cout << "Neff param.: " << std::endl;
	for (int i = 0; i < 8; i++)
	{
		std::cout << neff_test[i] << std::endl;
	}

	for (unsigned int i = 0; i < TRACSsim.size(); i++)
	{
	    delete TRACSsim[i];
	}
	std::quick_exit(1);
	return 0;
}

//This function will be called from a thread
void call_from_thread(int tid) {
	// every thread instantiates a new TRACSInterface object

	mtx.lock();
	std::cout << "Thread with tid " << tid << " is INSIDE the critical section "<< std::endl;
	TRACSsim[tid] = new TRACSInterface(fnm);
	TRACSsim[tid]->set_tcount(tid);
	if(tid==0)
	{
		i_ramo_array.clear();
		TRACSsim[tid]->resize_array();
		TRACSsim[tid]->write_header(tid);
		TRACSsim.resize(num_threads);
	}
	std::cout << "Thread with tid " << tid << " is OUTSIDE the critical section "<< std::endl;
	mtx.unlock();
	std::cout << "Thread with tid " << tid << " simulating ramo current - drifting "<< std::endl;
	TRACSsim[tid]->loop_on(tid);



}

