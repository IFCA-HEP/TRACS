#include "TRACSInterface.h"
//#include <iostream>
 #include <thread>
#include <mutex>          // std::mutex
#include "global.h"

//using namespace std;
//num_threads = 2;
//extern const int num_threads = 2;
std::mutex mtx;           // mutex for critical section
std::string fnm="Config.TRACS";
//num_threads = 8;//std::thread::hardware_concurrency();
int dummy_num_threads = 8;
void call_from_thread(int tid);
std::vector<TRACSInterface*> TRACSsim(num_threads);
//std::vector<TRACSInterface*> TRACSsim;
//std::vector<std::unique_ptr<TRACSInterface>> TRACSsim;
//num_threads = utilities::get_nthreads("Config.TRACS", num_threads);
//TRACSInterface *TRACSsim;
//TRACSInterface *TRACSsim[num_threads];// = new TRACSInterface( fnm ); //CORRECT
//TRACSInterface* TRACSsim = new TRACSInterface[num_threads];// = new TRACSInterface( fnm ); //CORRECT



      
/*
void resize_TI(size_t newSize)
{
    //TRACSInterface * TRACSsim = (TRACSInterface *) realloc(TRACSsim, newSize);
    //int* newArr = new int[newSize];
    TRACSInterface * TRACSsim_new[newSize];
    memcpy( TRACSsim_new, TRACSsim, num_threads * sizeof(TRACSInterface) );

    num_threads = newSize;
    delete [] TRACSsim;
    TRACSsim = TRACSsim_new;
}
 */     	
     

int main()
{
	std::cout << "Number of cores = " << std::thread::hardware_concurrency() <<std::endl;
	num_threads = std::thread::hardware_concurrency();
	//num_threads = 1;
	TRACSsim.resize(num_threads);
	//TRACSInterface *tp = NULL; // pointer
	//load up a vector
	//for (auto i = 0; i < num_threads; i++) {
    //TRACSsim.push_back(std::make_unique<TRACSInterface>(fnm));
	//}

	/*
	TRACSsim = new (malloc(sizeof(TRACSInterface))) TRACSInterface(fnm);
	TRACSInterface *MyDummy;
	num_threads = std::thread::hardware_concurrency();
	for(int i=1; i < num_threads; ++i)
	{
  		MyDummy = new(&TRACSsim[i]) TRACSInterface(fnm);
  	}
  	*/

	//num_threads = std::thread::hardware_concurrency();
	//Interface to TRACS
	//std::string fnm="Config.TRACS";
	//TRACSInterface TRACSsim(fnm);
	//TRACSInterface *TRACSsim = new TRACSInterface( fnm ); //CORRECT
/*	TRACSsim->set_zPos(30);
	TRACSsim->set_yPos(30);
	TRACSsim->set_vBias(400);
	TRACSsim->calculate_fields();
	TRACSsim->set_vBias(400);
	TRACSsim->set_trappingTime(std::numeric_limits<double>::max()-100);
	TRACSsim->set_neffType("Linear");
	TRACSsim->set_carrierFile("etct.carriers");
	TRACSsim->set_NeffParam({-25.,0.02,0.22,33.,0.,120.,220.,300});
	TRACSsim->calculate_fields();
	TRACSsim->simulate_ramo_current();
	TRACSsim->GetItRamo();
	TRACSsim->GetItRc();
	TRACSsim->GetItConv();
	*/
	//TRACSsim->loop_on("v");
	//TRACSsim->loop_on("y");
	//TRACSsim->loop_on("y","v","z");
	//TRACSsim->loop_on("v","y","z");
	//resize_TI(num_threads);
	std::thread t[num_threads];
 
         //Launch a group of threads
         for (int i = 0; i < num_threads; ++i) {
             t[i] = std::thread(call_from_thread, i);
         }
 
         std::cout << "Launched from the main\n";
 
         //Join the threads with the main thread
         for (int i = 0; i < num_threads; ++i) {
             t[i].join();
         }
    //write output to single file!
    TRACSsim[0]->write_to_file(0);

    
                // again
    // delete my class
	/*
	for(int i=0;i<dummy_num_threads;++i)
	{
   // Explicitly call the destructor for the placed object
  TRACSsim[i].~TRACSInterface();
	}
	free((void*)TRACSsim);
	TRACSsim = 0;
	*/
	//TRACSsim->loop_on("v","e");
	//TRACSsim->loop_on("v","y");
	//TRACSsim->loop_on("c");

	
	//std::cout<<"i_ramo:"<<TracsSim.GetItRamo()<<std::endl;
	//TracsSim.set_carrierFile("etct.carriers");
	//TRACSInterface* TracsSim = new TRACSInterface("Config.TRACS");
	
	return 0;
	
}


//This function will be called from a thread
  
      void call_from_thread(int tid) {
          //std::cout << "Launched by thread " << tid << std::endl;
      	 //std::string fnm="Config.TRACS";
      	std::cout << "Thread" << tid << std::endl;
      	mtx.lock();
      	TRACSsim[tid] = new TRACSInterface(fnm);
      	TRACSsim[tid]->set_tcount(tid);
      	if(tid==0)
      	{
      		i_ramo_array.clear();
			//i_conv_array.clear();
      		TRACSsim[tid]->resize_array();
      	}
       	TRACSsim[tid]->write_header(tid);
	    std::cout << "Made it t" << tid << std::endl;
	    mtx.unlock();
	    TRACSsim[tid]->loop_on(tid);
	    /*
	      	 if(tid==0){
	      	 		std::string fnm1="Config.TRACS";
	      	 	std::cout << "Thread1!" << std::endl;
	      	 	 mtx.lock();
	      	 	TRACSInterface *TRACSsim = new TRACSInterface(fnm1);
	      	 	std::cout << "Made it t1!" << std::endl;
	      	 	 mtx.unlock();
	         	TRACSsim->loop_on(1);
	      	 }
	      	 if(tid==1){
	      	 	std::string fnm2="Config.TRACS";
	      	 	std::cout << "Thread2!" << std::endl;
	      	 	//for(int i=0; i<100000000; i++);
	      	 	mtx.lock();
	      	 	TRACSInterface *TRACSsim2 = new TRACSInterface(fnm2);
	      	 	std::cout << "Made it t2!" << std::endl;
				mtx.unlock();
	      	 	TRACSsim2->loop_on(2);
	      	 }
	    */  	 
      	 }

