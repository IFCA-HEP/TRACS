#include "TRACSInterface.h"
//#include <iostream>
 
//using namespace std;


int main()
{
	//test
	std::cout << "It compiled, yaay!" << std::endl;
	//Interface to TRACS
	std::string fnm="Config.TRACS";
	TRACSInterface TRACSsim(fnm);
	TRACSsim.simulate_ramo_current();
	TRACSsim.GetItRamo();


	//TRACSInterface *TRACSsim = new TRACSInterface( fnm );
	//std::cout<<"i_ramo:"<<TracsSim.GetItRamo()<<std::endl;
	//TracsSim.set_carrierFile("etct.carriers");
	//TRACSInterface* TracsSim = new TRACSInterface("Config.TRACS");


	//test
	std::cout << "A bit slow!" << std::endl;

	
	//while(true);
	
	return 0;
	
}
