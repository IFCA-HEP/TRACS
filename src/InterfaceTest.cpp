#include "TRACSInterface.h"
//#include <iostream>
 
//using namespace std;


int main()
{
	//test
	std::cout << "It compiled, yaay!" << std::endl;
	//Interface to TRACS
	std::string fnm="Config.TRACS";
	//TRACSInterface TRACSsim(fnm);
	TRACSInterface *TRACSsim = new TRACSInterface( fnm );
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
	TRACSsim->loop_on("v","y","z");
	TRACSsim->loop_on("v","e");
	TRACSsim->loop_on("v","y");
	TRACSsim->loop_on("c");

	
	//std::cout<<"i_ramo:"<<TracsSim.GetItRamo()<<std::endl;
	//TracsSim.set_carrierFile("etct.carriers");
	//TRACSInterface* TracsSim = new TRACSInterface("Config.TRACS");


	//test
	std::cout << "A bit slow!" << std::endl;

	
	//while(true);
	
	return 0;
	
}
