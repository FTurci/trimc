// Francesco Turci f.turci@bristol.ac.uk 2019
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "../cubic_box.h"
#include "../particle.h"
#include "../system.h"
#include "../histogram.h"
#include <cmath>
#include <string>
#include <typeinfo>
#include <map>
using namespace std;
using namespace Eigen;


int test_threebody(){

	System <MultiBodyTabularInteraction> model (8.0, 8);
	std::string type = "fcc";
	model.crystal(type,1.2);

	std::vector<int> sizes {100,100,100};
	std::vector<double> los {0.0,0.0,0.0} ;
	std::vector<double> his {5.0,5.0,5.0} ;

	model.interaction.threebody.assign("../../Interactions/threebodypotential-L5-Pe60.0.txt",sizes,los,his );

	int errors = 0;
	double val1 =model.interaction.at(1.0,0.9,1.1);
	double val2 =model.interaction.at(1.2, 1.1,1.05);
	printf("val2 %g\n",val2);
	double val3 =model.interaction.at(1.05,1.1,1.0);

	if (val1>0) {
		std::cout<<"!!! NONZERO r1>r2 and r2<r3 : "<<val1<<std::endl;
		errors +=1;
	}
	if (val2>0) {
		std::cout<<"!!! NONZERO r1>r2>r3 : "<<val2<<std::endl;
		errors +=1;
	}
	if (val3>0) {
		std::cout<<"!!! NONZERO r1<r2 and r2>r3 : "<<val3<<std::endl;
		errors +=1;
	}

	return errors;
}

int test_compute(){
		System <MultiBodyTabularInteraction> model (8.0, 8);
	// std::string type = "fcc";
	// model.crystal(type,1.2);
	model.compute();
	int errors=0;
	return errors;
}
int main(int argc, char *argv[]){
	if(test_threebody()>0) std::cout<<"[TEST NOT PASSED]: test_threebody"<<std::endl;
	else std::cout<<"[TEST PASSED]: test_threebody"<<std::endl;
	
	test_compute();
}
