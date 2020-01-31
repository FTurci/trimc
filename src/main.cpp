// Francesco Turci f.turci@bristol.ac.uk 2019
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "cubic_box.h"
#include "particle.h"
#include "system.h"
#include "histogram.h"
#include <cmath>
#include <string>
#include <typeinfo>
#include <map>
#include "CLI11.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){
	CLI::App app("Monte-Carlo simulations woth two and threebody interactions.");
	// default values
    std::string filename = "default.xyz";
    double size = 8;
    int npart = 500;
	int totalsweeps = 100000;

	bool cell = 1;
	bool threebody = 0;
	std::string Pe = "10.0";

	std::string type = "fcc";
	

	double Temperature =1.0;
	double beta = 1.0/Temperature;
	double epsilon =1.0;
	double sigma = 1.0;
	double rcut =2.5;
	bool shift = true;

	double rho = -1;

    app.add_option("-f,--file", filename, "path to  xyz file");

    app.add_option("-s,--size", size, "box size");
    app.add_option("-p,--Pe", Pe, "Peclet");
    app.add_option("-r,--rho", rho, "number density");
    app.add_option("-t,--threebody", threebody, "threebody interaction");
    CLI11_PARSE(app, argc, argv);
	Histogram histo(0,2*size,600);


	System <MultiBodyTabularInteraction> model (size, npart);
	model.cell = cell;

	double rhi= 2.5;
	model.interaction.twobody.assign("../Interactions/depletion-Pe"+Pe+"L4.0.txt",118,9.399999999999999467e-01, 2.5);
	// model.interaction.twobody.assign("depletion-Pe10.0L4.0.txt",118,9.399999999999999467e-01, 2.5);

	std::vector<int> sizes {100,100,100};
	std::vector<double> los {0.0,0.0,0.0} ;
	std::vector<double> his {5.0,5.0,5.0} ;

	model.interaction.threebody.assign("../Interactions/threebodypotential-minustwo-L5-Pe"+Pe+".txt",sizes,los,his );
	// model.interaction.threebody.assign("threebodypotential.txt",sizes,los,his );
	double range = rhi;
	

	if (rho>0)	model.crystal(type,rho);

	// model.box.sides[2]*=2.0;
	// model.read_last_from_xyz("rho1.2-Pe10.0.xyz");

	model.build_cells(range);

	system(("rm "+filename).c_str());
	model.dump(filename);
	printf("Box size %g\n", model.box.sides[0]);
	
	// warmup without threebody
	int warmup =10000;

	for (int sweep = 0;sweep < warmup;sweep++){
			std::cout<<sweep<<std::endl;
			model.sweep(beta, 0.2,0);
		}
	model.dump(filename);
	for (int sweep = 0;sweep < totalsweeps;sweep++){
			std::cout<<sweep<<std::endl;
			model.sweep(beta, 0.2,threebody);
		// model.local_MC_step(beta, threebody);

		if(sweep%100==0)
			{
			// double E = model.compute_energy()/model.num_particles;
			// std::cerr<<sweep <<"\t"<<E<<std::endl;
			model.dump(filename);
		}

	}
	printf("cell %s\n", model.cell?"true":"false");


	
}
