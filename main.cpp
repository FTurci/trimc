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

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){
	// container of strings
	std::vector<std::string> all_args(argv, argv + argc);

	double size = atof(argv[1]);
	int npart = atoi(argv[2]);
	int totalsweeps = atoi(argv[3]);
	bool cell = bool(atoi(argv[4]));
	bool threebody = bool(atoi(argv[5]));
	std::string filename =all_args[6];

	std::string Pe = all_args[7];

	unsigned long long nsteps = 50000000;
	Histogram histo(0,2*size,600);

	double Temperature =1.0;
	double beta = 1.0/Temperature;
	double epsilon =1.0;
	double sigma = 1.0;
	double rcut =2.5;
	bool shift = true;
	System <MultiBodyTabularInteraction> model (size, npart);
	model.cell = cell;

	double rhi= 2.5;
	model.interaction.twobody.assign("Interactions/depletion-Pe"+Pe+"L4.0.txt",118,9.399999999999999467e-01, 2.5);
	// model.interaction.twobody.assign("depletion-Pe10.0L4.0.txt",118,9.399999999999999467e-01, 2.5);

	std::vector<int> sizes {100,100,100};
	std::vector<double> los {0.0,0.0,0.0} ;
	std::vector<double> his {4.0,4.0,4.0} ;

	model.interaction.threebody.assign("Interactions/threebodypotential-L5-Pe"+Pe+".txt",sizes,los,his );
	// model.interaction.threebody.assign("threebodypotential.txt",sizes,los,his );
	double range = rhi;
	
	std::string type = "fcc";
	model.crystal(type,0.6);
	// model.read_last_from_xyz("small.xyz");
	
	model.box.sides[2]*=2.0;

	model.build_cells(range);

	system(("rm "+filename).c_str());
	model.dump(filename);
	printf("New box size %g\n", model.box.sides[0]);
	
	for (int sweep = 0;sweep < totalsweeps;sweep++){
			std::cout<<sweep<<std::endl;
			model.sweep(beta, 0.2,threebody);
		// model.local_MC_step(beta, threebody);

		if(sweep%10==0)
			{
			// double E = model.compute_energy()/model.num_particles;
			// std::cerr<<sweep <<"\t"<<E<<std::endl;
			model.dump(filename);
		}

	}
	printf("cell %s\n", model.cell?"true":"false");


	
}
