#ifndef _SYSTEM_
#define _SYSTEM_ 
#include <random>
#include <vector>
#include "particle.h"
#include "cubic_box.h"
#include <string>
#include "interactions.h"
#include "cells.h"
#include <vector>
#include "system.h"
#include <Eigen/Dense>
#include <random>
#include <fstream>
#include "utils.h"
#include "structures.h"


using namespace Eigen;



template <class T> 
class System {
	public:
		cubicBox box;
		bool cell;

		int num_particles;
		double max_disp2=0;

		std::vector<Particle> particles;
		std::vector<Particle> copy_particles;
		std::vector<double> distances;


		CellList cells; 

		T interaction;
		double scale;

		double energy;
		//random generator and uniform distribution
    	std::mt19937 gen;
		std::uniform_real_distribution<> dis;

		std::uniform_real_distribution<> dice;

		std::uniform_int_distribution<> picker;

		System(){};
		System(double side, int nparticles);

		double	minimum_distance();
		void randomise();

		inline void move(unsigned int);

		void copy();
		void restore();

		void restore(unsigned int);
		inline double compute_energy();

		void build_cells(double interaction_range);

		double particle_energy(const Particle& p, bool threebody);
		double pseudo_particle_energy(const Particle& p,bool pseudo=false);

		inline void local_MC_step(double beta, bool threebody=false);
		inline void sweep (double beta,double stepsize=0.1,bool threebody=false);

		void anneal(unsigned int cycles, double high, double low);
		// void read_interaction_table(std::string filename, double rc);

		void dump(std::string filename);

		void crystal(std::string type,double rho, int take=-1);
		void read_last_from_xyz(std::string filename);
};



template <class T> System <T>::System(double side, int nparticles){
	// build box
	std::vector<double> c (3,0);
	cubicBox b(side,c);

	this->box = b;
	scale = 0.1;
	//initalise random number generator
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	gen.seed(0);//rd());
	
	dis = std::uniform_real_distribution<> (box.centre[0]-side*0.5, box.centre[0]+side*0.5);
	
	// set the dice
	dice = std::uniform_real_distribution<> (0,1);


	printf(":: Appending %d particles...", nparticles);
	for (int i=0; i<nparticles; ++i){
		// particle with three random numbers
		Particle p(i,dis(gen),dis(gen),dis(gen));
		this->particles.push_back(p);
		this->copy_particles.push_back(p);
	}
	printf(" done.\n");
	this-> num_particles = this->particles.size();

	picker = std::uniform_int_distribution<> (0,num_particles-1);


}

template <class T> 
double System<T>::minimum_distance(){

	Vector3d dr;
	double minimum = this->box.side*2;
	double distance;

	std::vector<double> dists;

	for(int i=0; i<num_particles-1;++i){
		for(int j=i+1; j<num_particles;++j){
			double d=0;
			double dk;
			dr = particles[i].pos-particles[j].pos;

			// // periodic boundaries
			for (int k=0; k<3;++k){
				if (dr[k]>0.5*this->box.side) dr[k]-=this->box.side;
				if (dr[k]<-0.5*this->box.side) dr[k]+=this->box.side;
			}
			distance = dr.norm();
			dists.push_back(distance);
			if (distance<minimum) minimum = distance;
		}

	}
	this->distances=dists;
	return minimum;
}

template <class T> 
void System<T>::build_cells(double interaction_range){
	// Hard-coded 3 dimensions
	cells.setDimension(this->box.dimension);
	printf(":: Interaction range %g\n", interaction_range);
	cells.initialise(this->box.sides,this->box.centre, interaction_range);

	// Associate particles with cells
	for (int i = 0; i <num_particles;++i){
		// printf(":: particle %d\n", i);
		particles[i].cell = cells.getCell(particles[i]);
		// printf(":: cell %d\n",particles[i].cell);
	// 	// Update cell list.
        cells.initCell(particles[i].cell, particles[i]);
	}

}
template <class T> 
void System<T>::randomise(){
	for(int i=0; i<num_particles;++i){
		for(int k=0; k<3;++k)
			this->particles[i].pos[k] = this->dis(this->gen);
	}

}
template <class T> 
void System<T>::copy(){
	for(int i=0; i<num_particles;++i){
		for(int k=0; k<3;++k)
			this->copy_particles[i].pos[k] = this->particles[i].pos[k] ;
	}

}

template <class T> 
void System<T>::restore(unsigned int i){
		for(int k=0; k<3;++k)
			this->particles[i].pos[k] = this->copy_particles[i].pos[k] ;
}
template <class T> 
void System<T>::restore(){
	for(int i=0; i<num_particles;++i){
		for(int k=0; k<3;++k)
			this->particles[i].pos[k] = this->copy_particles[i].pos[k] ;
	}

}

template <class T> 
void System<T>::move(unsigned int index){
	for(int k=0; k<box.dimension;++k){
		// store the  position before the move
		copy_particles[index].pos[k] = particles[index].pos[k];

		particles[index].pos[k]+=scale*(dice(gen)-0.5)*2;

		if (particles[index].pos[k]>box.centre[k]+box.sides[k]*0.5)
			particles[index].pos[k]-=box.sides[k];
		else if (particles[index].pos[k]<box.centre[k]-box.sides[k]*0.5)
			particles[index].pos[k]+=box.sides[k];

	}

}
template <class T> 
void System<T>::local_MC_step(double beta,bool threebody){
	// compute energy

	
	double oldenergy, newenergy;
	unsigned int pick = picker(this->gen);

	// printf("old pos %g %g %g\n", particles[pick].pos[0],particles[pick].pos[1],particles[pick].pos[2]);

	if(cell) {
		// printf("computing old\n");
		oldenergy = particle_energy(particles[pick],threebody);
		// printf("oldenergy %g\n",oldenergy);
	}

	else
	 oldenergy =compute_energy();

	move(pick); 
	// printf("new pos %g %g %g\n", particles[pick].pos[0],particles[pick].pos[1],particles[pick].pos[2]);

	if (cell) newenergy = particle_energy(particles[pick],threebody);
	else
	 newenergy = compute_energy();

	energy = newenergy;

	double deltaE  = newenergy-oldenergy;
	
	// printf("E: o %g  n %g Δ  %g\n", oldenergy,newenergy, deltaE);



	if (dice(this->gen)< exp(-beta*deltaE)){
		// printf("accept\n");
	}
	else{
		// printf("    reject\n");
		energy =oldenergy;
		restore(pick);
	}

	if (cell){
		    // Calculate the particle's cell index.
    	unsigned int newCell = cells.getCell(particles[pick]);
	    // Update cell lists if necessary.
	    if (particles[pick].cell != newCell){
	        cells.updateCell(newCell, particles[pick], particles);
	    }
	    particles[pick].energy = energy;
	}
}

template <class T> 
void System<T>::sweep(double beta, double step_size,bool threebody){
	scale= step_size;
	for(int i=0; i<num_particles; i++){
		local_MC_step(beta,threebody);
	}

}

template <class T> 
double System<T>::compute_energy(){

	Vector3d dr;

	double d,d2;
	double e=0;

	for(int i=0; i<num_particles-1;++i){
		// printf("i %d of %d\n",i, num_particles);
		for(int j=i+1; j<num_particles;++j){
			// printf("j %d\n",j);
			for(int k=0; k<3;++k){			
				dr[k] = this->particles[i].pos[k] - this->particles[j].pos[k] ;
				if (dr[k]>0.5*this->box.sides[k]) dr[k]-=this->box.sides[k];
				if (dr[k]<-0.5*this->box.sides[k]) dr[k]+=this->box.sides[k];
			}

			d2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
			// printf("    i %d j %d d %g\n",i,j,d);
			e+= this->interaction.at2(d2);

			}
	
	}

	return e;
}

template <class T> 
double System<T>::particle_energy(const Particle &p, bool threebody){

	Vector3d dr;
	double dr2;
	double e=0;
	int i = p.index;
	int cid = cells.getCell(p);

	// to help the 3-body calculation
	std::vector<int> neighs(100);
	std::vector<double> distances(100);
	std::vector<double> triangle(3);
	int num_neighs=0;
	double tb;

	// loop over all cells (including cid)
	for (int c=0; c<cells[cid].neighbours.size(); ++c){
		int cell = cells[cid].neighbours[c];
		// printf("here we are %d particle %d tally %d\n",cell,i,cells[c].tally);

		for (int n=0; n<cells[cell].tally;++n){
			int j = cells[cell].particles[n];
			// printf("here i %d j %d\n",i,j);
			if (i!=j){

				for(int k=0; k<3;++k)
				{			
				dr[k] = this->particles[i].pos[k] - this->particles[j].pos[k] ;
				if (dr[k]>0.5*this->box.sides[k]) dr[k]-=this->box.sides[k];
				if (dr[k]<-0.5*this->box.sides[k]) dr[k]+=this->box.sides[k];
				}

				dr2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
				e+= this->interaction.twobody.at2(dr2);

				if (threebody){
					if(interaction.twobody.distance<interaction.twobody.xhi){

						neighs[num_neighs]=j;
						// printf("r2 %g %g\n",sqrt(dr2),interaction.twobody.distance);
						distances[num_neighs]=interaction.twobody.distance;
						num_neighs++;
						
					}

				}

			}		
		}
	}
	// now get the threebody term from all triples {i j1 j2}, where j1 and j2 are neighbours of i. In particular, obtain triangles whose sides are sorted

	if (threebody)
	{
		int j1,j2;
		double e1,e2,e3;
		if (num_neighs>=2){
			for(int n1=0;n1<num_neighs-1;n1++){
					j1 = neighs[n1];
					for(int n2=n1+1;n2<num_neighs;++n2){
						j2 = neighs[n2];
						// reassign (needed because I am sorting inplace)
						triangle[0]=distances[n1];
						triangle[1]=distances[n2];
						//j1-j2
						for(int k=0; k<3;++k)
						{			
							dr[k] = this->particles[j1].pos[k] - this->particles[j2].pos[k] ;
							if (dr[k]>0.5*this->box.sides[k]) dr[k]-=this->box.sides[k];
							if (dr[k]<-0.5*this->box.sides[k]) dr[k]+=this->box.sides[k];
						}

						triangle[2] = sqrt1(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]); 

						e3 = this->interaction.at(triangle[2]);
						e2 = this->interaction.at(triangle[1]);
						e1 = this->interaction.at(triangle[0]);

						std::stable_sort (triangle.begin(), triangle.end());
						// printf("tri after %g %g %g\n", triangle[0],triangle[1],triangle[2]);

						// actual 3body interaction
						// if (triangle[2]<interaction.twobody.xhi && triangle[0]>0.97)
						if (triangle[2]<2.5 && triangle[0]>1.0)
						{ 

								// remove two bodies
								e-=e1;
								e-=e2;
								e-=e3;
								tb=this->interaction.at(triangle[0], triangle[1], triangle[2]);
								e+=tb;

				
							}


					}

			}


	

		}
	}
	// printf("end body\n\n");
	// printf("energy %g\n",e);

	return e;
}

template <class T> 
void System<T>::anneal(unsigned int cycles, double high, double low){
	double temp,disp;
	disp=0.1;
	for (unsigned int c=0; c<cycles;c++){
		temp = high+c*(low-high)/cycles;
		disp *= 0.999;
		sweep(temp,disp,true);
		printf("%g %g\n",temp,disp);
	}


}

template <class T> 
void System<T>::dump(std::string filename){
	std::ofstream fout(filename.c_str(), std::fstream::out | std::fstream::app);

	fout<<num_particles<<"\nAtoms\n";
	for (int i=0; i<num_particles;i++){
		fout<<"A "<<particles[i].pos[0]<<" "<<particles[i].pos[1]<<" "<<particles[i].pos[2]<<"\n";
	}

	fout.close();

}
template <class T> 
void System<T>::crystal(std::string type,double rho, int take){
	if (type=="fcc"){

		printf(":: Building FCC structure\n");
		FccStructure structure;

		structure.build(box,rho,take);
		particles.clear();
		for (int i; i<structure.r.size();i++){
			Particle p(i,structure.r[i][0],structure.r[i][1],structure.r[i][2]);
			particles.push_back(p);
			copy_particles.push_back(p);
		}
	num_particles = this->particles.size();
	picker = std::uniform_int_distribution<> (0,num_particles-1);
	printf(":: Number of particles %d\n", num_particles);
	}

}

template <class T> 
void System<T>::read_last_from_xyz(std::string filename){
	std::ifstream fin;

	int N;
	fin.open(filename);
	fin >>N;
	fin.close();
	// restart the reading
	fin.open(filename);
	std::vector<std::vector<double> > positions(N,std::vector<double> (3, 0));
	
	double a ,b,c;
	char type;
	char dummy[256];

	int counter = 0;
	while (fin.good()){
		for (int i=0; i<N+2; i++){
			if (i<2){ 
				fin>>dummy;
				if (fin.eof()) break;
			}
			else{
				fin>>type>>a>>b>>c;
				positions[i-2][0]=a;
				positions[i-2][1]=b;
				positions[i-2][2]=c;
			}
		}
		counter++;
	}
	fin.close();

	printf("Read the xyz file: last conf\n");
	
	this->particles.clear();
	this->copy_particles.clear();
	for(int i= 0; i<N;i++){
		Particle p(i,positions[i][0],positions[i][1],positions[i][2]);
		this->particles.push_back(p);
		this->copy_particles.push_back(p);	
	}
	this-> num_particles = this->particles.size();
	picker = std::uniform_int_distribution<> (0,num_particles-1);
}

// void Tabular::read_interaction_table(std::string filename, double rc){
// 	double r,u;
// 	std::ifstream fin;

// 	fin.open(filename.c_str(), std::ifstream::in);
// 	while(fin.good()){
// 		fin>>r>>u;
// 		printf("r %g, u %g",r,u);
		
// 		this->interaction.rs.push_back(r);
// 		this->interaction.us.push_back(u);
// 	}
	
// 	this->interaction.rcut = rc;
// 	this->interaction.shortDistanceValue = 10*this->interaction.us[0];
// 	this->interaction.largeDistanceValue =0;
// 	this->interaction.rmin = this->interaction.rs[0];

// 	printf("Setup: %g %g %g\n",this->interaction.rcut ,this->interaction.shortDistanceValue ,this->interaction.largeDistanceValue);
// 	for(int i=0; i<interaction.rs.size(); i++){
// 		printf("%g %g \n",interaction.rs[i], interaction.us[i]);
// 	}
// }



// No need to call this TemporaryFunction() function,
// it's just to avoid link error.void TemporaryFunction ()

// class Tabular : public System{
// 	public:
// 		TabularInteraction interaction;
// 		// inherit constructors from System
// 		using System::System;


// 		void read_interaction_table(std::string filename, double rc);
// };

// class LennardJonesium : public System{
// 	public:
		
// 		// inherit constructors from System
// 		using System::System;

// 		LennardJonesium(double side, int npart,double eps, double s, double rc, bool shift_flag) :LennardJonesium(side, npart){
// 				interaction.assign( eps,  s,  rc,  shift_flag);
// 		};

// 		LJInteraction interaction;
// };
#endif

