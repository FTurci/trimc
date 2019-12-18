#ifndef _PARTICLE_H
#define _PARTICLE_H 

#include <Eigen/Dense>
using namespace Eigen;

class Particle{

	public:
		//! position
		Vector3d pos;
		//! Default constructor
		Particle();
		//! Constructor with position
		Particle(unsigned int i,Vector3d position);
		//! Constructor with position values
		Particle(unsigned int i,double a, double b, double c);


		unsigned int index; //! Particle index
		unsigned int cell; //! index of the containing cell
		unsigned int posCell;

		double energy;
};
#endif