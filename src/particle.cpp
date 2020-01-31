#include "particle.h"


Particle::Particle(unsigned int i,Vector3d position){
	index = i;
	pos = position;
}

Particle::Particle(unsigned int i,double a, double b, double c){
	index = i;
	pos[0]=a;
	pos[1]=b;
	pos[2]=c;
}