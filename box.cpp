#include "box.h"
#include <Eigen/Dense>
using namespace Eigen;


//! 3D constructor
Box::Box(double Lx,double Ly, double Lz, std::vector<double> c){
	// sidex = Lx;
	// sidey = Ly;
	// sidez = Lz;
	centre = c;

	dimension = 3;

	sides.resize(dimension);
	sides[0] = Lx;
	sides[1] = Ly;
	sides[2] = Lz;

}

//! 2D constructor
Box::Box(double Lx,double Ly,  std::vector<double> c){
	sidex = Lx;
	sidey = Ly;

	centre[0] = c[0];
	centre[1] = c[1];

	dimension = 2;
	sides.resize(dimension);

	sides[0] = Lx;
	sides[1] = Ly;
}