#ifndef _BOX_H
#define _BOX_H 
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;

class Box {
  public:

  	unsigned int dimension;

	double sidex;
  	double sidey;
  	double sidez;
  	
  	std::vector<double> centre;
  	std::vector<double> sides;


  	Box(){};
	

	Box(double Lx,double Ly, double Lz, std::vector<double> c);
	Box(double Lx,double Ly, std::vector<double> c);


};

#endif