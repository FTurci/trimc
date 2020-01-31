#ifndef _INTERACTIONS_
#define _INTERACTIONS_ 

#include <vector>
#include <cstdio>
#include <algorithm>
#include <string>
#include "linearInterpolation.h"



class TwoBodyTabularInteraction  {
	public:
		std::vector<double> x; 
		std::vector<double> potential;
		int sizex;
		// region
		double xlo,xhi;
		double dx;
		double distance;

		TwoBodyTabularInteraction(){};
		void assign(std::string filename,int size,double lo, double hi);

		FL::LinearInterpolator<double> Interpolator;

		double at(double r);
		double at2(double r2);

};

class ThreeBodyTabularInteraction  {
	public:
		std::vector<double> positions; 
		std::vector<double> data;
		// sizes of the 3d tensor
		int sizex, sizey,sizez;
		// region
		double xlo,xhi,ylo,yhi,zlo,zhi;
		double dx,dy,dz;
	
		double rcut=1.4;

		ThreeBodyTabularInteraction(){};

		FL::LinearInterpolator<double> Interpolator;

		void assign(std::string filename,std::vector<int> sizes,std::vector<double> los, std::vector<double> his);


		double at(double rlo,double rmid,double rhi);
		

};

class Interaction{

	public:
		Interaction(){};
		TwoBodyTabularInteraction twobody;
		ThreeBodyTabularInteraction threebody;

};

class MultiBodyTabularInteraction: public Interaction{
	public:

		MultiBodyTabularInteraction(){};
		
		
		double at2(double r);
		
		double at(double rlo, double rmid, double rhi);
		double at(double r);
};

class LJInteraction: public Interaction {
	public:

		LJInteraction(){};
		LJInteraction(double eps, double s, double rc, bool shift_flag);

		void assign(double eps, double s, double rc, bool shift_flag);
		double epsilon;
		double sigma,sigma2;
		double rcut,rcut2;
		double Ucut,U;
		bool shift;

		double rinv2,rinv6,rinv12;

		double at(double r);
		double at2(double r2);
};

class HSInteraction : public Interaction{
	public:

		HSInteraction(){};
		HSInteraction(double sigma);
		double sigma,sigma2;
		void assign(double sigma);
		double at(double r);
		double at2(double r2);
		double pseudo(double r2);
		
};



#endif