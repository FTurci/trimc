#include "interactions.h"
#include <cmath>
#include <limits>
#include <fstream>
#include <iterator>
#include <iostream>
#include "utils.h"

// ### SHORTCUTS
double INF = std::numeric_limits<double>::infinity();

inline long int idx(int x,int y,int z,int sizex,int sizey,int sizez){
	// printf("x,y,z %d %d %d | %d %d %d\n", x,y,z, sizex,sizey,sizez);
	return z*(sizex*sizey)+sizey*x+y;
}

// ### TABULATED POTENTIALS

void TwoBodyTabularInteraction::assign(std::string filename,int size,double lo, double hi){
	// read data: 1d vector in output
	std::ifstream is(filename.c_str());
	if (!is.good() ){
		printf("[ERROR:TwoBodyTabularInteraction::assign] File does not exist!\n");
		exit(0);
	}

  	std::istream_iterator<double> start(is), end;
 	std::vector<double> data(start, end);
 	is.close();
 	// /split the two columns in range and potential
 	if (x.size()==0){
	 	for (int i = 0; i < data.size(); i+=2)
	 	{
	 		this->x.push_back(data.at(i));
	 		this->potential.push_back(data.at(i+1));

	 		// printf("= %g %g | %g %g \n", data.at(i), data.at(i+1), x.back(), potential.back());
	 	}
 	} else printf("[ERROR:TwoBodyTabularInteraction::assign] Pre-existing x range!\n");

 	this->sizex = size;
 	this->xlo = lo;
 	this->xhi = hi;
 	this->dx = (hi-lo)/size;
}

double TwoBodyTabularInteraction::at(double r){
	if (r<=xlo)		 {
		double slope = (potential[1]-potential[0])/dx;
		double c = potential[0]-slope*x[0];
		// extrapolate linearly
		// return slope*x[0]+c;
		return 4./(r*r*r*r*r*r*r*r*r*r*r*r);
		// printf("[ERROR:TwoBodyTabularInteraction::at] Evaluating undefined potential!\n");
		// exit(0);
	}
	else if (r>=xhi) 	return 0;//this->potential.back();
	else{

		
		//using the syntax of the linearinterolation library
		// constructing the evaluation point
		FL::point<1,double> p;
		p.coords[0]=r;
		// get the index of the point to the left
		int idx0=floor((p.coords[0]-xlo)/dx);
		//next index
		int idx1 = idx0+1;
		// construct the two points
		FL::point<1,double> p0,p1;
		p0.coords[0] = x[idx0];

		p0.val = potential[idx0];

		p1.coords[0] = x[idx1];
		p1.val = potential[idx1];

		// interpolation	
		FL::point<1,double> result= this->Interpolator.Linear(p, p0,p1);

		// printf("%g %g\n",r,result.val);
		// std::cerr<<r<<" "<<result.val<<std::endl;

		if (result.val>100){

			printf("Errore %g %g\n", p0.coords[0],r);
		}
		return result.val;
	}
}

double TwoBodyTabularInteraction::at2(double r2){

	distance = sqrt1(r2);  
	return at(distance);
}


void ThreeBodyTabularInteraction::assign(std::string filename,std::vector<int> sizes,std::vector<double> los, std::vector<double> his){
	// read data: 1d vector in output
	std::ifstream is(filename.c_str());
  	std::istream_iterator<double> start(is), end;
 	std::vector<double> input_data(start, end);
 	this->data = input_data;
 	is.close();


 	printf(":: The threebody data contain %d elements.\n", data.size());


 	this->sizex = sizes[0];
 	this->sizey = sizes[1];
 	this->sizez = sizes[2];

 	printf(":: The shape is %d %d %d\n", sizex,sizey,sizez);

 	
 	this->xlo = los[0];
 	this->ylo = los[1];
 	this->zlo = los[2];

 	this->xhi = his[0];
 	this->yhi = his[1];
 	this->zhi = his[2];

 	printf(":: The region is (%g,%g) x (%g,%g) x (%g,%g)\n", xlo,xhi, ylo,yhi,zlo,zhi);

 	this->dx = (xhi-xlo)/sizex;
 	this->dy = (yhi-ylo)/sizey;
 	this->dz = (zhi-zlo)/sizez;

 	printf(":: The deltas are is %g %g %g\n", dx,dy,dz);

 	this->rcut = std::max(std::max(xhi,yhi),zhi);
}
double ThreeBodyTabularInteraction::at(double rlo,double rmid,double rhi){
	// use 3d interpolation
	
	double xpos,ypos,zpos;
	
	xpos = rlo;
	ypos = rmid;
	zpos = rhi;

	FL::point<3,double> p;

	p.coords[0]=xpos;
	p.coords[1]=ypos;
	p.coords[2]=zpos;

	// use x for the rows, y for the columns and z for the slices

	// find slice

	FL::point<3, double> corners[8];
	// first point
	int zindex0 = floor( (zpos-zlo)/dz);
	int yindex0 = floor( (ypos-ylo)/dy);
	int xindex0 = floor( (xpos-xlo)/dx);

	// printf("pos %g %g %.5f\n", xpos, ypos,zpos);	

	corners[0].coords[0] = xlo+xindex0*dx;
	corners[0].coords[1] = ylo+yindex0*dy;
	corners[0].coords[2] = zlo+zindex0*dz;
	// printf("corners done %d \n",idx( xindex0, yindex0, zindex0, sizex, sizey, sizez));	
	corners[0].val = data.at(idx( xindex0, yindex0, zindex0, sizex, sizey, sizez));

	int zindex1 = zindex0;
	int yindex1 = yindex0;
	int xindex1 = xindex0+1;

	corners[1].coords[0] = xlo+xindex1*dx;
	corners[1].coords[1] = ylo+yindex1*dy;
	corners[1].coords[2] = zlo+zindex1*dz;
	corners[1].val = data.at(idx( xindex1, yindex1, zindex1, sizex, sizey, sizez));

	int zindex2 = zindex0;
	int yindex2 = yindex0+1;
	int xindex2 = xindex0;


	corners[2].coords[0] = xlo+xindex2*dx;
	corners[2].coords[1] = ylo+yindex2*dy;
	corners[2].coords[2] = zlo+zindex2*dz;
	corners[2].val = data.at(idx( xindex2, yindex2, zindex2, sizex, sizey, sizez));

	int zindex3 = zindex0;
	int yindex3 = yindex0+1;
	int xindex3 = xindex0+1;

	corners[3].coords[0] = xlo+xindex3*dx;
	corners[3].coords[1] = ylo+yindex3*dy;
	corners[3].coords[2] = zlo+zindex3*dz;
	corners[3].val = data.at(idx( xindex3, yindex3, zindex3, sizex, sizey, sizez));

	int zindex4 = zindex0+1;
	int yindex4 = yindex0;
	int xindex4 = xindex0;



	corners[4].coords[0] = xlo+xindex4*dx;
	corners[4].coords[1] = ylo+yindex4*dy;
	corners[4].coords[2] = zlo+zindex4*dz;
	corners[4].val = data.at(idx( xindex4, yindex4, zindex4, sizex, sizey, sizez));


	int zindex5 = zindex0+1;
	int yindex5 = yindex0;
	int xindex5 = xindex0+1;


	corners[5].coords[0] = xlo+xindex5*dx;
	corners[5].coords[1] = ylo+yindex5*dy;
	corners[5].coords[2] = zlo+zindex5*dz;
	corners[5].val = data.at(idx( xindex5, yindex5, zindex5, sizex, sizey, sizez));


	int zindex6 = zindex0+1;
	int yindex6 = yindex0+1;
	int xindex6 = xindex0;

	corners[6].coords[0] = xlo+xindex6*dx;
	corners[6].coords[1] = ylo+yindex6*dy;
	corners[6].coords[2] = zlo+zindex6*dz;
	corners[6].val = data.at(idx( xindex6, yindex6, zindex6, sizex, sizey, sizez));

	int zindex7 = zindex0+1;
	int yindex7 = yindex0+1;
	int xindex7 = xindex0+1;


	corners[7].coords[0] = xlo+xindex7*dx;
	corners[7].coords[1] = ylo+yindex7*dy;
	corners[7].coords[2] = zlo+zindex7*dz;
	corners[7].val = data.at(idx( xindex7, yindex7, zindex7, sizex, sizey, sizez));

	FL::point<3,double> result= Interpolator.Trilinear(p, corners);
	// printf("3 body interaction %g\n", result.val);
	return result.val;
}

double MultiBodyTabularInteraction::at2(double r2){ 
	return twobody.at2(r2);
}

double MultiBodyTabularInteraction::at(double  rlo, double rmid, double rhi){ 
	return threebody.at( rlo,  rmid,  rhi);
}

double MultiBodyTabularInteraction::at(double  r){ 
	return twobody.at(r);
}

// #### Lennard-Jones

LJInteraction::LJInteraction(double eps, double s, double rc, bool shift_flag){

	epsilon = eps;
	sigma = s;
	rcut = rc;

	shift = shift_flag;
	Ucut = 4.0*epsilon*(pow(sigma/rcut,12)-pow(sigma/rcut,6));
}

void LJInteraction::assign(double eps, double s, double rc, bool shift_flag){

	epsilon = eps;
	sigma = s;
	sigma2 = s*s;
	rcut = rc;
	rcut2 = rc*rc;

	shift = shift_flag;
	if(shift)
        Ucut = 4.0*epsilon*(pow(sigma/rcut,12)-pow(sigma/rcut,6));
    else Ucut=0;
}

double LJInteraction::at(double r){
	// printf(" -- lj--\n");

	U=0;
	if (r<rcut){
		U = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
		if (shift)  U-=Ucut;
	}
	
	return  U;
}

double LJInteraction::at2(double r2){

	U=0;
	if (r2<rcut2){
		rinv2 = sigma2/r2;
		rinv6 = rinv2*rinv2*rinv2;
		rinv12 = rinv6*rinv6;
		// printf(" -- lj--\n");
		U = 4*epsilon*(rinv12-rinv6);
	    U-=Ucut;
	}	
	return  U;
}

// ### HARD SPHERES

HSInteraction::HSInteraction(double s){
	sigma = s;
	sigma2 = s*s;
}

void HSInteraction::assign(double s){

	sigma = s;
	sigma2 = s*s;
}

double HSInteraction::at(double r){

	if (r>sigma) return 0;
	else return INF;
}

double HSInteraction::at2(double r2){

		if (r2>sigma2) return 0;
		else {
			return INF;
		}
			
	
}
double HSInteraction::pseudo(double r2){
		if (r2< 1.020408163){
			double ur= pow(r2,25);

			double e= 134.552662342*(1.0/(ur)-sqrt(r2)/ur);
			// printf("e %g\n",e);
			return e;
			}
		else return 0;

}

