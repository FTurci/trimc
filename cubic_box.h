
#ifndef _CUBICBOX_
#define _CUBICBOX_ 

#include "box.h"
#include <vector>
using namespace Eigen;

class cubicBox : public Box {

	public:
		double side;

		cubicBox() : Box(){

		}

		cubicBox(double s,std::vector<double> c) : Box(s,s,s,c)  {
			side = s;
			printf(":: Cubix box with side %g and centre %g %g %g.\n",side, centre[0], centre[1], centre[2]);
		}
};

#endif