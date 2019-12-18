#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include <Eigen/Dense>
#include <vector>
using namespace Eigen;


bool check_edge(Vector3d position, Vector3d edge){
	if (position[0]<edge[0] && position[1]<edge[1] && position[2]<edge[2]) {
		return true;
	}
	else return false;
}
class Structure{
	public:
		Structure(){};
		unsigned int rows,columns,stacks;
		std::vector<Vector3d> r;

		// void latticebuilder(int X, int Y, int Z, int *xx, int *yy, int *zz, int type)
		// {
		// 	int v = 2, x, y, z, xmax, ymax,arg;
		// 	float a, b, c;
		// 	xmax = 2*X; 
		// 	ymax = 2*Y; 
		// 	zmax = 2*Z;
		// 	for (z = 0; z <= zmax; z++) 
		// 		{ 
		// 			for (y = 0; y <= ymax; y++) 
		// 			{ 

		// 				for (x = 0; x <= xmax; x++)
		// 				zmax, arg;
		// 				{
		// 				a = x%v; b = y%v; c = z%v;
		// 				 arguments to produce BCC (8), SC (6) and FCC (12) lattices  
		// 				if (type == 8)
		// 					arg = (!a && !b && !c)||(a && b && c);
		// 				else if (type == 6)
		// 					arg = (!a && !b && !c); 
		// 				else
		// 					arg = (!a && !b && !c) ||(a && b && !c) ||(!a && b && c) ||(a && !b && c);
		// 				if (arg) 
		// 					{ 
		// 						*xx++ = x-X; 
		// 						*yy++ = y-Y; *zz++ = z-Z;
							
		// 					} 
		// 				}
		// 			}	 
		// 		}
			
		// }


};

class BccStructure: public Structure{

	public:
		using Structure::Structure;

		void build(Box &box,double rho,int take=-1){

			double constant = pow(2/rho,1./3.);

			rows = int(floor(box.sides[0]/constant));
			columns = int(floor(box.sides[1]/constant));
			stacks = int(floor(box.sides[2]/constant));

			double c0_x = -constant, c0_y = -constant, c0_z = -constant;
			double  c1_x, c1_y, c1_z = -constant/2;

			
			int maxn =  2*rows*columns*stacks ;

		


			for (int j = 0; j < stacks; j ++)   //for loops to produce the 3d points
			{		
					Vector3d v0,v1;

					c0_z += constant;
					c1_z += constant;
					v0[2] = c0_z-box.sides[2]*0.5+box.centre[2];
					v1[2] = c1_z-box.sides[2]*0.5+box.centre[2];

					c0_x = -constant;
				for (int i = 0; i < columns; i ++)
				{	
					c0_x += constant;
					c1_x = c0_x + constant/2;
					
					v0[0] = c0_x-box.sides[0]*0.5+box.centre[0];
					v1[0] = c1_x-box.sides[0]*0.5+box.centre[0];

					c0_y = -constant;
					for (int n = 0; n < rows; n ++)
					{
						c0_y += constant;
						c1_y = c0_y + constant/2;
						v0[1] = c0_y-box.sides[1]*0.5+box.centre[1];
						v1[1] = c1_y-box.sides[1]*0.5+box.centre[1];
						r.push_back(v0);
						r.push_back(v1);	
					}
				}
				}
		
			if(take!=-1) r.erase(r.begin()+take, r.end());	
		}

	
};

class FccStructure: public Structure{

	public:
		using Structure::Structure;

		void build(Box &box,double rho,int take=-1){

			double constant = pow(4/rho,1./3.);

			double c0_x = -constant, c0_y = -constant, c0_z = -constant;
			double  c1_x, c1_y, c1_z = -constant/2;

			rows = int(ceil(box.sides[0]/constant));
			columns = int(ceil(box.sides[1]/constant));
			stacks = int(ceil(box.sides[2]/constant));

			int maxn =  4*rows*columns*stacks ;
			printf(":: FCC %d x %d x %d, total %d particles\n",rows, columns,stacks, maxn );
			printf(":: FCC lattice constant %g\n", constant);

				// change box sides accordingly
			box.sides[0] = rows*constant;
			box.sides[1] = columns*constant;
			box.sides[2] = stacks*constant;

			printf(":: Adjusted box %g %g %g.\n", box.sides[0], box.sides[1],box.sides[2]);
			for (int j = 0; j < stacks; j ++) // for loops to print the atoms coordinates
			{	
				Vector3d v0,v1,v2,v3;
		        c0_z += constant;
		        c1_z += constant;
			    c0_x = -constant;
				
				for (int i = 0; i < columns; i ++)
				{
				        c0_x += constant;
				        c1_x = c0_x + constant/2;
				        // x0_values[i] = c0_x;
				        // x1_values[i] = c1_x;

				        c0_y = -constant;
				        for (int n = 0; n < rows; n ++)
				        {
				                c0_y += constant;
				                c1_y = c0_y + constant/2;
				                // y0_values[n] = c0_y;
				                // y1_values[n] = c1_y;
		                		
		                		v0[0]=c0_x-box.sides[0]*0.5+box.centre[0];
		                		v0[1]=c0_y-box.sides[1]*0.5+box.centre[1];
		                		v0[2]=c0_z-box.sides[2]*0.5+box.centre[2];

		                		v1[0]=c1_x-box.sides[0]*0.5+box.centre[0];
		                		v1[1]=c1_y-box.sides[1]*0.5+box.centre[1];
		                		v1[2]=c0_z-box.sides[2]*0.5+box.centre[2];

		                		v2[0]=c0_x-box.sides[0]*0.5+box.centre[0];
		                		v2[1]=c1_y-box.sides[1]*0.5+box.centre[1];
		                		v2[2]=c1_z-box.sides[2]*0.5+box.centre[2];

		                		v3[0]=c1_x-box.sides[0]*0.5+box.centre[0];
		                		v3[1]=c0_y-box.sides[1]*0.5+box.centre[1];
		                		v3[2]=c1_z-box.sides[2]*0.5+box.centre[2];

		                		Vector3d edge(3);
		                		for (int k =0; k<3;k++) edge[k]=box.sides[k]*0.5+box.centre[k];

		                		if (check_edge(v0,edge)) r.push_back(v0);
								if (check_edge(v1,edge)) r.push_back(v1);
								if (check_edge(v2,edge)) r.push_back(v2);
								if (check_edge(v3,edge)) r.push_back(v3);
		        		}
				}
			
				
				
			}
			printf(" Total r %d\n",r.size());
			if(take!=-1) r.erase(r.begin()+take, r.end());
			
		}


	
};


#endif // _STRUCTURES_