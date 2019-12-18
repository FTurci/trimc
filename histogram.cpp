#include "histogram.h"
#include <fstream>
#include  <cmath>

Histogram::Histogram(double min,double max, int numbins){
		double span = max-min;
		dx = span/(double)numbins;
		nbins=numbins;

		xmin = min;
		xmax = max;

		edges.push_back(min);

		for(int i=0;i<numbins;++i){
			counts.push_back(0);
			density.push_back(0);
			edges.push_back(edges[i]+dx);
			centres.push_back(edges[i]+0.5*dx);
		}
}

void Histogram::accumulate(double value){
	// find closest bin
	if (value>xmax) printf(":: Value above range: not counted.\n");
	if (value<xmin) printf(":: Value below range: not counted.\n");

	double v = value-xmin;
	int index = floor(v/dx);
	counts[index]++;

}

void Histogram::save(std::string filename){

	std::ofstream fout(filename.c_str(), std::ios::out | std::ios::trunc);

	fout<<"# Histrogram table:\n# centres counts"<<std::endl;
	for(int i=0; i<counts.size(); ++i){
		fout<<centres[i]<<"\t"<<counts[i]<<"\t"<<density[i]<<std::endl;
	}
	fout.close();
}

void Histogram::pdf(){
	double totalcount=0;
	for (int k; k<counts.size(); ++k){
		totalcount+=counts[k];
	}
	for (int k; k<counts.size(); ++k){
		density[k]=counts[k]/dx/totalcount;
	}
}
