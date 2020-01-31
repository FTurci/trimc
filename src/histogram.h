#ifndef _HISTOGRAM_
#define _HISTOGRAM_ 
#include <string>
#include <vector>

class Histogram{

	public:
		double dx;
		double nbins;

		double xmin;
		double xmax;

		std::vector<double> counts;
		std::vector<double> density;
		std::vector<double> edges;
		std::vector<double> centres;


		Histogram(){};
		Histogram(double min,double max, int numbins);

		void accumulate(double value);
		void save(std::string filename);
		void pdf();

};
#endif