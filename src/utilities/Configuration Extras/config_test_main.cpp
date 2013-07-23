/* basic unit test jig for gridpack::utilities::Configuration */

#include <iostream>
#include <string>
using namespace std;
#ifdef USE_MPI
#include <mpi.h>
#endif 
#include "configuration.hpp"
using namespace gridpack::utilities;

int main(int argc, char ** argv) {
#ifdef USE_MPI
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	cout << "Hello from rank " << world_rank << endl;
#endif
	Configuration config;
#ifdef USE_MPI
	if(world_rank == 0) {
#endif
		string src = "Documentation\\config_test.xml";
		cout << "Reading configuration " << src << endl;	
		config.enable_logging(&cout);
#ifdef USE_MPI
		if(config.open(src, MPI_COMM_WORLD)) {
#else
		if(config.open(src)) {
#endif
			// select a value by path name
		
			string start = config.get("Configuration.Time.Start", "DefaultStart");
			cout << "Start " << start << endl;
			// select a value by path name, example for missing value 
			string start0 = config.get("Configuration.Time.Start0", "DefaultStart");
			cout << "Start0 " << start0 << endl;
			// select a subtree by path 
			Configuration::Cursor * time = config.get_cursor("Configuration.Time");
			double step = -1.0;
			if(time->get("Step", & step)) // true if definition found
    	  		cout << "Step " << step << endl;

			std::vector<double> dv ;
			if(config.get("Configuration.Option.DefaultVelocity", &dv)) {
				cout << "DefaultVelocity: ";
				for(double e : dv) 
					cout << e << " ";
				cout << endl;
			}
		}
#ifdef USE_MPI
	}
	else {
		config.enable_logging(&cout);
		if(config.initialize(MPI_COMM_WORLD)) {
			string start = config.get("Configuration.Time.Start", "DefaultStart");
			cout << "Start[" << world_rank << "] " << start << endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	return 0;

}
