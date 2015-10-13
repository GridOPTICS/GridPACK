/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_module_h_
#define _pf_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "pf_factory_module.hpp"

namespace gridpack {
namespace powerflow {

// Structs that are used for some applications

struct pathBranch{
  int fromBus;
  int toBus;
  std::string branchID;
};

struct stressBus{
  int busID;
  std::string genID;
  double participation;
};

struct pathStress{
  std::string name;
  std::vector<pathBranch> path;
  std::vector<stressBus> sourceArea;
  std::vector<stressBus> sinkArea;
};

// Contingency types
enum ContingencyType{Generator, Branch};

// Struct that is used to define a collection of contingencies

struct Contingency
{
  int p_type;
  std::string p_name;
  // Line contingencies
  std::vector<int> p_from;
  std::vector<int> p_to;
  std::vector<std::string> p_ckt;
  // Status of line before contingency
  std::vector<bool> p_saveLineStatus;
  // Generator contingencies
  std::vector<int> p_busid;
  std::vector<std::string> p_genid;
  // Status of generator before contingency
  std::vector<bool> p_saveGenStatus;
};

// Calling program for powerflow application

class PFAppModule
{
  public:
    /**
     * Basic constructor
     */
    PFAppModule(void);

    /**
     * Basic destructor
     */
    ~PFAppModule(void);

    /**
     * Read in and partition the powerflow network. The input file is read
     * directly from the Powerflow block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a PFNetwork object. This should not have any
     * buses or branches defined on it.
     * @param config point to open configuration file
     */
    void readNetwork(boost::shared_ptr<PFNetwork> &network,
                     gridpack::utility::Configuration *config);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Execute the iterative solve portion of the application using a hand-coded
     * Newton-Raphson solver
     * @return false if an error was caught in the solution algorithm
     */
    bool solve();

    /**
     * Execute the iterative solve portion of the application using a library
     * non-linear solver
     * @return false if an error was caught in the solution algorithm
     */
    bool nl_solve();

    /**
     * Write out results of powerflow calculation to standard output
     * Separate calls for writing only data from buses or branches
     * @param signal tell underlying write what records to print
     */
    void write();
    void writeBus(const char* signal);
    void writeBranch(const char* signal);

    /**
     * Redirect output from standard out
     * @param filename name of file to write results to
     */
    void open(const char *filename);
    void close();

    /**
     * Print string. This can be used to direct output to the file opened using
     * the open command
     * @param buf string to be printed
     */
    void print(const char *buf);

    /**
     * Save results of powerflow calculation to data collection objects
     */
    void saveData();
    
    /**
     * Set a contingency
     * @param event data describing location and type of contingency
     * @return false if location of contingency is not found in network
     */
    bool setContingency(Contingency &event);

    /**
     * Return system to the state before the contingency
     * @param event data describing location and type of contingency
     * @return false if location of contingency is not found in network
     */
    bool unSetContingency(Contingency &event);

    /**
     * Check to see if there are any voltage violations in the network
     * @param minV maximum voltage limit
     * @param maxV maximum voltage limit
     * @return true if no violations found
     */
    bool checkVoltageViolations(double Vmin, double Vmax);
    bool checkVoltageViolations(int area, double Vmin, double Vmax);

    /**
     * Set "ignore" parameter on all buses with violations so that subsequent
     * checks are not counted as violations
     * @param minV maximum voltage limit
     * @param maxV maximum voltage limit
     */
    void ignoreVoltageViolations(double Vmin, double Vmax);

    /**
     * Clear "ignore" parameter on all buses
     */
    void clearVoltageViolations();

    /**
     * Check to see if there are any line overload violations in
     * the network
     * @return true if no violations found
     */
    bool checkLineOverloadViolations();
    bool checkLineOverloadViolations(int area);

    /**
     * Set "ignore" paramter on all lines with violations so that subsequent
     * checks are not counted as violations
     */
    void ignoreLineOverloadViolations();

    /**
     * Clear "ignore" parameter on all lines
     */
    void clearLineOverloadViolations();

    /**
     * Reset voltages to values in network configuration file
     */
    void resetVoltages();
  private:

    // pointer to network
    boost::shared_ptr<PFNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<PFFactoryModule> p_factory;

    // maximum number of iterations
    int p_max_iteration;

    // convergence tolerance
    double p_tolerance;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<PFNetwork> > p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<PFNetwork> > p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;
};

} // powerflow
} // gridpack
#endif
