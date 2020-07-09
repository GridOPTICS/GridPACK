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
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/configuration/configuration.hpp"
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
     * Reinitialize calculation from data collections
     */
    void reload();

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
    void writeCABus();
    void writeCABranch();
    void writeHeader(const char *msg);
    std::vector<std::string> writeBusString(const char *signal = NULL);
    std::vector<std::string> writeBranchString(const char *signal = NULL);

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
     * Export final configuration to PSS/E formatted file
     * @param filename name of file to store network configuration
     */
    void exportPSSE33(std::string &filename);
    
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
     * Set voltage limits on all buses
     * @param Vmin lower bound on voltages
     * @param Vmax upper bound on voltages
     */
    void setVoltageLimits(double Vmin, double Vmax);

    /**
     * Check to see if there are any voltage violations in the network.
     * @param area only check for violations in specified area
     * @return true if no violations found
     */
    bool checkVoltageViolations();
    bool checkVoltageViolations(int area);

    /**
     * Set "ignore" parameter on all buses with violations so that subsequent
     * checks are not counted as violations
     */
    void ignoreVoltageViolations();

    /**
     * Clear "ignore" parameter on all buses
     */
    void clearVoltageViolations();

    /**
     * Check to see if there are any line overload violations in
     * the network. The last call checks for overloads on specific lines.
     * @param area only check for violations in specified area
     * @param bus1 original index of "from" bus for branch
     * @param bus2 original index of "to" bus for branch
     * @param tags line IDs for individual lines
     * @param violations true if violation detected on branch, false otherwise
     * @return true if no violations found
     */
    bool checkLineOverloadViolations();
    bool checkLineOverloadViolations(int area);
    bool checkLineOverloadViolations(std::vector<int> &bus1,
        std::vector<int> &bus2, std::vector<std::string> &tags,
        std::vector<bool> &violations);

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
     * Check to see if there are any Q limit violations in the network
     * @param area only check for violations in specified area
     * @return true if no violations found
     */
    bool checkQlimViolations();
    bool checkQlimViolations(int area);

    /**
     * Clear changes that were made for Q limit violations and reset
     * system to its original state
     */
    void clearQlimViolations();

    /**
     * Reset voltages to values in network configuration file
     */
    void resetVoltages();

    /**
     * Scale generator real power. If zone less than 1 then scale all
     * generators in the area.
     * @param scale factor to scale real power generation
     * @param area index of area for scaling generation
     * @param zone index of zone for scaling generation
     */
    void scaleGeneratorRealPower(double scale, int area, int zone);

    /**
     * Scale load power. If zone less than 1 then scale all
     * loads in the area.
     * @param scale factor to scale load real power
     * @param area index of area for scaling load
     * @param zone index of zone for scaling load
     */
    void scaleLoadPower(double scale, int area, int zone);

    /**
     * Return the total real power load for all loads in the zone. If zone
     * less than 1, then return the total load real power for the area
     * @param area index of area
     * @param zone index of zone
     * @return total load
     */
    double getTotalLoadRealPower(int area, int zone);

    /**
     * Return the current real power generation and the maximum and minimum total
     * power generation for all generators in the zone. If zone is less than 1
     * then return values for all generators in the area
     * @param area index of area
     * @param zone index of zone
     * @param total total real power generation
     * @param pmin minimum allowable real power generation
     * @param pmax maximum available real power generation
     */
    void getGeneratorMargins(int area, int zone, double *total, double *pmin,
        double *pmax);

    /**
     * Reset power of loads and generators to original values
     */
    void resetPower();

    /**
     * Write real time path rating diagnostics
     * @param src_area generation area
     * @param src_zone generation zone
     * @param load_area load area
     * @param load_zone load zone
     * @param gen_scale scale factor for generation
     * @param load_scale scale factor for loads
     * @param file name of file containing diagnostics
     */
    void writeRTPRDiagnostics(int src_area, int src_zone, int load_area,
        int load_zone, double gen_scale, double load_scale, const char *file);

    /**
     * Get strings documenting contingency failures. Strings are *not*
     * terminated with a carriage return
     */
    std::vector<std::string> getContingencyFailures();

    /**
     * User rate B parameter for line overload violations
     * @param flag if true, use RATEB parameter
     */
    void useRateB(bool flag);

    /**
     * Suppress all output from power flow module
     * @param flag if true, suppress printing
     */
    void suppressOutput(bool flag);

#ifdef USE_GOSS
    /**
     * Set GOSS client if one already exists
     * @param client pointer to existing GOSS client
     * @param simID simulation ID
     */
    void setGOSSClient(gridpack::goss::GOSSClient &client, std::string simID);
#endif
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

    // qlim enforce flag
    int p_qlim;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<PFNetwork> > p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<PFNetwork> > p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;

    // string containing current contingency name
    std::string p_contingency_name;

    // Flag to suppress all printing to standard out
    bool p_no_print;

#ifdef USE_GOSS
    gridpack::goss::GOSSClient p_goss_client;

    std::string p_simID;
#endif
};

} // powerflow
} // gridpack
#endif
