/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app_module.hpp
 * @author Yousu Chen, Bruce Palmer
 * @date   1/23/2015
 * @Last modified 1/23/2015
 *
 * @brief
 * @update Yousu Chen
 *         Adding functions of bad data dection, chi-square testing 
 * @date   2025-03-05
 *         Adding more functions to handle measurements more efficiently
 * @date   2025-04-02
 *
 *
 */
// -------------------------------------------------------------

#ifndef _se_app_module_h_
#define _se_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "se_factory_module.hpp"

namespace gridpack {
namespace state_estimation {

// Calling program for state estimation application
// Calling program for state estimation application
class SEAppModule
{
  public:
    /**
     * Basic constructor
     */
    //SEAppModule(gridpack::parallel::Communicator comm);
    SEAppModule(void);

    /**
     * Basic destructor
     */
    ~SEAppModule(void);

    /**
     * Get list of measurements from external file
     * @param cursor pointer to contingencies in input deck
     * @return vector of measurements
     */
    std::vector<gridpack::state_estimation::Measurement> getMeasurements(
        gridpack::utility::Configuration::ChildCursors measurements);

    /**
     * Read in and partition the network. The input file is read
     * directly from the state_estimation block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a SENetwork object. This should not have any
     * buses or branches defined on it.
     * @param config pointer to open configuration file
     */
    void readNetwork(boost::shared_ptr<SENetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Assume that SENetwork already exists and just cache an internal pointer
     * to it. This routine does not call the partition function. Also read in
     * simulation parameters from configuration file
     * @param network pointer to a complete SENetwork object.
     * @param config pointer to open configuration file
     */
    void setNetwork(boost::shared_ptr<SENetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Read branch and bus measurements. These will come from a separate file.
     * The name of this file comes from the input configuration file.
     */
    void readMeasurements(void);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Solve the state estimation problem
     */
    void solve();

    /**
     * Write final results of state estimation calculation to standard output
     */
    void write();

    /**
     * Save results of state estimation calculation to data collection objects
     */
    void saveData();

    bool hasConverged() const { return p_converged; }

    /**
     * Perform pre-check for suspicious measurements
     */
    void preCheckMeasurements(void);

    /**
     * Perform a targeted check for bad measurements on Bus 8
     */
    void debugMapper();

    /**
     * Add virtual measurements to enforce voltage magnitude constraints
     * @param vmin minimum allowed voltage magnitude
     * @param vmax maximum allowed voltage magnitude
     * @param deviation standard deviation for virtual measurements
     */
    void addVoltageLimitMeasurements(double vmin = 0.9, double vmax = 1.1, 
                                    double deviation = 0.001);

    /**
     * Identify PV buses in the network and their connections
     * Used to properly handle voltage constraints at generator buses
     */
    void identifyPVBusConstraints();

    /**
     * Check for potential measurement inconsistencies
     * Identifies cases where measurements may conflict with physical constraints
     */
    void checkMeasurementConsistency();

    /**
     * Apply proper treatment for PV bus voltage measurements
     * Ensures PV bus voltages are treated as constraints rather than regular measurements
     */
    void handlePVBusVoltages();
    
    /**
     * Apply special handling for voltage angle (VA) measurements
     * Ensures angle measurements at slack buses are treated as constraints
     */
    void handleVAMeasurements();

    // Adjust weights of bad measurements by increasing their sigmas
    void adjustWeights(const std::vector<int>& badIndices);


    private:

    // pointer to network
    boost::shared_ptr<SENetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<SEFactoryModule> p_factory;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<SENetwork> >
      p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<SENetwork> >
      p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;

    // maximum number of iterations
    int p_max_iteration;

    // convergence tolerance
    double p_tolerance;

    // void detectBadData(void);
    std::vector<int> detectBadData(void);
    double p_bad_data_threshold;
    double p_lastChiSquareValue; // Store last chi-square value for reporting

    double getChiSquareThreshold(int dof, double confidence);
    bool p_converged;

    // Helper function to access and modify the sigma of a measurement by index
    double& getMeasurementSigma(int idx);

    // Structure to store information about each bad data iteration
    struct BadDataIterInfo {
        std::vector<int> badIndices;        // Newly identified bad measurements in this iteration
        std::vector<int> allBadIndices;     // All bad measurements identified so far
        double chiSquareValue;
        int iterationNumber;
    };
    
    // Structure to store PV/Slack bus constraint information
    struct PVBusConstraint {
        int busIndex;          // Original bus index
        double voltageValue;   // Fixed voltage value
        bool isPVConnection;   // Whether this bus has connections to non-PV buses
        bool isSlackBus;       // Whether this is a slack bus (vs. PV bus)
    };
    
    // Member variables to store measurement sigmas and bad measurement indices
    std::vector<double> p_measurementSigmas;
    std::vector<int> p_badMeasurementIndices;  // Current bad measurement indices
    std::vector<int> p_allBadMeasurementIndices; // All bad measurement indices across iterations
    std::vector<BadDataIterInfo> p_badDataIterationInfo; // Store information about bad data iterations
    std::vector<PVBusConstraint> p_pvBusConstraints; // Store information about PV bus constraints
};

} // state estimation
} // gridpack
#endif
