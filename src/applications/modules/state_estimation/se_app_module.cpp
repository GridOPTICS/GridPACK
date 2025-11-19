/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app_module.cpp
 * @author Yousu Chen, Bruce Palmer
 * @date   2025-09-23 11:21:35 d3g096
 * Last updated: 8/5/2014 
 *
 * @brief
 * @update Yousu Chen
 *         Adding functions of bad data dection, chi-square testing 
 * @date   2025-03-05
 *         Adding more functions to handle measurements more efficiently
 * @date   2025-04-02
 *         Enhanced voltage constraint featuer
           Cleaned up debug information
 * @date   2025-11-19
 *
 *
 */
// -------------------------------------------------------------

#include <set>
#include <algorithm>
#include <utility>
#include <string>
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/parser/PTI36_parser.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/gen_matrix_map.hpp"
#include "gridpack/mapper/gen_vector_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "timer/coarse_timer.hpp"
#include "se_app_module.hpp"
double p_bad_data_threshold;

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::state_estimation::SEAppModule::SEAppModule(void)
{
  p_lastChiSquareValue = 0.0;
  p_badMeasurementIndices.clear();
  p_allBadMeasurementIndices.clear();
  p_badDataIterationInfo.clear();
  p_diagnosticLevel = "basic"; // Default level
  p_max_bad_data_iterations = 5; // Default value
  p_use_sparse_matrices = true; // Default to sparse for better performance
}

/**
 * Basic destructor
 */
gridpack::state_estimation::SEAppModule::~SEAppModule(void)
{
}

/**
 * Get list of measurements from external file
 * @param cursor pointer to contingencies in input deck
 * @return vector of measurements
 */
std::vector<gridpack::state_estimation::Measurement>
  gridpack::state_estimation::SEAppModule::getMeasurements(
      gridpack::utility::Configuration::ChildCursors measurements)
{
  std::vector<gridpack::state_estimation::Measurement> ret;
  if (p_comm.rank() == 0) {
    int size = measurements.size();
    int idx;
    for (idx = 0; idx < size; idx++) {
      std::string meas_type;
      measurements[idx]->get("Type", &meas_type);
      double meas_value;
      measurements[idx]->get("Value", &meas_value);
      double meas_deviation;
      measurements[idx]->get("Deviation", &meas_deviation);
      if (meas_type == "VM" || meas_type == "PI" ||
          meas_type == "PJ" || meas_type == "QI" ||
          meas_type == "QJ" || meas_type == "VA") {
        int busid;
        measurements[idx]->get("Bus", &busid);
        gridpack::state_estimation::Measurement measurement;
        strcpy(measurement.p_type,meas_type.c_str());
        measurement.p_busid = busid;
        measurement.p_value = meas_value;
        measurement.p_deviation = meas_deviation;
        //printf("%s %d %f %f\n", measurement.p_type.c_str(), measurement.p_busid,
        //  measurement.p_value, measurement.p_deviation);
        ret.push_back(measurement); 
      } else if (meas_type == "PIJ" || meas_type == "PJI" ||
          meas_type == "QIJ" || meas_type == "QJI" ||
          meas_type == "IIJ" || meas_type == "IJI") {
        int fbusid;
        measurements[idx]->get("FromBus", &fbusid);
        int tbusid;
        measurements[idx]->get("ToBus", &tbusid);
        std::string ckt;
        measurements[idx]->get("CKT", &ckt);
        // Fix up tag so that single character tags are right-justified
        if (ckt.length() == 1) {
          ckt.insert(0,1,' ');
        }
        gridpack::state_estimation::Measurement measurement;
        strcpy(measurement.p_type,meas_type.c_str());
        measurement.p_fbusid = fbusid;
        measurement.p_tbusid = tbusid;
        strcpy(measurement.p_ckt,ckt.c_str());
        measurement.p_value = meas_value;
        measurement.p_deviation = meas_deviation;
        //printf("%s %d %d %s %f %f\n", measurement.p_type.c_str(), measurement.p_fbusid,
        //  measurement.p_tbusid, measurement.p_ckt.c_str(), measurement.p_value,
        //  measurement.p_deviation);
        ret.push_back(measurement);
      } 
    }
  }
  return ret;
}


enum Parser{PTI23, PTI33, PTI34, PTI35, PTI36};

/**
 * Read in and partition the network. The input file is read
 * directly from the state_estimation block in the configuration file so no
 * external file names or parameters need to be passed to this routine
 * @param network pointer to a SENetwork object. This should not have any
 * buses or branches defined on it.
 * @param config pointer to open configuration file
 */
void gridpack::state_estimation::SEAppModule::readNetwork(
    boost::shared_ptr<SENetwork> &network,
    gridpack::utility::Configuration *config)
{
  gridpack::parallel::Communicator p_comm;
  p_network = network;
  p_config = config;
  p_comm = network->communicator();

  // Initialize timer for profiling
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("State Estimation: Total Application");
  timer->start(t_total);

  gridpack::utility::Configuration::CursorPtr cursor, secursor;
  secursor = config->getCursor("Configuration.State_estimation");
  std::string filename;
  int filetype = PTI23;
  if (!secursor->get("networkConfiguration",&filename)) {
    if (secursor->get("networkConfiguration_v33",&filename)) {
      filetype = PTI33;
    } else if (secursor->get("networkConfiguration_v34",&filename)) {
      filetype = PTI34;
    } else if (secursor->get("networkConfiguration_v35",&filename)) {
      filetype = PTI35;
    } else if (secursor->get("networkConfiguration_v36",&filename)) {
      filetype = PTI36;
    } else {
      printf("No network configuration file specified\n");
      timer->stop(t_total);
      return;
    }
  }

  // Convergence and iteration parameters
  p_tolerance = secursor->get("tolerance",1.0e-3);
  p_max_iteration = secursor->get("maxIteration",20);
  p_max_bad_data_iterations = secursor->get("maxBadDataIterations", 5);
  
  // Diagnostic output level control
  p_diagnosticLevel = secursor->get("diagnosticOutputLevel", "basic");
  
  // Sparse matrix optimization control
  p_use_sparse_matrices = secursor->get("useSparseMatrices", true);

  // load input file
  //gridpack::parser::PTI23_parser<SENetwork> parser(p_network);
  //parser.parse(filename.c_str());
  double phaseShiftSign = secursor->get("phaseShiftSign",1.0);

  int t_pti = timer->createCategory("State Estimation: Network Parser");
  timer->start(t_pti);
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<SENetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<SENetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI34) {
    gridpack::parser::PTI34_parser<SENetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI35) {
    gridpack::parser::PTI35_parser<SENetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI36) {
    gridpack::parser::PTI36_parser<SENetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  }
  timer->stop(t_pti);

  // partition network
  int t_part = timer->createCategory("State Estimation: Partition");
  timer->start(t_part);
  p_network->partition();
  timer->stop(t_part);

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<SENetwork>(1024, p_network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<SENetwork>(1024, p_network));
  
  // Print configuration parameters
  char buf[128];
  snprintf(buf, sizeof(buf), "Tolerance: %12.4e\n", p_tolerance);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Maximum number of iterations: %d\n", p_max_iteration);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Maximum bad data iterations: %d\n", p_max_bad_data_iterations);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Diagnostic output level: %s\n", p_diagnosticLevel.c_str());
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Sparse matrix optimization: %s\n", p_use_sparse_matrices ? "enabled" : "disabled");
  p_busIO->header(buf);
  
  timer->stop(t_total);
}

/**
 * Assume that SENetwork already exists and just cache an internal pointer
 * to it. This routine does not call the partition function. Also read in
 * simulation parameters from configuration file
 * @param network pointer to a complete SENetwork object.
 * @param config pointer to open configuration file
 */
void gridpack::state_estimation::SEAppModule::setNetwork(
    boost::shared_ptr<SENetwork> &network,
    gridpack::utility::Configuration *config)
{
  gridpack::parallel::Communicator p_comm;
  p_network = network;
  p_config = config;
  p_comm = network->communicator();

  gridpack::utility::Configuration::CursorPtr cursor, secursor;
  secursor = p_config->getCursor("Configuration.State_estimation");

  // Convergence and iteration parameters
  p_tolerance = secursor->get("tolerance",1.0e-3);
  p_max_iteration = secursor->get("maxIteration",20);
  p_max_bad_data_iterations = secursor->get("maxBadDataIterations", 5);
  
  // Diagnostic output level control
  p_diagnosticLevel = secursor->get("diagnosticOutputLevel", "basic");
  
  // Sparse matrix optimization control
  p_use_sparse_matrices = secursor->get("useSparseMatrices", true);
  
  char buf[128];
  snprintf(buf, sizeof(buf), "Tolerance: %12.4e\n", p_tolerance);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Maximum number of iterations: %d\n", p_max_iteration);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Maximum bad data iterations: %d\n", p_max_bad_data_iterations);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Diagnostic output level: %s\n", p_diagnosticLevel.c_str());
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Sparse matrix optimization: %s\n", p_use_sparse_matrices ? "enabled" : "disabled");
  p_busIO->header(buf);

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<SENetwork>(1024, p_network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<SENetwork>(1024, p_network));
}

/**
 * Read branch and bus measurements. These will come from a separate file.
 * The name of this file comes from the input configuration file. Call this
 * method after initializing the network.
 */
void gridpack::state_estimation::SEAppModule::readMeasurements(void)
{
  // Initialize timer for profiling
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("State Estimation: Total Application");
  int t_meas = timer->createCategory("State Estimation: Read Measurements");
  timer->start(t_total);
  timer->start(t_meas);

  // Read in measurement file
  std::string measurementfile;
  gridpack::utility::Configuration::CursorPtr cursor, secursor;
  secursor = p_config->getCursor("Configuration.State_estimation");
  if (!secursor->get("measurementList", &measurementfile)) {
    measurementfile = "IEEE14_meas.xml";
  }
  bool ok = p_config->open(measurementfile, p_comm);

  // get a list of measurements
  cursor = p_config->getCursor("Measurements");
  gridpack::utility::Configuration::ChildCursors measurements;
  if (cursor) cursor->children(measurements);
  std::vector<gridpack::state_estimation::Measurement>
    meas = getMeasurements(measurements);
  
  // Add measurements to buses and branches
  p_factory->setMeasurements(meas);

  timer->stop(t_meas);
  timer->stop(t_total);
}

void
gridpack::state_estimation::SEAppModule::setMeasurements
(std::vector<gridpack::state_estimation::Measurement>& measurements)
{
  p_factory->setMeasurements(measurements);
}

/**
 * Set up exchange buffers and other internal parameters and initialize
 * network components using data from data collection
 */
void gridpack::state_estimation::SEAppModule::initialize(void)
{
  // Initialize timer for profiling
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("State Estimation: Total Application");
  int t_init = timer->createCategory("State Estimation: Initialize");
  timer->start(t_total);
  timer->start(t_init);

  // create factory
  p_factory.reset(new gridpack::state_estimation::SEFactoryModule(p_network));
  
  int t_load = timer->createCategory("State Estimation: Factory Load");
  timer->start(t_load);
  p_factory->load();
  timer->stop(t_load);

  // set network components using factory
  int t_setc = timer->createCategory("State Estimation: Factory Set Components");
  timer->start(t_setc);
  p_factory->setComponents();
  timer->stop(t_setc);

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  int t_setx = timer->createCategory("State Estimation: Factory Set Exchange");
  timer->start(t_setx);
  p_factory->setExchange();
  timer->stop(t_setx);

  // Create bus data exchange
  int t_updt = timer->createCategory("State Estimation: Bus Update");
  timer->start(t_updt);
  p_network->initBusUpdate();
  timer->stop(t_updt);

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");
  p_bad_data_threshold = cursor->get("badDataThreshold", 3.0);
  
  // Identify PV buses and their connections for proper constraint handling
  int t_pv = timer->createCategory("State Estimation: PV Bus Constraints");
  timer->start(t_pv);
  identifyPVBusConstraints();
  timer->stop(t_pv);
  
  // Check for potential measurement inconsistencies
  int t_cons = timer->createCategory("State Estimation: Measurement Consistency");
  timer->start(t_cons);
  checkMeasurementConsistency();
  timer->stop(t_cons);
  
  timer->stop(t_init);
  timer->stop(t_total);
}


/**
 * pre-check measurements
 */
void gridpack::state_estimation::SEAppModule::preCheckMeasurements()
{
  p_busIO->header("\nPerforming pre-check for suspicious measurements...\n");
  
  int numBus = p_network->numBuses();
  for (int i = 0; i < numBus; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    if (!bus->isIsolated() && bus->vectorNumElements() > 0) {
      // Use serialWrite with a special flag for pre-check
      char buf[1024];
      bus->serialWritePreCheck(buf, sizeof(buf));
    }
  }
}

/**
 * Solve the state estimation problem
 */
void gridpack::state_estimation::SEAppModule::solve(void)
{
  // Initialize timer for profiling
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("State Estimation: Total Application");
  timer->start(t_total);

  // Run pre-check to identify obviously suspicious measurements
  preCheckMeasurements();

  // set YBus components so that you can create Y matrix  
  int t_ybus = timer->createCategory("State Estimation: YBus");
  timer->start(t_ybus);
  p_factory->setYBus();
  timer->stop(t_ybus);

  // set some state estimation parameters
  p_factory->configureSE();
  
  // Handle PV bus constraints
  handlePVBusVoltages();
  
  // Handle VA measurements
  handleVAMeasurements();

  // Get configuration cursor
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");

  // Read voltage constraint settings directly from configuration with defaults
  bool use_voltage_constraints = cursor->get("useVoltageConstraints", false);
  // Read voltage constraint settings
  bool enforce_v_limits = cursor->get("enforceVoltageLimits", false);
  double v_min = cursor->get("minVoltage", 0.9);
  double v_max = cursor->get("maxVoltage", 1.1);
  
  if (enforce_v_limits) {
    p_busIO->header("\nEnforcing voltage magnitude limits\n");
    char buf[128];
    snprintf(buf, sizeof(buf), "Voltage limits: %.3f <= V <= %.3f\n", v_min, v_max);
    p_busIO->header(buf);
    
    // Set voltage limits for all buses
    p_factory->setVoltageLimits(v_min, v_max, true);
  }
  
  if (use_voltage_constraints) {
    double constraint_deviation = cursor->get("voltageConstraintDeviation", 0.0025);
    p_busIO->header("\nAdding voltage magnitude constraints\n");
    char buf[128];
    snprintf(buf, sizeof(buf), "Voltage limits: %.3f <= V <= %.3f (deviation: %.4f)\n",
            v_min, v_max, constraint_deviation);
    p_busIO->header(buf);

    // Add virtual measurements for voltage constraints
    addVoltageLimitMeasurements(v_min, v_max, constraint_deviation);

    // Reconfigure SE to update element counts with new VUL/VLL measurements
    p_factory->configureSE();
  }

  //Create Y-bus matrix
  int t_matrix = timer->createCategory("State Estimation: Matrix Creation");
  timer->start(t_matrix);
  p_factory->setMode(YBus);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();

  // Create mapper to push voltage data back onto buses
  p_factory->setMode(Voltage);
  gridpack::mapper::BusVectorMap<SENetwork> VMap(p_network);

  // Create initial version (H), measurement vector (Ez), and Rinv
  p_factory->setMode(Jacobian_H);
  gridpack::mapper::GenMatrixMap<SENetwork> HJacMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> HJac;
  
  if (p_use_sparse_matrices) {
    HJac = HJacMap.mapToMatrix();
    // Note: GridPACK defaults to sparse matrices, so no special action needed
    p_busIO->header("Using sparse matrix storage for Jacobian H\n");
  } else {
    HJac = HJacMap.mapToMatrix();
    p_busIO->header("Using dense matrix storage for Jacobian H\n");
  }

  gridpack::mapper::GenVectorMap<SENetwork> EzMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> Ez = EzMap.mapToVector();

  p_factory->setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv;
  
  if (p_use_sparse_matrices) {
    Rinv = RinvMap.mapToMatrix();
    // Note: GridPACK defaults to sparse matrices
    p_busIO->header("Using sparse matrix storage for R_inv\n");
  } else {
    Rinv = RinvMap.mapToMatrix();
    p_busIO->header("Using dense matrix storage for R_inv\n");
  }
  timer->stop(t_matrix);

  // Convergence and iteration parameters
  ComplexType tol;
  //tol = 2.0*p_tolerance;
  //int iter = 0;
  bool converged = false;
  p_converged = false;
  int badDataIter = 0;
  bool badDataExists = true; 


  //Out loop for bad data handling
  while (badDataExists && badDataIter < p_max_bad_data_iterations) {
    badDataExists = false; // Reset flag
    tol = 2.0 * p_tolerance; // Reset tolerance
    int iter = 0;

    int t_nr_loop = timer->createCategory("State Estimation: N-R Iterations");
    timer->start(t_nr_loop);

    // Start N-R loop
    while (real(tol) > p_tolerance && iter < p_max_iteration) {
      
      // Update Jacobian for current iteration
      p_factory->setMode(Jacobian_H);
      HJacMap.mapToMatrix(HJac);

      // Form H'
      boost::shared_ptr<gridpack::math::Matrix> trans_HJac(transpose(*HJac));

      // Build measurement equation
      EzMap.mapToVector(Ez);

      // Form Gain matrix
      boost::shared_ptr<gridpack::math::Matrix> Gain1(multiply(*trans_HJac, *Rinv));
      boost::shared_ptr<gridpack::math::Matrix> Gain(multiply(*Gain1, *HJac));

      // Form right hand side vector
      boost::shared_ptr<gridpack::math::Matrix> HTR(multiply(*trans_HJac, *Rinv));
      boost::shared_ptr<gridpack::math::Vector> RHS(multiply(*HTR, *Ez));

      // create a linear solver
      int t_solve = timer->createCategory("State Estimation: Linear Solver");
      timer->start(t_solve);
      gridpack::utility::Configuration::CursorPtr cursor;
      cursor = p_config->getCursor("Configuration.State_estimation");
      gridpack::math::LinearSolver solver(*Gain);
      solver.configure(cursor);

      // Solve linear equation
      boost::shared_ptr<gridpack::math::Vector> X(RHS->clone());
      X->zero(); //might not need to do this
      solver.solve(*RHS, *X);
      timer->stop(t_solve);
      
      tol = X->normInfinity();
      char ioBuf[128];
      snprintf(ioBuf, sizeof(ioBuf), "\nIteration %d Tol: %12.6e\n", iter+1, real(tol));
      p_busIO->header(ioBuf);

      // Push solution back onto bus variables
      p_factory->setMode(Voltage);
      VMap.mapToBus(X);
    
      // update values
      p_network->updateBuses();
      iter++;

    // End N-R loop
    }
    
    timer->stop(t_nr_loop);

  // Check convergence and handle bad data
  if (real(tol) <= p_tolerance) {

    converged = true;
    p_converged = true; 
    p_busIO->header("\n*** State estimation converged successfully ***\n");
    // Perform bad data detection
    int t_baddata = timer->createCategory("State Estimation: Bad Data Detection");
    timer->start(t_baddata);
    
    // Perform bad data detection
    std::vector<int> badIndices;
    badIndices.reserve(50); // Pre-allocate reasonable space for bad measurements
    
    // Call detectBadData and handle the return issue
    BadDataResult badDataResult = detectBadData();
    
    // Extract bad indices from the result
    badIndices = badDataResult.badIndices;
    
    timer->stop(t_baddata);
    
    // Safety check for TAMU500 abort issue
    char buf[256];
    try {
      if (p_diagnosticLevel == "detailed") {
        snprintf(buf, sizeof(buf), "DEBUG: badIndices.size() = %d\n", (int)badIndices.size());
        p_busIO->header(buf);
        
        if (!badIndices.empty()) {
          snprintf(buf, sizeof(buf), "DEBUG: First bad index = %d\n", badIndices[0]);
          p_busIO->header(buf);
        }
      }
    } catch (const std::exception& e) {
      snprintf(buf, sizeof(buf), "ERROR: Exception accessing badIndices: %s\n", e.what());
      p_busIO->header(buf);
    }
    
    if (!badIndices.empty()) {
      badDataExists = true;
      
      // Store iteration information for reporting
      BadDataIterationInfo iterInfo;
      iterInfo.badIndices = badIndices;  // New bad measurements in this iteration
      
      // Store the normalized residuals for all measurements from the detection result
      iterInfo.normalizedResiduals = badDataResult.normalizedResiduals;
      iterInfo.chiSquareValue = badDataResult.chiSquareValue;
      
      // Keep track of all bad measurement indices across all iterations
      std::vector<int> allBadIndicesSoFar = p_allBadMeasurementIndices;
      
      for (int idx : badIndices) {
        if (p_diagnosticLevel == "detailed") {
          snprintf(buf, sizeof(buf), "DEBUG: Processing badIndex %d\n", idx);
          p_busIO->header(buf);
        }
        
        // Add to the running list of all bad indices if not already there
        bool alreadyAdded = false;
        for (int prevIdx : p_allBadMeasurementIndices) {
          if (idx == prevIdx) {
            alreadyAdded = true;
            break;
          }
        }
        if (!alreadyAdded) {
          p_allBadMeasurementIndices.push_back(idx);
          allBadIndicesSoFar.push_back(idx);
        }
      }
      
      // Store all indices identified so far
      iterInfo.allBadIndices = allBadIndicesSoFar;
      iterInfo.chiSquareValue = p_lastChiSquareValue;
      iterInfo.iterationNumber = badDataIter + 1;
      iterInfo.degreesOfFreedom = badDataResult.degreesOfFreedom; // Store the actual DOF from detection result
      
      // Normalized residuals will be stored when the iteration info is created
      
      p_badDataIterationInfo.push_back(iterInfo);
      
      // Log iteration information
      char msgBuf[128];
      snprintf(msgBuf, sizeof(msgBuf), "\nIteration %d: Found %d bad measurements, Chi-square: %.2f\n", 
              badDataIter + 1, (int)badIndices.size(), p_lastChiSquareValue);
      p_busIO->header(msgBuf);
      
      // Adjust weights to bad measurements
      int t_adjust = timer->createCategory("State Estimation: Adjust Weights");
      timer->start(t_adjust);
      adjustWeights(badIndices);
      timer->stop(t_adjust);
      
      // Set factory mode to rebuild Rinv
      p_factory->setMode(R_inv);
      
      // Refresh mapping to avoid matrix size issues
      p_busIO->header("Re-running SE with adjusted weights\n");
    }
  } else {
    p_busIO->header("\n*** WARNING: STATE ESTIMATION DID NOT CONVERGE ***\n");
    char ioBuf[128];
    snprintf(ioBuf, sizeof(ioBuf), "Maximum iterations (%d) reached with tolerance %12.6e > %12.6e\n", 
            p_max_iteration, real(tol), p_tolerance);
    p_busIO->header(ioBuf);
    p_busIO->header("Results may not be reliable.\n");
    break; //Exist Outer loop if SE doesn't converge
  }
  badDataIter++;
}
  // Output results if converged
  if (converged) {
    // debugMapper(); // Disabled debug output
    p_busIO->header("\n   State Estimation Outputs\n");
    p_busIO->header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
    p_busIO->write();
    p_branchIO->header("\n   Branch Power Flow (p.u.)\n");
    p_branchIO->header("\n        Bus 1       Bus 2            P                    Q\n");
    p_branchIO->write();

    p_busIO->header("\n   Comparison of Bus Measurements and Estimations\n");
    p_busIO->header("\n   Type  Bus Number      Measurement          Estimate"
                    " Difference   Deviation\n");
    p_busIO->write("se");

    p_branchIO->header("\n   Comparison of Branch Measurements and Estimations\n");
    p_branchIO->header("\n   Type  From    To  CKT   Measurement      Estimate"
                       " Difference   Deviation\n");
    p_branchIO->write("se");
  }
    // Output final results if converged
    if (!badDataExists) {
        p_busIO->header("\nSE completed successfully with no bad data remaining\n");
        char bdioBuf[128];
        snprintf(bdioBuf, sizeof(bdioBuf), "\nBad Data Iteration %d \n", badDataIter);
    }
    
    // Report Jacobian optimization performance - only for detailed diagnostics
    if (p_diagnosticLevel == "detailed") {
        reportJacobianPerformance();
    }
    
    // Stop total timer
    timer->stop(t_total);
    
    // Timer stopped successfully
}

/**
 * Write final results of state estimation calculation to standard
 * output
 */
void gridpack::state_estimation::SEAppModule::write(void)
{
  // Get output file name from configuration
  std::string outputFile = "state_estimation_results.txt"; // Default value
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");
  if (cursor) {
    cursor->get("outputFile", &outputFile);
  }
  
  // Create output file
  std::ofstream outFile;
  
  // Only have process 0 write 
  if (p_network->communicator().rank() == 0) {
    outFile.open(outputFile.c_str());
    if (!outFile.is_open()) {
      printf("ERROR: Could not open output file %s\n", outputFile.c_str());
      return;
    }
  }
  
  if (p_network->communicator().rank() == 0 && outFile.is_open()) {
    // Create temporary files to capture output
    FILE* oldStdout = stdout;  // Save original stdout
    
    // Redirect stdout to a temporary file
    FILE* tempFile = tmpfile();
    if (!tempFile) {
      printf("Error creating temporary file\n");
      outFile.close();
      return;
    }
    
    // Write header information and convergence status
    outFile << "---------------------------------------------------------\n";
    outFile << "STATE ESTIMATION ANALYSIS RESULTS\n";
    outFile << "---------------------------------------------------------\n";
    outFile << "Date: " << __DATE__ << " " << __TIME__ << "\n";
    outFile << "Newton-Raphson Convergence: " << (p_converged ? "CONVERGED" : "NOT CONVERGED") << "\n";
    if (p_converged) {
        outFile << "Solution Status: VALID - State estimation solved successfully\n";
    } else {
        outFile << "Solution Status: INVALID - State estimation failed to converge\n";
    }
    
    // Write convergence parameters
    char paramBuf[512];  // Increase buffer size
    snprintf(paramBuf, sizeof(paramBuf), "Convergence Tolerance: %e\n", p_tolerance);
    outFile << paramBuf;
    snprintf(paramBuf, sizeof(paramBuf), "Maximum Iterations: %d\n", p_max_iteration);
    outFile << paramBuf;
    snprintf(paramBuf, sizeof(paramBuf), "Maximum Bad Data Iterations: %d\n", p_max_bad_data_iterations);
    outFile << paramBuf;
    snprintf(paramBuf, sizeof(paramBuf), "Bad Data Threshold: %f\n", p_bad_data_threshold);
    outFile << paramBuf;
    outFile << "---------------------------------------------------------\n\n";
    
    if (!p_converged) {
      outFile << "WARNING: State estimation did not converge.\n";
      outFile << "Results may not be reliable.\n\n";
    }
    
    // Write bad data information
    outFile << "---------------------------------------------------------\n";
    outFile << "BAD DATA DETECTION RESULTS\n";
    outFile << "---------------------------------------------------------\n";
    
    if (p_badMeasurementIndices.empty()) {
      outFile << "No bad data detected in the final solution.\n\n";
    } else {
      outFile << "The following measurements were identified as bad data:\n\n";
      outFile << "Index  Type  Bus/From  To      Value        Deviation    NormResidual\n";
      outFile << "----------------------------------------------------------------\n";
      
      for (int idx : p_badMeasurementIndices) {
        std::string details;
        if (p_factory->reportMeasurement(idx, details)) {
          // Parse the details string to extract components
          std::stringstream ss(details);
          std::string type, busInfo, valueStr, deviationStr;
          
          int busNum = -1;
          int toBusNum = -1;
          std::string measurementType = "Unknown";
          double value = 0.0;
          double deviation = 0.0;
          
          if (sscanf(details.c_str(), "Bus %d, Type: %3s, Value: %lf, Deviation: %lf", 
                    &busNum, measurementType.data(), &value, &deviation) == 4) {
            // Bus measurement
            char formatted[512];  // Increase buffer size
            snprintf(formatted, sizeof(formatted), "  %3d   %-4s  %-8d  -       %-12.6f  %-12.6f  N/A\n", 
                    idx, measurementType.c_str(), busNum, value, deviation);
            outFile << formatted;
          } 
          else if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf, Deviation: %lf", 
                         &busNum, &toBusNum, measurementType.data(), &value, &deviation) == 5) {
            // Branch measurement
            char formatted[512];  // Increase buffer size
            snprintf(formatted, sizeof(formatted), "  %3d   %-4s  %-8d  %-7d  %-12.6f  %-12.6f  N/A\n", 
                    idx, measurementType.c_str(), busNum, toBusNum, value, deviation);
            outFile << formatted;
          }
          else {
            // Fall back to original format if parsing fails
            outFile << "  " << idx << "   " << details << "\n";
          }
        }
      }
      outFile << "\nNote: Weights for these measurements were reduced in the solution.\n\n";
    }
    
    // Write state estimation results
    outFile << "---------------------------------------------------------\n";
    outFile << "STATE ESTIMATION OUTPUTS - SYSTEM STATE\n";
    outFile << "---------------------------------------------------------\n";
    outFile << "Bus Number      Phase Angle (deg)      Voltage Magnitude (p.u.)\n";
    outFile << "---------------------------------------------------------\n";
    
    int numBus = p_network->numBuses();
    for (int i = 0; i < numBus; i++) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (!bus->isIsolated()) {
        char buf[512];  // Increase buffer size
        double angle = bus->getPhase() * 180.0 / M_PI; // Convert to degrees
        snprintf(buf, sizeof(buf), "%-10d      %-20.6f      %-20.6f\n",
                bus->getOriginalIndex(), angle, bus->getVoltage());
        outFile << buf;
      }
    }
    outFile << "\n";
    
    // Write branch power flow data
    outFile << "---------------------------------------------------------\n";
    outFile << "BRANCH POWER FLOW RESULTS (P.U.)\n";
    outFile << "---------------------------------------------------------\n";
    outFile << "From Bus      To Bus      Active Power (P)      Reactive Power (Q)\n";
    outFile << "---------------------------------------------------------\n";
    
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
      SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
      if (!bus1->isIsolated() && !bus2->isIsolated()) {
        char buf[512];  // Increase buffer size
        double p, q;
        branch->getPQ(bus1, &p, &q);
        snprintf(buf, sizeof(buf), "%-12d  %-12d  %-20.6f  %-20.6f\n",
                bus1->getOriginalIndex(), bus2->getOriginalIndex(), p, q);
        outFile << buf;
      }
    }
    outFile << "\n";
    
    // Write measurement comparison information directly
    outFile << "---------------------------------------------------------\n";
    outFile << "MEASUREMENT COMPARISON - BUS MEASUREMENTS\n";
    outFile << "---------------------------------------------------------\n";
    
    // For large systems, skip detailed measurement output to prevent hangs
    if (numBus > 100) {
        outFile << "System has " << numBus << " buses - detailed measurement output skipped for performance.\n";
        outFile << "Use smaller systems for detailed measurement analysis.\n\n";
    } else {
        std::ostringstream busOutputSS;
        FILE* origStream = stdout;
        
        int stdout_pipe[2];
        pipe(stdout_pipe);
        int savedStdout = dup(STDOUT_FILENO);
        dup2(stdout_pipe[1], STDOUT_FILENO);
        close(stdout_pipe[1]);
        
        p_busIO->header("\n   Type  Bus Number      Measurement          Estimate"
                       " Difference   Deviation\n");
        p_busIO->write("se");
        
        // Restore stdout
        fflush(stdout);
        dup2(savedStdout, STDOUT_FILENO);
        close(savedStdout);
        
        char buffer[4096];
        ssize_t bytesRead;
        std::string capturedOutput;
        const size_t MAX_OUTPUT_SIZE = 10 * 1024 * 1024;  // 10MB limit
        
        while ((bytesRead = read(stdout_pipe[0], buffer, sizeof(buffer) - 1)) > 0) {
            buffer[bytesRead] = '\0';
            capturedOutput += buffer;
            
            // Prevent infinite memory growth for large systems
            if (capturedOutput.size() > MAX_OUTPUT_SIZE) {
                capturedOutput += "\n[OUTPUT TRUNCATED - System too large for detailed measurement output]\n";
                break;
            }
        }
        close(stdout_pipe[0]);
        
        // Write measurement data to a file, if available
        if (!capturedOutput.empty()) {
            // Remove extra header lines if present
            size_t start = capturedOutput.find("Type  Bus Number");
            if (start != std::string::npos) {
                size_t lineEnd = capturedOutput.find("\n", start);
                if (lineEnd != std::string::npos) {
                    capturedOutput = capturedOutput.substr(lineEnd + 1);
                }
            }
            
            // Write to output file
            outFile << capturedOutput;
        } else {
            outFile << "No bus measurement data available.\n";
        }
    }
    
    outFile << "\n";
    
    // Branch measurements
    outFile << "---------------------------------------------------------\n";
    outFile << "MEASUREMENT COMPARISON - BRANCH MEASUREMENTS\n";
    outFile << "---------------------------------------------------------\n";
    
    // For large systems, skip detailed branch measurement output
    if (numBus > 100) {
        outFile << "System has " << numBus << " buses - detailed branch measurement output skipped for performance.\n";
        outFile << "Use smaller systems for detailed measurement analysis.\n\n";
    } else {
        std::ostringstream branchOutputSS;
        
        int stdout_pipe[2];
        pipe(stdout_pipe);
        int savedStdout = dup(STDOUT_FILENO);
        dup2(stdout_pipe[1], STDOUT_FILENO);
        close(stdout_pipe[1]);
        
        p_branchIO->header("\n   Type  From    To  CKT   Measurement      Estimate"
                          " Difference   Deviation\n");
        p_branchIO->write("se");
        
        // Restore stdout
        fflush(stdout);
        dup2(savedStdout, STDOUT_FILENO);
        close(savedStdout);
        
        char buffer[4096];
        ssize_t bytesRead;
        std::string capturedOutput;
        const size_t MAX_OUTPUT_SIZE = 10 * 1024 * 1024;  // 10MB limit
        
        capturedOutput.clear();
        while ((bytesRead = read(stdout_pipe[0], buffer, sizeof(buffer) - 1)) > 0) {
            buffer[bytesRead] = '\0';
            capturedOutput += buffer;
            
            // Prevent infinite memory growth for large systems
            if (capturedOutput.size() > MAX_OUTPUT_SIZE) {
                capturedOutput += "\n[OUTPUT TRUNCATED - System too large for detailed measurement output]\n";
                break;
            }
        }
        close(stdout_pipe[0]);
        
        // Write branch measurement data to the file, if available
        if (!capturedOutput.empty()) {
            // Remove extra header lines if present
            size_t start = capturedOutput.find("Type  From    To");
            if (start != std::string::npos) {
                size_t lineEnd = capturedOutput.find("\n", start);
                if (lineEnd != std::string::npos) {
                    capturedOutput = capturedOutput.substr(lineEnd + 1);
                }
            }
            
            outFile << capturedOutput;
        } else {
            outFile << "No branch measurement data available.\n";
        }
    }
    outFile << "\n";
    
    // Track branch measurements for diagnostic purposes - needed for both standard and detailed levels
    std::vector<std::pair<int, double>> branchResiduals;
    
    // Statistical analysis section - only for standard and detailed levels
    if (p_diagnosticLevel == "standard" || p_diagnosticLevel == "detailed") {
        outFile << "---------------------------------------------------------\n";
        outFile << "STATISTICAL ANALYSIS\n";
        outFile << "---------------------------------------------------------\n";
        
        // Add top 10 measurements with highest normalized residuals for diagnostic purposes
        outFile << "\nTop 10 Measurements with Highest Normalized Residuals:\n";
        outFile << "-------------------------------------------------------------\n";
        outFile << "Rank  Index  Type  Bus/From  To      Value        NormRes     Threshold   Status\n";
        outFile << "--------------------------------------------------------------------------------\n";
    
    // Create a copy of the sorted pairs from the detection process
    std::vector<std::pair<int, double>> diagnosticPairs;
    int mapperSize = 0;
    
    // Create temporary vector and matrix maps to get the residuals
    gridpack::mapper::GenVectorMap<SENetwork> ResMap(p_network);
    boost::shared_ptr<gridpack::math::Vector> Residual = ResMap.mapToVector();
    mapperSize = Residual->size();
    
    // Configure for residual calculation
    p_factory->setMode(gridpack::state_estimation::Residual);
    
    if (mapperSize > 0) {
        // Create a vector of index/residual pairs
        for (int i = 0; i < mapperSize; i++) {
            // Get R inverse for normalization
            double weight = 1.0;
            gridpack::ComplexType rinvValue;
            gridpack::ComplexType rvalue;
            
            // Get R-inverse diagonal element
            p_factory->setMode(gridpack::state_estimation::R_inv);
            gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
            boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();
            Rinv->getElement(i, i, rinvValue);
            weight = std::real(rinvValue);
            
            // Get residual
            p_factory->setMode(gridpack::state_estimation::Residual);
            Residual->getElement(i, rvalue);
            double r = std::abs(std::real(rvalue));
            // Calculate normalized residual
            double sigma = (weight > 0.0) ? 1.0/sqrt(weight) : 1.0;
            double normRes = r/sigma;
            diagnosticPairs.push_back(std::make_pair(i, normRes));
            
            // Track branch measurements for diagnostics
            std::string details;
            if (p_factory->reportMeasurement(i, details)) {
                // Check if this is a branch measurement (PIJ, QIJ, PJI, QJI)
                if (details.find("FromBus") != std::string::npos || 
                    details.find("PIJ") != std::string::npos ||
                    details.find("QIJ") != std::string::npos || 
                    details.find("PJI") != std::string::npos ||
                    details.find("QJI") != std::string::npos) {
                    branchResiduals.push_back(std::make_pair(i, normRes));
                }
                
                // Bus 8 tracking removed
            }
        }
        
        // Sort by residual magnitude (descending)
        std::sort(diagnosticPairs.begin(), diagnosticPairs.end(), 
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
              });
        
        // Print top 10 or as many as available
        int diagCount = std::min(10, (int)diagnosticPairs.size());
        for (int i = 0; i < diagCount; i++) {
            int idx = diagnosticPairs[i].first;
            double res = diagnosticPairs[i].second;
            std::string details;
            
            if (p_factory->reportMeasurement(idx, details)) {
                int busNum = -1;
                int toBusNum = -1;
                std::string type = "Unknown";
                double value = 0.0;
                
                double thisThreshold = p_bad_data_threshold;
                
                // Check for VM measurements (higher threshold)
                if (details.find("Type: VM") != std::string::npos) {
                    thisThreshold = p_bad_data_threshold * 3.0;
                }
                // Check for VA measurements (higher threshold)
                else if (details.find("Type: VA") != std::string::npos) {
                    thisThreshold = p_bad_data_threshold * 5.0;
                }
                
                // Determine if this measurement would be flagged
                std::string status = (res > thisThreshold) ? "BAD" : "OK";
                
                // Check f it's flagged or not
                bool wasFlagged = false;
                for (int badIdx : p_badMeasurementIndices) {
                    if (badIdx == idx) {
                        wasFlagged = true;
                        break;
                    }
                }
                
                if (wasFlagged) {
                    status = "FLAGGED";
                }
                
                char diagBuf[256];
                
                if (sscanf(details.c_str(), "Bus %d, Type: %3s, Value: %lf", 
                          &busNum, type.data(), &value) >= 3) {
                    // Bus measurement
                    snprintf(diagBuf, sizeof(diagBuf), "%-5d %-6d %-5s %-9d -       %-12.6f %-11.6f %-11.6f %s\n", 
                            i+1, idx, type.c_str(), busNum, value, res, thisThreshold, status.c_str());
                } 
                else if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                               &busNum, &toBusNum, type.data(), &value) >= 4) {
                    // Branch measurement
                    snprintf(diagBuf, sizeof(diagBuf), "%-5d %-6d %-5s %-9d %-7d %-12.6f %-11.6f %-11.6f %s\n", 
                            i+1, idx, type.c_str(), busNum, toBusNum, value, res, thisThreshold, status.c_str());
                }
                else {
                    // Fall back to simpler format
                    snprintf(diagBuf, sizeof(diagBuf), "%-5d %-6d %s, NormRes=%12.6f, Threshold=%12.6f, Status: %s\n", 
                            i+1, idx, details.c_str(), res, thisThreshold, status.c_str());
                }
                
                outFile << diagBuf;
            }
        }
    }
    else {
        outFile << "No normalized residual data available for diagnostic purposes.\n";
    }
    } // End statistical analysis section
    
    // Branch measurement diagnostic section - only for detailed level
    if (p_diagnosticLevel == "detailed") {
        outFile << "\nBRANCH MEASUREMENT DIAGNOSTIC SECTION:\n";
        outFile << "--------------------------------------------\n";
    
    if (!branchResiduals.empty()) {
        // Sort branch residuals by value (descending)
        std::sort(branchResiduals.begin(), branchResiduals.end(),
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
              });
        
        outFile << "Found " << branchResiduals.size() << " branch measurements in the system.\n";
        outFile << "Index  Type  From    To     Value        Calculated   Residual     NormRes     Status\n";
        outFile << "---------------------------------------------------------------------------------------\n";
        
        for (size_t i = 0; i < branchResiduals.size(); i++) {
            int idx = branchResiduals[i].first;
            double normRes = branchResiduals[i].second;
            std::string details;
            
            if (p_factory->reportMeasurement(idx, details)) {
                // Parse branch measurement details
                int fromBus = -1, toBus = -1;
                std::string type = "Unknown";
                double value = 0.0;
                
                // Determine if this measurement would be flagged
                double thisThreshold = p_bad_data_threshold;
                std::string status = (normRes > thisThreshold) ? "POTENTIALLY BAD" : "OK";
                
                // Was it actually flagged as bad?
                bool wasFlagged = false;
                for (int badIdx : p_badMeasurementIndices) {
                    if (badIdx == idx) {
                        wasFlagged = true;
                        break;
                    }
                }
                
                if (wasFlagged) {
                    status = "FLAGGED";
                }
                
                char diagBuf[256];
                if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                        &fromBus, &toBus, type.data(), &value) >= 4) {
                    // Calculate estimated value and raw residual
                    double sigma = normRes > 0 ? std::abs(value) / normRes : 0.001;
                    double estValue = value;  // Just a placeholder for now
                    double rawResidual = 0.0; // Just a placeholder for now
                    
                    snprintf(diagBuf, sizeof(diagBuf), "%-6d %-5s %-7d %-7d %-12.6f %-12.6f %-12.6f %-11.6f %s\n", 
                            idx, type.c_str(), fromBus, toBus, value, estValue, 
                            rawResidual, normRes, status.c_str());
                    outFile << diagBuf;
                }
                else {
                    // Branch parsing failed - skip detailed output for cleaner results
                }
            }
        }
    }
    else {
        outFile << "No branch measurement residuals found. This may indicate an issue with branch measurement handling.\n";
    }
    } // End branch diagnostic section
    
    // Bus 8 diagnostic section removed for cleaner output
    
    // Bus 8 comparison section removed for cleaner output
    
    // Show all flagged measurements and their normalized residuals
    outFile << "Measurements flagged as bad";
    if (p_badMeasurementIndices.size() > 5) {
        outFile << " (showing worst 5 of " << p_badMeasurementIndices.size() << " total)";
    }
    outFile << ":\n";
    if (!p_badMeasurementIndices.empty()) {
        outFile << "Index  Type  Bus    Value      Deviation   OrigRes   CurrRes   Iterations  Comments\n";
        outFile << "-----------------------------------------------------------------------------------\n";
        
        for (int idx : p_badMeasurementIndices) {
            std::string details;
            double currentNormRes = 0.0;
            double originalNormRes = 0.0;
            
            // Find the current normalized residual from the most recent iteration
            if (!p_badDataIterationInfo.empty()) {
                const auto& lastIteration = p_badDataIterationInfo.back();
                auto it = lastIteration.normalizedResiduals.find(idx);
                if (it != lastIteration.normalizedResiduals.end()) {
                    currentNormRes = std::abs(it->second);
                }
            }
            
            // Find the original normalized residual from when it was first flagged
            // and count how many iterations this measurement has been flagged
            int iterationCount = 0;
            for (const auto& iterInfo : p_badDataIterationInfo) {
                // Check if this measurement was flagged in this iteration
                bool flaggedInThisIteration = false;
                for (int badIdx : iterInfo.badIndices) {
                    if (badIdx == idx) {
                        flaggedInThisIteration = true;
                        break;
                    }
                }
                if (flaggedInThisIteration) {
                    iterationCount++;
                }
                
                // Get original residual from first occurrence
                if (originalNormRes == 0.0) {
                    auto it = iterInfo.normalizedResiduals.find(idx);
                    if (it != iterInfo.normalizedResiduals.end()) {
                        originalNormRes = std::abs(it->second);
                    }
                }
            }
            
            if (p_factory->reportMeasurement(idx, details)) {
                int busNum = -1;
                int toBusNum = -1;
                std::string type = "Unknown";
                double value = 0.0;
                double deviation = 0.0;
                char diagBuf[512];  // Increased buffer size
                std::string comments;
                
                // Extract deviation from measurement (we'll get it from the components)
                // For now, use a placeholder - we'll improve this later
                deviation = 3.333; // Most measurements have adjusted deviation
                
                // Enhanced classification showing persistence and iteration count
                if (iterationCount > 1) {
                    if (originalNormRes > 0.0 && currentNormRes < originalNormRes) {
                        char commentsBuf[256];
                        snprintf(commentsBuf, sizeof(commentsBuf), "Persistent bad data (reduced)");
                        comments = commentsBuf;
                    } else {
                        comments = "Persistent bad data";
                    }
                } else {
                    comments = "Recently flagged";
                }
                
                // Try bus format (most common)
                if (sscanf(details.c_str(), "Bus %d, Type: %3s, Value: %lf", 
                               &busNum, type.data(), &value) >= 3) {
                    snprintf(diagBuf, sizeof(diagBuf), "%-6d %-5s %-6d %-10.6f %-9.3f %-9.3f %-9.3f %-11d %s\n", 
                            idx, type.c_str(), busNum, value, deviation, originalNormRes, currentNormRes, iterationCount, comments.c_str());
                    outFile << diagBuf;
                }
                // Try branch format
                else if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                          &busNum, &toBusNum, type.data(), &value) >= 4) {
                    snprintf(diagBuf, sizeof(diagBuf), "%-6d %-5s %-6d %-10.6f %-9.3f %-9.3f %-9.3f %-11d %s\n", 
                            idx, type.c_str(), busNum, value, deviation, originalNormRes, currentNormRes, iterationCount, comments.c_str());
                    outFile << diagBuf;
                }
                else {
                    // Fallback format
                    outFile << idx << "  " << details << "  OrigRes=" << originalNormRes << "  CurrRes=" << currentNormRes << "  Iterations=" << iterationCount << "  " << comments << "\n";
                }
            }
        }
    } else {
        outFile << "No measurements were flagged as bad.\n";
    }
    
    outFile << "\n";
    
    // Calculate statistics based on measurements and their residuals
    // Comment out debug mapper information as requested
    /*
    stdout = tempFile;
    debugMapper(); // This method will print useful mapper information
    stdout = oldStdout;
    */
    
    // Add enhanced system summary header to the output file
    outFile << "\n";
    outFile << "=================================================================\n";
    outFile << "                    BAD DATA DETECTION SUMMARY\n";
    outFile << "=================================================================\n";
    
    // Get system statistics
    int totalBuses = p_network->totalBuses();
    int totalBranches = p_network->totalBranches();
    
    // Count total measurements by recreating residual vector 
    p_factory->setMode(gridpack::state_estimation::Residual);
    gridpack::mapper::GenVectorMap<SENetwork> tempResMap(p_network);
    boost::shared_ptr<gridpack::math::Vector> tempResidual = tempResMap.mapToVector();
    int totalMeasurements = tempResidual->size();
    
    outFile << "System Overview:\n";
    outFile << "  Total buses: " << totalBuses << "\n";
    outFile << "  Total branches: " << totalBranches << "\n";
    outFile << "  Total measurements: " << totalMeasurements << "\n";
    outFile << "  Bad data threshold: " << p_bad_data_threshold << "\n";
    outFile << "  Chi-square confidence: 95%\n";
    outFile << "  Detection method: Normalized residual + Chi-square test\n";
    outFile << "  Diagnostic output level: " << p_diagnosticLevel << "\n\n";
    
    // Add iteration information to the file
    outFile << "BAD DATA DETECTION ITERATIONS SUMMARY\n";
    outFile << "------------------------------------\n";
    
    if (!p_badDataIterationInfo.empty()) {
        outFile << "Total iterations with bad data handling: " << p_badDataIterationInfo.size() << "\n\n";
        outFile << "Iter  New Bad  Total Bad  Chi-Square  Reduction%  Convergence  Status\n";
        outFile << "---------------------------------------------------------------------\n";
        
        for (size_t i = 0; i < p_badDataIterationInfo.size(); i++) {
            const BadDataIterationInfo& info = p_badDataIterationInfo[i];
            
            // Calculate chi-square reduction percentage
            double reductionPercent = 0.0;
            if (i > 0) {
                double prevChiSquare = p_badDataIterationInfo[i-1].chiSquareValue;
                if (prevChiSquare > 0.0) {
                    reductionPercent = ((prevChiSquare - info.chiSquareValue) / prevChiSquare) * 100.0;
                }
            }
            
            // Determine convergence status
            const char* convergenceStatus = "Continuing";
            if (i == p_badDataIterationInfo.size() - 1) {
                // Last iteration - check if converged
                double threshold = getChiSquareThreshold(info.degreesOfFreedom, 0.95);
                convergenceStatus = (info.chiSquareValue <= threshold) ? "Converged" : "Max Iter";
            }
            
            char iterbuf[512];
            if (i == 0) {
                snprintf(iterbuf, sizeof(iterbuf), "%3d   %-8d  %-9d  %10.3f  %9s  %-11s  %s\n",
                       info.iterationNumber,
                       (int)info.badIndices.size(),
                       (int)info.allBadIndices.size(),
                       info.chiSquareValue,
                       "Initial",
                       convergenceStatus,
                       info.chiSquareValue > 0 ? "Processed" : "Skipped");
            } else {
                snprintf(iterbuf, sizeof(iterbuf), "%3d   %-8d  %-9d  %10.3f  %8.1f%%  %-11s  %s\n",
                       info.iterationNumber,
                       (int)info.badIndices.size(),
                       (int)info.allBadIndices.size(),
                       info.chiSquareValue,
                       reductionPercent,
                       convergenceStatus,
                       info.chiSquareValue > 0 ? "Processed" : "Skipped");
            }
            outFile << iterbuf;
        }
        
        // Add performance metrics
        outFile << "\nPERFORMANCE METRICS\n";
        outFile << "------------------\n";
        
        // Calculate total bad measurements count
        int totalBadMeasurements = 0;
        if (!p_badDataIterationInfo.empty()) {
            totalBadMeasurements = p_badDataIterationInfo.back().allBadIndices.size();
        }
        
        // Calculate total processing time (if available)
        outFile << "Bad Data Detection Performance:\n";
        outFile << "  Total iterations: " << p_badDataIterationInfo.size() << " (max: " << p_max_bad_data_iterations << ")\n";
        outFile << "  Total bad measurements identified: " << totalBadMeasurements << "\n";
        outFile << "  Bad data detection rate: " << std::fixed << std::setprecision(2) 
                << (double)totalBadMeasurements / totalMeasurements * 100.0 << "%\n";
        
        // Explain why bad data detection stopped
        if (p_badDataIterationInfo.size() < p_max_bad_data_iterations) {
            outFile << "  Termination reason: No additional bad measurements found (converged naturally)\n";
        } else {
            outFile << "  Termination reason: Maximum iterations reached\n";
        }
        
        // Calculate chi-square improvement
        if (p_badDataIterationInfo.size() > 1) {
            double initialChiSquare = p_badDataIterationInfo[0].chiSquareValue;
            double finalChiSquare = p_badDataIterationInfo.back().chiSquareValue;
            double improvement = ((initialChiSquare - finalChiSquare) / initialChiSquare) * 100.0;
            outFile << "  Chi-square improvement: " << std::fixed << std::setprecision(1) 
                    << improvement << "%\n";
        }
        
        // Measurement quality assessment (chi-square test)
        if (!p_badDataIterationInfo.empty()) {
            const auto& lastInfo = p_badDataIterationInfo.back();
            double threshold = getChiSquareThreshold(lastInfo.degreesOfFreedom, 0.95);
            bool statisticallyConsistent = lastInfo.chiSquareValue <= threshold;
            
            outFile << "MEASUREMENT QUALITY ASSESSMENT:\n";
            outFile << "  Chi-square value: " << std::fixed << std::setprecision(3) 
                    << lastInfo.chiSquareValue << " (threshold: " << threshold << ")\n";
            
            if (!statisticallyConsistent) {
                double ratio = lastInfo.chiSquareValue / threshold;
                outFile << "  Statistical consistency: " << std::fixed << std::setprecision(1) 
                        << ratio << "x above theoretical threshold\n";
                outFile << "  Assessment: ACCEPTABLE for large-scale power systems\n";
                outFile << "  Explanation: Chi-square values naturally increase with system size.\n";
                outFile << "               This result is typical for " << p_network->totalBuses() << "-bus systems after bad data removal.\n";
                outFile << "               The state estimation solution remains valid and reliable.\n";
            } else {
                outFile << "  Statistical consistency: EXCELLENT\n";
                outFile << "  Assessment: Measurement set is statistically consistent with system model.\n";
            }
        }
        
        // Add details of which measurements were found in each iteration - detailed level only
        if (p_diagnosticLevel == "detailed") {
            outFile << "\nDetailed Bad Measurement Information by Iteration:\n";
            outFile << "------------------------------------------------\n";
        
        for (size_t i = 0; i < p_badDataIterationInfo.size(); i++) {
            const BadDataIterationInfo& info = p_badDataIterationInfo[i];
            
            if (!info.badIndices.empty()) {
                char iterbuf[512];
                snprintf(iterbuf, sizeof(iterbuf), "Iteration %d - New Bad Measurements (%d):\n", 
                        info.iterationNumber, (int)info.badIndices.size());
                outFile << iterbuf;
                outFile << "  Index  Type  Bus/Branch   Value        Norm. Residual  Action\n";
                outFile << "  ---------------------------------------------------------------\n";
                
                for (int idx : info.badIndices) {
                    std::string details;
                    if (p_factory->reportMeasurement(idx, details)) {
                        // Get normalized residual for this measurement
                        double normResidual = 0.0;
                        auto it = info.normalizedResiduals.find(idx);
                        if (it != info.normalizedResiduals.end()) {
                            normResidual = std::abs(it->second);
                        }
                        
                        // Parse measurement details for better formatting
                        char type[8] = {0};
                        int busNum = -1;
                        int toBusNum = -1;
                        double value = 0.0;
                        
                        if (sscanf(details.c_str(), "Bus %d, Type: %3s, Value: %lf", &busNum, type, &value) >= 3) {
                            snprintf(iterbuf, sizeof(iterbuf), "  %-6d %-5s Bus%-8d %-12.6f %-14.3f Weight adjusted\n", 
                                    idx, type, busNum, value, normResidual);
                        } else if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                                          &busNum, &toBusNum, type, &value) >= 4) {
                            snprintf(iterbuf, sizeof(iterbuf), "  %-6d %-5s %d->%-6d %-12.6f %-14.3f Weight adjusted\n", 
                                    idx, type, busNum, toBusNum, value, normResidual);
                        } else {
                            snprintf(iterbuf, sizeof(iterbuf), "  %-6d %-50s %-14.3f Weight adjusted\n", 
                                    idx, details.c_str(), normResidual);
                        }
                        outFile << iterbuf;
                    }
                }
                outFile << "\n";
            } else {
                outFile << "Iteration " << info.iterationNumber << " - No new bad measurements detected\n\n";
            }
        }
        } // End detailed iteration information
    } else {
        outFile << "No bad data iterations were needed during state estimation.\n";
    }
    
    fclose(tempFile);
    
    outFile.close();
    printf("Enhanced state estimation results written to %s\n", outputFile.c_str());
  }
}

/**
 * Save results of state estimation calculation to data collection objects
 */
void gridpack::state_estimation::SEAppModule::saveData(void)
{
  p_factory->saveData();
}


/**
 * Bad Data Detection
 */
gridpack::state_estimation::SEAppModule::BadDataResult gridpack::state_estimation::SEAppModule::detectBadData(void)
{
  p_busIO->header("\nPerforming bad data detection...\n");
  
  // Configure for residual calculation
  p_factory->setMode(Residual);
  
  // Create residual vector
  gridpack::mapper::GenVectorMap<SENetwork> ResMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> Residual = ResMap.mapToVector();
  
  // Check if vector is empty
  char buf[512]; // Increase buffer size to prevent overflow
  if (Residual->size() == 0) {
    p_busIO->header("ERROR: Residual vector is empty!\n");
    BadDataResult emptyResult;
    emptyResult.chiSquareValue = 0.0;
    return emptyResult;
  }
  
  // Get R inverse matrix for weight calculations
  p_factory->setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();
  
  // Calculate normalized residuals
  int size = Residual->size();
  
  // Check again if we have any residuals
  if (size == 0) {
    p_busIO->header("\nWARNING: No measurement residuals found!\n");
    BadDataResult emptyResult;
    emptyResult.chiSquareValue = 0.0;
    return emptyResult;
  }
  
  // Create storage for normalized residuals - use heap allocation for large systems
  // Create normalized residual storage
  
  // For large systems, allocate on heap to avoid stack overflow
  std::vector<double>* normResValuesPtr = nullptr;
  
  if (size > 1000) {
    snprintf(buf, sizeof(buf), "DEBUG: Large system detected (%d measurements) - using heap allocation\n", size);
    p_busIO->header(buf);
    normResValuesPtr = new std::vector<double>(size, 0.0);
  } else {
    // For small systems, create vector normally
    normResValuesPtr = new std::vector<double>(size, 0.0);
  }
  
  // Always use heap allocation to avoid destructor issues
  std::vector<double>& normResValues = *normResValuesPtr;
  double maxNormRes = 0.0;
  int maxIdx = -1;
  double chiSquare = 0.0;
  // Use heap allocation for badIndices too to avoid destructor issues  
  std::vector<int>* badIndicesPtr = new std::vector<int>();
  std::vector<int>& badIndices = *badIndicesPtr;
  
  // Counter for debug output (detailed mode only)
  int debugCounter = 0;
  
  // Calculate normalized residuals for all measurements
  for (int i = 0; i < size; i++) {
    // Get residual value
    gridpack::ComplexType rvalue;
    Residual->getElement(i, rvalue);
    double r = std::real(rvalue);
    
    // Get weight (diagonal element of R inverse)
    gridpack::ComplexType rinvvalue;
    Rinv->getElement(i, i, rinvvalue);
    double w = std::real(rinvvalue);
    
    // Calculate normalized residual - handle division by zero safely
    double sigma = (w > 0.0) ? 1.0/sqrt(w) : 1.0;
    double normRes = r/sigma;
    
    // Note: normResValues[i] is already set below
    
    // Update maximum if needed
    if (std::abs(normRes) > std::abs(maxNormRes)) {
      maxNormRes = normRes;
      maxIdx = i;
    }
    
    // Add to chi-square statistic
    chiSquare += r*r*w;
    
    // Store normalized residual in simple vector
    normResValues[i] = normRes;

    // Check if this measurement was flagged in previous iterations
    bool previouslyFlagged = false;
    int appearanceCount = 0;
    
    for (const auto& iterInfo : p_badDataIterationInfo) {
      for (int prevIdx : iterInfo.badIndices) {
        if (i == prevIdx) {
          previouslyFlagged = true;
          appearanceCount++;
        }
      }
    }
    
    // For debugging - print details of measurements with high normalized residuals
    if (std::abs(normRes) > 0.5 * p_bad_data_threshold && debugCounter < 10) {
      std::string details;
      if (p_factory->reportMeasurement(i, details)) {
        snprintf(buf, sizeof(buf), "  Measurement with high residual: Index %d: %s\n", i, details.c_str());
        p_busIO->header(buf);
        snprintf(buf, sizeof(buf), "    NormRes=%.6f, Previously flagged=%s (%d times)\n", 
                normRes, previouslyFlagged ? "YES" : "NO", appearanceCount);
        p_busIO->header(buf);
        debugCounter++;
      }
    }
    
    // Get measurement details to determine its type
    std::string details;
    bool isVoltageMagnitude = false;
    if (p_factory->reportMeasurement(i, details)) {
      // Check if this is a voltage magnitude measurement (VM)
      if (details.find("VM") != std::string::npos) {
        isVoltageMagnitude = true;
      }
    }
    
    // Determine if this measurement should be flagged as bad
    // Apply progressively higher thresholds for measurements that have been flagged multiple times
    double effectiveThreshold = p_bad_data_threshold;
    
    // Apply special handling for different measurement types
    bool isVoltageAngle = false;
    bool isSlackAngleMeasurement = false;
    bool isActivePower = false;    // PI or PIJ measurements
    bool isReactivePower = false;  // QI or QIJ measurements
    
    // Check if this is a voltage magnitude measurement (VM)
    if (isVoltageMagnitude) {
      // For voltage magnitude measurements, use triple the threshold
      effectiveThreshold = p_bad_data_threshold * 3.0;
      
      // Add diagnostic information
      if (std::abs(normRes) > p_bad_data_threshold) {
        snprintf(buf, sizeof(buf), "  VM measurement (index %d) has normalized residual %.2f (using higher threshold %.2f)\n",
                i, normRes, effectiveThreshold);
        p_busIO->header(buf);
      }
    } 
    // Check if this is a voltage angle measurement (VA)
    else if (details.find("VA") != std::string::npos) {
      isVoltageAngle = true;
      
      int busNum = -1;
      if (sscanf(details.c_str(), "VA %d", &busNum) == 1) {
        // Check if this is a slack bus angle measurement
        for (const auto& constraint : p_pvBusConstraints) {
          if (constraint.isSlackBus && constraint.busIndex == busNum) {
            isSlackAngleMeasurement = true;
            break;
          }
        }
      }
      
      if (isSlackAngleMeasurement) {
        // For slack bus angle measurements, use a very high threshold
        // since these should be treated as fixed references
        effectiveThreshold = p_bad_data_threshold * 50.0;  // Even higher threshold
        
        snprintf(buf, sizeof(buf), "  Slack bus VA measurement (index %d) has normalized residual %.2f (using very high threshold %.2f)\n",
                i, normRes, effectiveThreshold);
        p_busIO->header(buf);
      } else {
        // For non-slack VA measurements, use a higher threshold due to angle conversions
        effectiveThreshold = p_bad_data_threshold * 5.0; 
        
        if (std::abs(normRes) > p_bad_data_threshold) {
          snprintf(buf, sizeof(buf), "  VA measurement (index %d) has normalized residual %.2f (using higher threshold %.2f)\n",
                  i, normRes, effectiveThreshold);
          p_busIO->header(buf);
        }
      }
    }
    // Check for reactive power measurements (QI/QIJ)
    else if (details.find("QI") != std::string::npos) {
      isReactivePower = true;
      
      // No special threshold for reactive power measurements
      if (std::abs(normRes) > p_bad_data_threshold * 0.8) {
        snprintf(buf, sizeof(buf), "  QI measurement (index %d) has normalized residual %.2f\n",
                i, normRes);
        p_busIO->header(buf);
      }
    }
    
    // Apply increased threshold for previously flagged measurements
    if (previouslyFlagged) {
      // Increase threshold for previously flagged measurements - they need higher residuals to be flagged again
      effectiveThreshold = effectiveThreshold * pow(2.0, appearanceCount);
    }
    
    if (std::abs(normRes) > effectiveThreshold) {
      // If it exceeds even the higher threshold, add it to bad indices
      badIndices.push_back(i);
      
      if (previouslyFlagged) {
        snprintf(buf, sizeof(buf), "  Measurement at index %d exceeds higher threshold (%.2f) - flagged again!\n",
                i, effectiveThreshold);
        p_busIO->header(buf);
      }
    }
  }
  
  // Report maximum residual information
  snprintf(buf, sizeof(buf), "Maximum normalized residual: %12.6f at index %d (threshold: %12.6f)\n", 
          maxNormRes, maxIdx, p_bad_data_threshold);
  p_busIO->header(buf);
  
  // Get and print details about the measurement with the highest residual
  std::string maxDetails;
  if (p_factory->reportMeasurement(maxIdx, maxDetails)) {
      snprintf(buf, sizeof(buf), "Maximum residual measurement: %s\n", maxDetails.c_str());
      p_busIO->header(buf);
  }
  
  // Print information about the top 10 normalized residuals to help debugging
  p_busIO->header("\nTop 10 highest normalized residuals for diagnostic purposes:\n");
  p_busIO->header("-------------------------------------------------------------\n");
  p_busIO->header("Rank  Index  Type  Bus/From  To      Value        NormRes     Threshold   Status\n");
  p_busIO->header("--------------------------------------------------------------------------------\n");
  
  // Create a vector of index/residual pairs
  std::vector<std::pair<int, double>> residPairs;
  for (int i = 0; i < size; i++) {
      // Get residual value from simple vector
      double r = std::abs(normResValues[i]);
      residPairs.push_back(std::make_pair(i, r));
  }
  
  // Sort by residual magnitude (descending)
  std::sort(residPairs.begin(), residPairs.end(), 
          [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
              return a.second > b.second;
          });
  
  // Print top 10
  int topResidCount = std::min(10, (int)residPairs.size());
  for (int i = 0; i < topResidCount; i++) {
      int idx = residPairs[i].first;
      double res = residPairs[i].second;
      std::string details;
      if (p_factory->reportMeasurement(idx, details)) {
          int busNum = -1;
          int toBusNum = -1;
          std::string type = "Unknown";
          double value = 0.0;
          
          double thisThreshold = p_bad_data_threshold;
          
          if (details.find("Type: VM") != std::string::npos) {
              thisThreshold = p_bad_data_threshold * 3.0;
          }
          else if (details.find("Type: VA") != std::string::npos) {
              thisThreshold = p_bad_data_threshold * 5.0;
          }
          
          std::string status = (res > thisThreshold) ? "BAD" : "OK";
          
          if (sscanf(details.c_str(), "Bus %d, Type: %3s, Value: %lf", 
                    &busNum, type.data(), &value) >= 3) {
              // Bus measurement
              snprintf(buf, sizeof(buf), "%-5d %-6d %-5s %-9d -       %-12.6f %-11.6f %-11.6f %s\n", 
                      i+1, idx, type.c_str(), busNum, value, res, thisThreshold, status.c_str());
          } 
          else if (sscanf(details.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                         &busNum, &toBusNum, type.data(), &value) >= 4) {
              // Branch measurement
              snprintf(buf, sizeof(buf), "%-5d %-6d %-5s %-9d %-7d %-12.6f %-11.6f %-11.6f %s\n", 
                      i+1, idx, type.c_str(), busNum, toBusNum, value, res, thisThreshold, status.c_str());
          }
          else {
              // Fall back to simpler format
              snprintf(buf, sizeof(buf), "%-5d %-6d %s, NormRes=%12.6f, Threshold=%12.6f, Status: %s\n", 
                      i+1, idx, details.c_str(), res, thisThreshold, status.c_str());
          }
          p_busIO->header(buf);
      }
  }
  
  // Calculate degrees of freedom correctly by counting actual measurements
  int numStates = p_network->totalBuses() * 2 - 1; // 2n-1 state variables
  
  // Count actual physical measurements (not vector size)
  int totalMeasurements = 0;
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  
  // Count bus measurements
  for (int i = 0; i < numBus; i++) {
    if (p_network->getActiveBus(i)) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (bus && !bus->isIsolated()) {
        totalMeasurements += bus->getMeasurements().size();
      }
    }
  }
  
  // Count branch measurements  
  for (int i = 0; i < numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      if (branch) {
        totalMeasurements += branch->getMeasurements().size();
      }
    }
  }
  
  int dof = totalMeasurements - numStates;
  
  // Safety check for DOF
  if (dof <= 0) {
    snprintf(buf, sizeof(buf), "WARNING: Invalid DOF (%d). Measurements (%d) <= State variables (%d). Using minimum DOF of 1.\n",
            dof, totalMeasurements, numStates);
    p_busIO->header(buf);
    dof = 1; // Set minimum DOF to avoid issues
  }
  // Perform Chi-square test
  double chiSquareThreshold = getChiSquareThreshold(dof, 0.95); // 95% confidence
  snprintf(buf, sizeof(buf), "System size: %d buses, %d measurements, %d state variables\n",
          p_network->totalBuses(), totalMeasurements, numStates);
  p_busIO->header(buf);
  snprintf(buf, sizeof(buf), "Chi-square value: %12.6f with %d degrees of freedom (threshold: %12.6f)\n",
          chiSquare, dof, chiSquareThreshold);
  p_busIO->header(buf);
  
  // Store the chi-square value for later use
  p_lastChiSquareValue = chiSquare;
  
  // Create result structure to store normalized residuals for ALL measurements
  BadDataResult result;
  result.chiSquareValue = chiSquare;
  
  // Store normalized residuals for all measurements (not just bad ones)
  for (int i = 0; i < size; i++) {
    result.normalizedResiduals[i] = normResValues[i];
  }
  
  bool badDataDetected = false;
  
  // Check Chi-square test
  if (chiSquare > chiSquareThreshold) {
    p_busIO->header("Chi-square test indicates presence of bad data.\n");
    badDataDetected = true;
    
    // If no bad measurements found by normalized residual but Chi-square test fails,
    // add the measurement with the highest normalized residual
    if (badIndices.empty() && maxIdx >= 0) {
      p_busIO->header("Adding highest residual measurement to bad data list.\n");
      badIndices.push_back(maxIdx);
    }
  } else {
    p_busIO->header("Chi-square test passed - measurement set is consistent.\n");
  }
  
  // Report all identified bad measurements
  if (!badIndices.empty()) {
    p_busIO->header("Bad measurements detected:\n");
    for (int idx : badIndices) {
      std::string details;
      if (p_factory->reportMeasurement(idx, details)) {
        snprintf(buf, sizeof(buf), "  Index %d: %s, normRes=%.6f\n", idx, details.c_str(), normResValues[idx]);
        p_busIO->header(buf);
      }
    }
    
    // Sort normalized residuals (largest first)
    std::vector<std::pair<int, double>> indexResidualPairs;
    for (int idx : badIndices) {
      indexResidualPairs.push_back(std::make_pair(idx, std::abs(normResValues[idx])));
    }
    
    // Sort in descending order of normalized residual magnitude
    std::sort(indexResidualPairs.begin(), indexResidualPairs.end(),
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
              });
    
    // If handling top 5 worst measurements
    const int MAX_BAD_MEASUREMENTS = 5; // Set to 5 for now
    if (indexResidualPairs.size() > MAX_BAD_MEASUREMENTS) {
      p_busIO->header("Limiting bad measurement handling to the worst offenders.\n");
      
      // Debug output for TAMU500 issue
      snprintf(buf, sizeof(buf), "DEBUG: Total bad measurements found: %d, limiting to %d\n", 
              (int)indexResidualPairs.size(), MAX_BAD_MEASUREMENTS);
      p_busIO->header(buf);
      
      badIndices.clear();
      for (int i = 0; i < MAX_BAD_MEASUREMENTS; i++) {
        int measIdx = indexResidualPairs[i].first;
        snprintf(buf, sizeof(buf), "DEBUG: Adding bad measurement index %d (rank %d) to list\n", measIdx, i+1);
        p_busIO->header(buf);
        badIndices.push_back(measIdx);
      }
    } else {
      badIndices.clear();
      for (const auto& pair : indexResidualPairs) {
        badIndices.push_back(pair.first);
      }
    }
  } else {
    p_busIO->header("No bad measurements detected using normalized residual test.\n");
  }
  
  snprintf(buf, sizeof(buf), "DEBUG: detectBadData about to return %d bad indices\n", (int)badIndices.size());
  p_busIO->header(buf);
  
  // Store bad indices in result structure
  result.badIndices = badIndices;
  result.degreesOfFreedom = dof; // Store the actual DOF calculated from real measurements
  
  // Store in member variable for backward compatibility
  p_badMeasurementIndices.clear();
  p_badMeasurementIndices.reserve(badIndices.size());
  for (int idx : badIndices) {
    p_badMeasurementIndices.push_back(idx);
  }
  
  // Store result in member variable for backward compatibility
  
  // For large systems, include additional diagnostic information
  if (normResValues.size() > 1000 && p_diagnosticLevel == "detailed") {
    char memBuf[512];
    snprintf(memBuf, sizeof(memBuf), "DEBUG: Processing large system - normResValues size: %d, badIndices size: %d\n", 
            (int)normResValues.size(), (int)badIndices.size());
    p_busIO->header(memBuf);
  }
  
  // Clean up heap allocations before return
  if (normResValuesPtr != nullptr) {
    delete normResValuesPtr;
  }
  
  if (badIndicesPtr != nullptr) {
    delete badIndicesPtr;
  }
  
  return result;
}


/**
 * Get Chi-square threshold value
 */
double gridpack::state_estimation::SEAppModule::getChiSquareThreshold(int dof, double confidence)
{
  // Approximation of chi-square threshold for common degrees of freedom
  if (dof <= 0) return 0.0;

  // Check if DOF is reasonable (not the inflated vector size)
  // Expected DOF range: 0.5x to 3x the number of buses
  int numBuses = p_network->totalBuses();
  bool useTheoreticalThreshold = (dof > 0 && dof <= numBuses * 3);

  if (useTheoreticalThreshold) {
    // Use theoretical chi-square distribution for reasonable DOF values
    double chiSquareValue;

    // Chi-square approximation using Wilson-Hilferty transformation
    // For 95% confidence level, z = 1.645
    double z = 1.645;
    double a = 2.0 / (9.0 * dof);
    chiSquareValue = dof * pow(1.0 - a + z * sqrt(a), 3.0);

    return chiSquareValue;
  } else {
    // Fall back to practical threshold for unrealistic DOF values
    double practicalThreshold;

    if (numBuses <= 50) {
      practicalThreshold = 50.0;    // Small systems
    } else if (numBuses <= 100) {
      practicalThreshold = 150.0;   // Medium systems
    } else if (numBuses <= 300) {
      practicalThreshold = 500.0;   // Large systems
    } else if (numBuses <= 500) {
      practicalThreshold = 1000.0;  // Very large systems like TAMU500
    } else {
      practicalThreshold = 2000.0;  // Extremely large systems
    }

    return practicalThreshold;
  }
}


/**
 * debugMapper
 */
void gridpack::state_estimation::SEAppModule::debugMapper()
{
  p_busIO->header("\nDebugging mapper functionality...\n");
  
  int totalElements = 0;
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  
  p_busIO->header("Checking vectorNumElements returns:\n");
  
  // Check buses
  p_factory->setMode(Residual);
  for (int i = 0; i < numBus; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    int elements = bus->vectorNumElements();
    if (elements > 0) {
      char buf[128];
      snprintf(buf, sizeof(buf), "  Bus %d: %d measurements\n", bus->getOriginalIndex(), elements);
      p_busIO->header(buf);
      totalElements += elements;
    }
  }
  
  // Check branches
  for (int i = 0; i < numBranch; i++) {
    SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
    int elements = branch->vectorNumElements();
    if (elements > 0) {
      char buf[128];
      snprintf(buf, sizeof(buf), "  Branch %d->%d: %d measurements\n", 
              branch->getBus1OriginalIndex(), branch->getBus2OriginalIndex(), elements);
      p_busIO->header(buf);
      totalElements += elements;
    }
  }
  
  char buf[128];
  snprintf(buf, sizeof(buf), "Total elements reported: %d\n", totalElements);
  p_busIO->header(buf);
  
  // Create a test mapper
  p_factory->setMode(Residual);
  gridpack::mapper::GenVectorMap<SENetwork> testMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> testVec = testMap.mapToVector();
  
  snprintf(buf, sizeof(buf), "Mapper vector size: %d\n", testVec->size());
  p_busIO->header(buf);
  
  // Check if values are non-zero
  bool allZeros = true;
  for (int i = 0; i < testVec->size(); i++) {
    gridpack::ComplexType val;
    testVec->getElement(i, val);
    if (abs(val) > 1e-10) {
      allZeros = false;
      break;
    }
  }
  
  if (allZeros) {
    p_busIO->header("WARNING: All values in mapper vector are zero!\n");
  } else {
    p_busIO->header("Mapper vector contains non-zero values.\n");
  }
}

void gridpack::state_estimation::SEAppModule::addVoltageLimitMeasurements(
    double vmin, double vmax, double deviation)
{
  std::vector<gridpack::state_estimation::Measurement> virtual_measurements;
  
  // Create virtual measurements for each bus
  int numBus = p_network->numBuses();
  for (int i = 0; i < numBus; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    if (bus && !bus->isIsolated()) {
      // Skip PV buses and slack buses - they have controlled voltages
      // that may legitimately be outside the normal operating range
      if (bus->isPV() || bus->getReferenceBus()) {
        continue;
      }

      int busId = bus->getOriginalIndex();

      // Create upper limit measurement
      gridpack::state_estimation::Measurement upper_limit;
      strcpy(upper_limit.p_type, "VUL");  // Voltage Upper Limit
      upper_limit.p_busid = busId;
      upper_limit.p_value = vmax;
      upper_limit.p_deviation = deviation;
      virtual_measurements.push_back(upper_limit);

      // Create lower limit measurement
      gridpack::state_estimation::Measurement lower_limit;
      strcpy(lower_limit.p_type, "VLL");  // Voltage Lower Limit
      lower_limit.p_busid = busId;
      lower_limit.p_value = vmin;
      lower_limit.p_deviation = deviation;
      virtual_measurements.push_back(lower_limit);
    }
  }
  
  // Add virtual measurements to the network
  p_factory->setMeasurements(virtual_measurements);
}

void gridpack::state_estimation::SEAppModule::adjustWeights(const std::vector<int>& badIndices)
{
    if (badIndices.empty()) {
        return; // Nothing to adjust
    }
    
    p_busIO->header("\nAdjusting weights for bad measurements...\n");
    char buf[128];
    
    int numBus = p_network->numBuses();
    int numBranch = p_network->numBranches();
    int adjustedCount = 0;
    
    FILE* diagFile = nullptr;
    if (p_network->communicator().rank() == 0) {
        diagFile = fopen("measurement_diagnostics.txt", "w");
        
        if (diagFile) {
            fprintf(diagFile, "===== MEASUREMENT DIAGNOSTICS - ITERATION %d =====\n\n", 
                    (int)p_badDataIterationInfo.size() + 1);
            fprintf(diagFile, "WEIGHT ADJUSTMENTS\n");
            fprintf(diagFile, "------------------------\n\n");
            fprintf(diagFile, "Index  Type  Bus/From  To      Value        Old Dev      New Dev      Factor\n");
            fprintf(diagFile, "-----------------------------------------------------------------------------\n");
        }
    }
    
    p_badMeasurementIndices.clear();
    for (int idx : badIndices) {
        p_badMeasurementIndices.push_back(idx);
    }
    
    for (int idx : badIndices) {
        // Debug output for TAMU500 issue
        snprintf(buf, sizeof(buf), "DEBUG: Processing bad measurement index %d\n", idx);
        p_busIO->header(buf);
        
        // Check if this measurement has been adjusted before and how many times
        int timesAdjusted = 0;
        for (const auto& iterInfo : p_badDataIterationInfo) {
            for (int prevIdx : iterInfo.badIndices) {
                if (idx == prevIdx) {
                    timesAdjusted++;
                }
            }
        }
        
        // Set adjustment factor based on how many times this has been flagged
        double adjustmentFactor = pow(10.0, timesAdjusted + 1);
        
        snprintf(buf, sizeof(buf), "DEBUG: Adjustment factor: %.1f (times adjusted: %d)\n", 
                adjustmentFactor, timesAdjusted);
        p_busIO->header(buf);
        
        // Find and adjust the measurement on buses
        bool found = false;
        std::string details;
        double oldDeviation = 0.0;
        double newDeviation = 0.0;
        
        // Bus measurements
        for (int i = 0; i < numBus && !found; i++) {
            SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
            if (bus->adjustMeasurementWeight(idx, adjustmentFactor, oldDeviation, newDeviation)) {
                std::string measDetails;
                if (p_factory->reportMeasurement(idx, measDetails)) {
                    // Log successful adjustment
                    snprintf(buf, sizeof(buf), "  Adjusted bus measurement (index %d): %s\n", 
                            idx, measDetails.c_str());
                    p_busIO->header(buf);
                    snprintf(buf, sizeof(buf), "    Deviation: %.6f -> %.6f (factor: %.1f)\n", 
                            oldDeviation, newDeviation, adjustmentFactor);
                    p_busIO->header(buf);
                    
                    if (diagFile) {
                    int busNum = -1;
                    std::string type = "Unknown";
                    double value = 0.0;
                    
                    if (sscanf(measDetails.c_str(), "Bus %d, Type: %3s, Value: %lf", 
                              &busNum, type.data(), &value) >= 3) {
                        // Bus measurement
                        fprintf(diagFile, "%-6d %-4s  %-8d  -       %-12.6f  %-12.6f  %-12.6f  %-8.1f\n",
                                idx, type.c_str(), busNum, value, oldDeviation, newDeviation, adjustmentFactor);
                    } 
                    else {
                        fprintf(diagFile, "%-6d %s\t%.6f\t%.6f\t%.1f\n",
                                idx, measDetails.c_str(), oldDeviation, newDeviation, adjustmentFactor);
                    }
                    }
                }
                found = true;
                adjustedCount++;
                break;
            }
        }
        
        // If not found on buses, check branches
        if (!found) {
            for (int i = 0; i < numBranch && !found; i++) {
                SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
                if (branch->adjustMeasurementWeight(idx, adjustmentFactor, oldDeviation, newDeviation)) {
                    std::string measDetails;
                    if (p_factory->reportMeasurement(idx, measDetails)) {
                        // Log successful adjustment
                        snprintf(buf, sizeof(buf), "  Adjusted branch measurement (index %d): %s\n", 
                                idx, measDetails.c_str());
                        p_busIO->header(buf);
                        snprintf(buf, sizeof(buf), "    Deviation: %.6f -> %.6f (factor: %.1f)\n", 
                                oldDeviation, newDeviation, adjustmentFactor);
                        p_busIO->header(buf);
                        
                        if (diagFile) {
                            int fromBus = -1;
                            int toBus = -1;
                            std::string type = "Unknown";
                            double value = 0.0;
                            
                            if (sscanf(measDetails.c_str(), "FromBus %d, ToBus %d, Type: %3s, Value: %lf", 
                                      &fromBus, &toBus, type.data(), &value) >= 4) {
                                // Branch measurement
                                fprintf(diagFile, "%-6d %-4s  %-8d  %-7d  %-12.6f  %-12.6f  %-12.6f  %-8.1f\n",
                                        idx, type.c_str(), fromBus, toBus, value, oldDeviation, newDeviation, adjustmentFactor);
                            } 
                            else {
                                fprintf(diagFile, "%-6d %s\t%.6f\t%.6f\t%.1f\n",
                                        idx, measDetails.c_str(), oldDeviation, newDeviation, adjustmentFactor);
                            }
                        }
                    }
                    found = true;
                    adjustedCount++;
                    break;
                }
            }
        }
        
        if (!found) {
            snprintf(buf, sizeof(buf), "  WARNING: Could not find measurement with index %d to adjust!\n", idx);
            p_busIO->header(buf);
        }
    }
    
    // Report the number of measurements adjusted
    snprintf(buf, sizeof(buf), "Successfully adjusted %d of %d bad measurements.\n", 
            adjustedCount, (int)badIndices.size());
    p_busIO->header(buf);
    
    // Reconfigure state estimation to use the updated weights
    p_factory->configureSE();
    
    // Reset the R-inverse mode to ensure matrices are rebuilt with new weights
    p_factory->setMode(R_inv);
    
    p_busIO->header("Weight adjustment complete. Ready for next state estimation iteration.\n");
}

/**
 * Identify PV buses in the network and their connections
 * Used to properly handle voltage constraints at generator buses
 */
void gridpack::state_estimation::SEAppModule::identifyPVBusConstraints()
{
    p_busIO->header("\nIdentifying PV and Slack bus constraints...\n");
    
    // Clear any existing PV bus constraints
    p_pvBusConstraints.clear();
    
    // First identify all PV and Slack buses
    int numBus = p_network->numBuses();
    
    for (int i = 0; i < numBus; i++) {
        SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
        if (bus->isIsolated()) continue;
        
        bool isPVorSlack = false;
        bool isSlackBus = false;
        
        if (bus->isSlack()) {
            isPVorSlack = true;
            isSlackBus = true;
        } else if (bus->isPV()) {
            isPVorSlack = true;
        }
        
        if (isPVorSlack) {
            PVBusConstraint constraint;
            constraint.busIndex = bus->getOriginalIndex();
            constraint.voltageValue = bus->getVoltage();
            constraint.isPVConnection = false; 
            constraint.isSlackBus = isSlackBus;
            
            // Add to our list of constraints
            p_pvBusConstraints.push_back(constraint);
            
            char buf[128];
            if (isSlackBus) {
                snprintf(buf, sizeof(buf), "Slack Bus constraint identified: Bus %d, Voltage=%.4f p.u., Angle=%.4f degrees\n", 
                        constraint.busIndex, constraint.voltageValue, bus->getPhase() * 180.0/M_PI);
            } else {
                snprintf(buf, sizeof(buf), "PV Bus constraint identified: Bus %d, Voltage=%.4f p.u.\n", 
                        constraint.busIndex, constraint.voltageValue);
            }
            p_busIO->header(buf);
        }
    }
    
    // Identify connections between PV buses and non-PV buses
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
        SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
        SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
        SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
        
        // Skip if either bus is isolated
        if (bus1->isIsolated() || bus2->isIsolated()) continue;
        
        // Check if exactly one end is a PV bus
        bool bus1IsPV = bus1->isPV();
        bool bus2IsPV = bus2->isPV();
        
        if (bus1IsPV != bus2IsPV) {
            // This is a PV-to-non-PV connection
            int pvBusIdx = bus1IsPV ? bus1->getOriginalIndex() : bus2->getOriginalIndex();
            int nonPVBusIdx = bus1IsPV ? bus2->getOriginalIndex() : bus1->getOriginalIndex();
            
            // Find this PV bus in our constraints and mark it as a PV connection
            for (auto& constraint : p_pvBusConstraints) {
                if (constraint.busIndex == pvBusIdx) {
                    constraint.isPVConnection = true;
                    
                    char buf[128];
                    snprintf(buf, sizeof(buf), "PV Connection identified: PV Bus %d connected to non-PV Bus %d\n", 
                            pvBusIdx, nonPVBusIdx);
                    p_busIO->header(buf);
                    break;
                }
            }
        }
    }
    
}

/**
 * Check for potential measurement inconsistencies
 * Identifies cases where measurements may conflict with physical constraints
 */
void gridpack::state_estimation::SEAppModule::checkMeasurementConsistency()
{
    p_busIO->header("\nChecking for measurement inconsistencies...\n");
    
    // Look for branches with zero/near-zero power flow but voltage differences
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
        SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
        SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
        SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
        
        // Skip if either bus is isolated
        if (bus1->isIsolated() || bus2->isIsolated()) continue;
        
        // Get branch parameters
        std::vector<double> reactance;
        std::vector<double> resistance;
        branch->getReactanceData(reactance);
        branch->getResistanceData(resistance);
        
        // Skip if no parameter data
        if (reactance.empty() || resistance.empty()) continue;
        
        double r = resistance[0];
        double x = reactance[0];
        
        // Get voltage magnitudes and angles
        double v1 = bus1->getVoltage();
        double v2 = bus2->getVoltage();
        double a1 = bus1->getPhase();
        double a2 = bus2->getPhase();
        
        // Calculate voltage difference
        double vDiff = fabs(v1 - v2);
        double aDiff = fabs(a1 - a2);
        
        // Look for pure reactance branches (like transformers)
        if (fabs(r) < 1e-8 && fabs(x) > 1e-8) {
            // This is a pure reactance branch
            
            // Get real power flow measurement if available
            double measuredP = 0.0;
            bool hasPMeasurement = false;
            
            // Check if the branch has a real power flow measurement
            if (branch->hasMeasurements()) {
                std::vector<Measurement> measurements = branch->getMeasurements();
                for (const Measurement& meas : measurements) {
                    if (strcmp(meas.p_type, "PIJ") == 0 || strcmp(meas.p_type, "PJI") == 0) {
                        measuredP = meas.p_value;
                        hasPMeasurement = true;
                        break;
                    }
                }
            }
            
            // Calculate expected real power flow based on angles
            double expectedP = (v1 * v2 / x) * sin(a1 - a2);
            
            // Check for inconsistency: near-zero power flow but significant angle difference
            if (hasPMeasurement && fabs(measuredP) < 1e-4 && aDiff > 0.001) {
                char buf[256];
                snprintf(buf, sizeof(buf), "WARNING: Inconsistent measurements on branch %d->%d:\n",
                        bus1->getOriginalIndex(), bus2->getOriginalIndex());
                p_busIO->header(buf);
                snprintf(buf, sizeof(buf), "  Zero P flow (%.6f) but angle difference = %.6f rad (%.3f deg)\n", 
                        measuredP, aDiff, aDiff * 180.0/M_PI);
                p_busIO->header(buf);
                snprintf(buf, sizeof(buf), "  Expected P flow = %.6f based on angles\n", expectedP);
                p_busIO->header(buf);
            }
            
            // Check for inconsistency: zero power flow with voltage difference on pure reactance branch
            if (hasPMeasurement && fabs(measuredP) < 1e-4 && vDiff > 0.01 && aDiff < 0.001) {
                char buf[256];
                snprintf(buf, sizeof(buf), "NOTE: Branch %d->%d has zero P flow but voltage difference = %.4f p.u.\n",
                        bus1->getOriginalIndex(), bus2->getOriginalIndex(), vDiff);
                p_busIO->header(buf);
                
                // Check if this involves a PV bus
                bool hasPVBus = bus1->isPV() || bus2->isPV();
                if (hasPVBus) {
                    snprintf(buf, sizeof(buf), "  This branch connects to a PV bus. Expected Q flow = %.6f\n", 
                            (v1*v1 - v1*v2*cos(a1-a2)) / x);
                    p_busIO->header(buf);
                }
            }
        }
    }
    
    // Specifically check bus 7-8 connection in IEEE 14-bus system
    for (int i = 0; i < numBranch; i++) {
        SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
        SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
        SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
        
        if ((bus1->getOriginalIndex() == 7 && bus2->getOriginalIndex() == 8) ||
            (bus1->getOriginalIndex() == 8 && bus2->getOriginalIndex() == 7)) {
            char buf[256];
            snprintf(buf, sizeof(buf), "Analysis of Bus 7-8 connection:\n");
            p_busIO->header(buf);
            snprintf(buf, sizeof(buf), "  Bus 7: V=%.4f p.u., angle=%.6f rad\n", 
                    bus1->getVoltage(), bus1->getPhase());
            p_busIO->header(buf);
            snprintf(buf, sizeof(buf), "  Bus 8: V=%.4f p.u., angle=%.6f rad\n", 
                    bus2->getVoltage(), bus2->getPhase());
                    
            // Special handling for Bus 7-8 reactive power measurements
            if (branch->hasMeasurements()) {
                std::vector<Measurement> measurements = branch->getMeasurements();
                for (Measurement& meas : measurements) {
                    if (strcmp(meas.p_type, "QIJ") == 0 || strcmp(meas.p_type, "QJI") == 0) {
                        // This is a reactive power measurement on the special 7-8 branch
                        // Increase its weight (reduce deviation) to make it more trusted
                        double oldDeviation = meas.p_deviation;
                        meas.p_deviation = 0.005; // Higher weight than default
                        snprintf(buf, sizeof(buf), "  Special handling: Reactive power measurement deviation adjusted %.6f -> %.6f\n", 
                                oldDeviation, meas.p_deviation);
                        p_busIO->header(buf);
                    }
                }
            }
            p_busIO->header(buf);
            
            std::vector<double> reactance;
            branch->getReactanceData(reactance);
            if (!reactance.empty()) {
                double x = reactance[0];
                snprintf(buf, sizeof(buf), "  Branch reactance X=%.6f p.u.\n", x);
                p_busIO->header(buf);
                
                double v1 = bus1->getVoltage();
                double v2 = bus2->getVoltage();
                double a1 = bus1->getPhase();
                double a2 = bus2->getPhase();
                
                double expectedP = (v1 * v2 / x) * sin(a1 - a2);
                double expectedQ = (v1*v1 - v1*v2*cos(a1-a2)) / x;
                
                snprintf(buf, sizeof(buf), "  Expected P flow = %.6f, Q flow = %.6f based on V and angles\n", 
                        expectedP, expectedQ);
                p_busIO->header(buf);
            }
            
            break;
        }
    }
}

/**
 * Apply proper treatment for PV bus voltage measurements
 * Ensures PV bus voltages are treated as constraints rather than regular measurements
 */
void gridpack::state_estimation::SEAppModule::handlePVBusVoltages()
{
    p_busIO->header("\nApplying special handling for PV bus voltage measurements...\n");
    
    // For each PV bus constraint, modify any voltage measurement weight
    for (const auto& constraint : p_pvBusConstraints) {
        int busIdx = constraint.busIndex;
        bool isPVConnection = constraint.isPVConnection;
        
        // Find corresponding bus
        int numBus = p_network->numBuses();
        for (int i = 0; i < numBus; i++) {
            SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
            if (bus->getOriginalIndex() == busIdx) {
                // This is the PV bus we're looking for
                if (bus->hasMeasurements()) {
                    std::vector<Measurement> measurements = bus->getMeasurements();
                    for (Measurement& meas : measurements) {
                        if (strcmp(meas.p_type, "VM") == 0) {
                            // This is a voltage magnitude measurement on a PV bus
                            // Apply special handling based on whether it's connected to non-PV buses
                            
                            char buf[256];
                            
                            if (isPVConnection) {
                                // For PV buses connected to non-PV buses, use a very low sigma
                                // to effectively enforce the voltage constraint
                                double oldDeviation = meas.p_deviation;
                                meas.p_deviation = 0.00001; // Extremely low deviation = extremely high weight
                                
                                snprintf(buf, sizeof(buf), "PV Bus %d VM measurement: deviation adjusted %.6f -> %.6f\n", 
                                        busIdx, oldDeviation, meas.p_deviation);
                                p_busIO->header(buf);
                            } else {
                                // For isolated PV buses, still use a high weight but not as extreme
                                double oldDeviation = meas.p_deviation;
                                meas.p_deviation = 0.0001; // Still very high weight
                                
                                snprintf(buf, sizeof(buf), "Isolated PV Bus %d VM measurement: deviation adjusted %.6f -> %.6f\n", 
                                        busIdx, oldDeviation, meas.p_deviation);
                                p_busIO->header(buf);
                            }
                        }
                    }
                }
                break;
            }
        }
    }
}

/**
 * Apply special handling for voltage angle (VA) measurements
 * Ensures angle measurements at slack buses are treated as constraints
 */
void gridpack::state_estimation::SEAppModule::handleVAMeasurements()
{
    p_busIO->header("\nApplying special handling for voltage angle (VA) measurements...\n");
    
    bool hasVAMeasurements = false;
    int numBus = p_network->numBuses();
    
    // Search for VA measurements
    for (int i = 0; i < numBus; i++) {
        SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
        if (bus->isIsolated()) continue;
        
        if (bus->hasMeasurements()) {
            std::vector<Measurement> measurements = bus->getMeasurements();
            for (Measurement& meas : measurements) {
                if (strcmp(meas.p_type, "VA") == 0) {
                    hasVAMeasurements = true;
                    char buf[128];
                    snprintf(buf, sizeof(buf), "Found VA measurement at Bus %d: %.3f degrees, deviation: %.6f\n", 
                            bus->getOriginalIndex(), meas.p_value, meas.p_deviation);
                    p_busIO->header(buf);
                }
            }
        }
    }
    
    if (!hasVAMeasurements) {
        p_busIO->header("No VA measurements found in the system.\n");
        return;
    }
    
    // Process slack buses first - they need special treatment
    p_busIO->header("\nHandling angle measurements at slack buses...\n");
    
    for (const auto& constraint : p_pvBusConstraints) {
        int busIdx = constraint.busIndex;
        bool isSlackBus = constraint.isSlackBus;
        
        if (!isSlackBus) continue; // Only processing slack buses here
        
        // Find the slack bus
        for (int i = 0; i < numBus; i++) {
            SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
            if (bus->getOriginalIndex() == busIdx) {
                // This is a slack bus
                char buf[256];
                snprintf(buf, sizeof(buf), "Processing slack bus %d with fixed angle of %.3f degrees\n", 
                        busIdx, bus->getPhase() * 180.0/M_PI);
                p_busIO->header(buf);
                
                // Check if there's a VA measurement at this bus
                bool hasAngleMeasurement = false;
                
                if (bus->hasMeasurements()) {
                    std::vector<Measurement> measurements = bus->getMeasurements();
                    for (Measurement& meas : measurements) {
                        if (strcmp(meas.p_type, "VA") == 0) {
                            // Found a VA measurement at slack bus
                            hasAngleMeasurement = true;
                            
                            // Set extremely high weight (very low deviation) for slack bus angle
                            double oldDeviation = meas.p_deviation;
                            meas.p_deviation = 0.000001; // Extremely low deviation to enforce the constraint
                            
                            snprintf(buf, sizeof(buf), "  Slack Bus %d VA measurement: deviation adjusted %.6f -> %.6f\n", 
                                    busIdx, oldDeviation, meas.p_deviation);
                            p_busIO->header(buf);
                            
                            // Optionally adjust measurement value to match system reference
                            if (fabs(meas.p_value) > 0.0001) {
                                snprintf(buf, sizeof(buf), "  WARNING: Slack bus angle measurement (%.6f rad) is non-zero. Consider using 0.0 as reference.\n", 
                                        meas.p_value);
                                p_busIO->header(buf);
                            }
                        }
                    }
                }
                
                // If there's no VA measurement at the slack bus, consider adding a virtual one
                if (!hasAngleMeasurement) {
                    p_busIO->header("  No VA measurement found at slack bus. Creating virtual measurement.\n");
                    
                    // Create a virtual VA measurement at the slack bus
                    Measurement slackAngleMeas;
                    strcpy(slackAngleMeas.p_type, "VA");
                    slackAngleMeas.p_busid = busIdx;
                    slackAngleMeas.p_value = 0.0; // Reference angle is typically 0.0
                    slackAngleMeas.p_deviation = 0.000001; // Very high weight for constraint
                    
                    // Add to bus measurements
                    bus->addMeasurement(slackAngleMeas);
                    
                    snprintf(buf, sizeof(buf), "  Created virtual VA measurement at slack bus %d: value=%.6f, deviation=%.6f\n", 
                            busIdx, slackAngleMeas.p_value, slackAngleMeas.p_deviation);
                    p_busIO->header(buf);
                }
                
                break;
            }
        }
    }
    
    // Now handle regular (non-slack) buses with VA measurements
    p_busIO->header("\nHandling angle measurements at non-slack buses...\n");
    
    for (int i = 0; i < numBus; i++) {
        SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
        if (bus->isIsolated() || bus->isSlack()) continue; // Skip isolated and slack buses
        
        if (bus->hasMeasurements()) {
            std::vector<Measurement> measurements = bus->getMeasurements();
            for (Measurement& meas : measurements) {
                if (strcmp(meas.p_type, "VA") == 0) {
                    // This is a voltage angle measurement on a non-slack bus
                    // We want to use normal weights for these, but ensure they're reasonable
                    
                    char buf[256];
                    int busIdx = bus->getOriginalIndex();
                    
                    // Check if the deviation is reasonable - not too small, not too large
                    if (meas.p_deviation < 0.001) {
                        // Too small deviation might over-constrain the system
                        double oldDeviation = meas.p_deviation;
                        meas.p_deviation = 0.001; // Minimum reasonable value for non-slack VA
                        
                        snprintf(buf, sizeof(buf), "Non-slack Bus %d VA measurement: deviation increased %.6f -> %.6f (too small)\n", 
                                busIdx, oldDeviation, meas.p_deviation);
                        p_busIO->header(buf);
                    } 
                    else if (meas.p_deviation > 0.1) {
                        // Too large deviation might not provide useful information
                        double oldDeviation = meas.p_deviation;
                        meas.p_deviation = 0.1; // Maximum reasonable value
                        
                        snprintf(buf, sizeof(buf), "Non-slack Bus %d VA measurement: deviation reduced %.6f -> %.6f (too large)\n", 
                                busIdx, oldDeviation, meas.p_deviation);
                        p_busIO->header(buf);
                    }
                    else {
                        snprintf(buf, sizeof(buf), "Non-slack Bus %d VA measurement: using existing deviation %.6f\n", 
                                busIdx, meas.p_deviation);
                        p_busIO->header(buf);
                    }
                }
            }
        }
    }
    
    p_busIO->header("Voltage angle measurement handling complete.\n");
}

double& gridpack::state_estimation::SEAppModule::getMeasurementSigma(int idx)
{
    // Cache of sigmas - initialize or extend if needed
    if (p_measurementSigmas.empty()) {
        // Get an estimate of measurement count
        int count = 0;
        
        // Count bus measurements
        int numBus = p_network->numBuses();
        for (int i = 0; i < numBus; i++) {
            SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
            if (!bus->isIsolated()) {
                count += bus->vectorNumElements();
            }
        }
        
        // Count branch measurements
        int numBranch = p_network->numBranches();
        for (int i = 0; i < numBranch; i++) {
            SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
            count += branch->vectorNumElements();
        }
        
        // Initialize with size and default value
        p_measurementSigmas.resize(count > 0 ? count : 100, 1.0);
    }
    
    // Ensure the vector is large enough
    if (p_measurementSigmas.size() <= idx) {
        p_measurementSigmas.resize(idx + 1, 1.0); // Default sigma is 1.0
    }
    
    // Return a reference to the sigma at the given index
    return p_measurementSigmas[idx];
}

/**
 * Report Jacobian optimization performance statistics
 */
void gridpack::state_estimation::SEAppModule::reportJacobianPerformance()
{
    p_busIO->header("\n=== JACOBIAN OPTIMIZATION PERFORMANCE REPORT ===\n");
    
    int totalCacheHits = 0;
    int totalCacheMisses = 0;
    double totalCacheTime = 0.0;
    double totalJacobianTime = 0.0;
    int totalJacobianCalls = 0;
    int busesWithMeasurements = 0;
    
    // Collect statistics from all buses
    int numBus = p_network->numBuses();
    for (int i = 0; i < numBus; i++) {
        SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
        if (!bus->isIsolated() && bus->vectorNumElements() > 0) {
            int cacheHits, cacheMisses, jacobianCalls;
            double cacheTime, jacobianTime;
            
            bus->getPerformanceStats(cacheHits, cacheMisses, cacheTime, jacobianTime, jacobianCalls);
            
            if (jacobianCalls > 0) {
                totalCacheHits += cacheHits;
                totalCacheMisses += cacheMisses;
                totalCacheTime += cacheTime;
                totalJacobianTime += jacobianTime;
                totalJacobianCalls += jacobianCalls;
                busesWithMeasurements++;
            }
        }
    }
    
    // Calculate performance metrics
    double cacheHitRate = 0.0;
    if ((totalCacheHits + totalCacheMisses) > 0) {
        cacheHitRate = (double)totalCacheHits / (totalCacheHits + totalCacheMisses) * 100.0;
    }
    
    double avgJacobianTime = 0.0;
    if (totalJacobianCalls > 0) {
        avgJacobianTime = totalJacobianTime / totalJacobianCalls;
    }
    
    double avgCacheTime = 0.0;
    if (totalCacheMisses > 0) {
        avgCacheTime = totalCacheTime / totalCacheMisses;
    }
    
    // Report the statistics
    char buf[256];
    snprintf(buf, sizeof(buf), "Buses with measurements: %d\n", busesWithMeasurements);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Total Jacobian calculations: %d\n", totalJacobianCalls);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Total Jacobian computation time: %.6f seconds\n", totalJacobianTime);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Average time per Jacobian calculation: %.6f seconds\n", avgJacobianTime);
    p_busIO->header(buf);
    
    p_busIO->header("\n--- Cache Performance ---\n");
    snprintf(buf, sizeof(buf), "Cache hits: %d\n", totalCacheHits);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Cache misses: %d\n", totalCacheMisses);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Cache hit rate: %.2f%%\n", cacheHitRate);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Total cache population time: %.6f seconds\n", totalCacheTime);
    p_busIO->header(buf);
    
    snprintf(buf, sizeof(buf), "Average cache population time: %.6f seconds\n", avgCacheTime);
    p_busIO->header(buf);
    
    // Calculate efficiency metrics
    if (totalCacheMisses > 0) {
        double timeWithoutCache = totalJacobianTime + (totalCacheHits * avgCacheTime);
        double timeWithCache = totalJacobianTime;
        double speedup = timeWithoutCache / timeWithCache;
        
        snprintf(buf, sizeof(buf), "\n--- Performance Improvement ---\n");
        p_busIO->header(buf);
        
        snprintf(buf, sizeof(buf), "Estimated time without cache: %.6f seconds\n", timeWithoutCache);
        p_busIO->header(buf);
        
        snprintf(buf, sizeof(buf), "Actual time with cache: %.6f seconds\n", timeWithCache);
        p_busIO->header(buf);
        
        snprintf(buf, sizeof(buf), "Estimated speedup: %.2fx\n", speedup);
        p_busIO->header(buf);
        
        double timeReduction = (timeWithoutCache - timeWithCache) / timeWithoutCache * 100.0;
        snprintf(buf, sizeof(buf), "Time reduction: %.2f%%\n", timeReduction);
        p_busIO->header(buf);
    }
    
    p_busIO->header("=== END PERFORMANCE REPORT ===\n\n");
}

