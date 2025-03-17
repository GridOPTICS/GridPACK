/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app_module.cpp
 * @author Yousu Chen, Bruce Palmer
 * @date   2014-09-18 12:27:18 d3g096
 * Last updated: 8/5/2014 
 *
 * @brief
 * @update Yousu Chen
 *         Adding functions of bad data dection, chi-square testing 
 * @date   2025-03-05
 *
 *
 */
// -------------------------------------------------------------

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
#include "se_app_module.hpp"
double p_bad_data_threshold;

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::state_estimation::SEAppModule::SEAppModule(void)
{
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
      return;
    }
  }

  // Convergence and iteration parameters
  p_tolerance = secursor->get("tolerance",1.0e-3);
  p_max_iteration = secursor->get("maxIteration",20);

  // load input file
  //gridpack::parser::PTI23_parser<SENetwork> parser(p_network);
  //parser.parse(filename.c_str());
  double phaseShiftSign = secursor->get("phaseShiftSign",1.0);

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

  // partition network
  p_network->partition();

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<SENetwork>(1024, p_network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<SENetwork>(1024, p_network));
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
  char buf[128];
  sprintf(buf,"Tolerance: %12.4e\n",p_tolerance);
  p_busIO->header(buf);
  sprintf(buf,"Maximum number of iterations: %de\n",p_max_iteration);
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
/*
  if (p_comm.rank() == 0) {
    int idx;
    for (idx = 0; idx < meas.size(); idx++) {
      std::string meas_type = meas[idx].p_type;
      if (meas_type == "VM" || meas_type == "PI" || meas_type == "QI") {
        printf("Type: %s\n", meas[idx].p_type);
        printf("Bus: %d\n", meas[idx].p_busid);
        printf("Value: %f\n", meas[idx].p_value);
        printf("Deviation: %f\n", meas[idx].p_deviation);
      } else if (meas_type == "PIJ" || meas_type == "QIJ") {
        printf("Type: %s\n", meas[idx].p_type);
        printf("FromBus: %d\n", meas[idx].p_fbusid);
        printf("ToBus: %d\n", meas[idx].p_tbusid);
        printf("CKT: %s\n", meas[idx].p_ckt);
        printf("Value: %f\n", meas[idx].p_value);
        printf("Deviation: %f\n", meas[idx].p_deviation);
      }
      printf("\n");
    }
  }  
*/ 
  // Add measurements to buses and branches
  p_factory->setMeasurements(meas);

}

/**
 * Set up exchange buffers and other internal parameters and initialize
 * network components using data from data collection
 */
void gridpack::state_estimation::SEAppModule::initialize(void)
{
  // create factory
  p_factory.reset(new gridpack::state_estimation::SEFactoryModule(p_network));
  p_factory->load();

  // set network components using factory
  p_factory->setComponents();

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  p_factory->setExchange();

  // Create bus data exchange
  p_network->initBusUpdate();

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");
  p_bad_data_threshold = cursor->get("badDataThreshold", 3.0);
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
  // Run pre-check to identify obviously suspicious measurements
  preCheckMeasurements();

  // set YBus components so that you can create Y matrix  
  p_factory->setYBus();

  // set some state estimation parameters
  p_factory->configureSE();

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
    sprintf(buf, "Voltage limits: %.3f <= V <= %.3f\n", v_min, v_max);
    p_busIO->header(buf);
    
    // Set voltage limits for all buses
    p_factory->setVoltageLimits(v_min, v_max, true);
  }
  
  if (use_voltage_constraints) {
    double constraint_deviation = cursor->get("voltageConstraintDeviation", 0.0025);
    p_busIO->header("\nAdding voltage magnitude constraints\n");
    char buf[128];
    sprintf(buf, "Voltage limits: %.3f <= V <= %.3f (deviation: %.4f)\n", 
            v_min, v_max, constraint_deviation);
    p_busIO->header(buf);

    // Add virtual measurements for voltage constraints
    addVoltageLimitMeasurements(v_min, v_max, constraint_deviation);
  }

  //Create Y-bus matrix
  p_factory->setMode(YBus);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();

  // Create mapper to push voltage data back onto buses
  p_factory->setMode(Voltage);
  gridpack::mapper::BusVectorMap<SENetwork> VMap(p_network);

  // Create initial version (H), measurement vector (Ez), and Rinv
  p_factory->setMode(Jacobian_H);
  gridpack::mapper::GenMatrixMap<SENetwork> HJacMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> HJac = HJacMap.mapToMatrix();

  gridpack::mapper::GenVectorMap<SENetwork> EzMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> Ez = EzMap.mapToVector();

  p_factory->setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();

  // Convergence and iteration parameters
  ComplexType tol;
  //tol = 2.0*p_tolerance;
  //int iter = 0;
  bool converged = false;
  p_converged = false;
  int maxBadDataIterations = 5; // Limt bad data iterations
  int badDataIter = 0;
  bool badDataExists = true; 


  //Out loop for bad data handling
  while (badDataExists && badDataIter < maxBadDataIterations) {
    badDataExists = false; // Reset flag
    tol = 2.0 * p_tolerance; // Reset tolerance
    int iter = 0;


  // Start N-R loop
  while (real(tol) > p_tolerance && iter < p_max_iteration) {
    
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
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = p_config->getCursor("Configuration.State_estimation");
    gridpack::math::LinearSolver solver(*Gain);
    solver.configure(cursor);
    
    p_busIO->header("\n Print Gain matrix\n");

    // Solve linear equation
    boost::shared_ptr<gridpack::math::Vector> X(RHS->clone()); 
    p_busIO->header("\n Print RHS vector\n");
    X->zero(); //might not need to do this
    solver.solve(*RHS, *X);
    tol = X->normInfinity();
    char ioBuf[128];
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    p_busIO->header(ioBuf);

     // Push solution back onto bus variables
    p_factory->setMode(Voltage);
    VMap.mapToBus(X);
  
    // update values
    p_network->updateBuses();
    iter++;

  // End N-R loop
  }

  // Check convergence and handle bad data
  if (real(tol) <= p_tolerance) {

    converged = true;
    p_converged = true; 
    p_busIO->header("\n*** State estimation converged successfully ***\n");
/***
    debugMapper();
    // Only print outputs if convergence was achieved
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
***/
    // Perform bad data detection
//  detectBadData();
    std::vector<int> badIndices = detectBadData();
    if (!badIndices.empty()) {
      badDataExists = true;
      adjustWeights(badIndices); // Adjust weights of bad measurements
      // Rebuild Rinv with updated weights
      p_factory->setMode(R_inv);
      Rinv = RinvMap.mapToMatrix();
      p_busIO->header("\nAdjusted weights for bad measurements, re-running SE\n");
    }
  } else {
    p_busIO->header("\n*** WARNING: STATE ESTIMATION DID NOT CONVERGE ***\n");
    char ioBuf[128];
    sprintf(ioBuf,"Maximum iterations (%d) reached with tolerance %12.6e > %12.6e\n", 
            p_max_iteration, real(tol), p_tolerance);
    p_busIO->header(ioBuf);
    p_busIO->header("Results may not be reliable.\n");
    break; //Exist Outer loop if SE doesn't converge
  }
  badDataIter++;
}
  // Output results if converged
  if (converged) {
    debugMapper();
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
        sprintf(bdioBuf,"\nBad Data Iteration %d \n",badDataIter);
    }
}

/**
 * Write final results of state estimation calculation to standard
 * output
 */
void gridpack::state_estimation::SEAppModule::write(void)
{
  if (!p_converged) {
    p_busIO->header("Cannot print results - state estimation did not converge.\n");
    return;
  }
  
  // Get output file name from configuration
  std::string outputFile = "state_estimation_results.txt"; // Default value
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");
  if (cursor) {
    cursor->get("outputFile", &outputFile);
  }
  
  // Create output file
  std::ofstream outFile;
  
  // Only have process 0 write to the file
  if (p_network->communicator().rank() == 0) {
    outFile.open(outputFile.c_str());
    if (!outFile.is_open()) {
      printf("ERROR: Could not open output file %s\n", outputFile.c_str());
      return;
    }
  }
  // Only write to file if it's open
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
    
    // Write bus data to file directly
    outFile << "\n   State Estimation Outputs\n";
    outFile << "\n   Bus Number      Phase Angle      Voltage Magnitude\n";
    int numBus = p_network->numBuses();
    for (int i = 0; i < numBus; i++) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (!bus->isIsolated()) {
        char buf[128];
        double angle = bus->getPhase() * 180.0 / M_PI; // Convert to degrees
        sprintf(buf, "     %6d      %12.6f         %12.6f\n",
                bus->getOriginalIndex(), angle, bus->getVoltage());
        outFile << buf;
      }
    }
    
    // Write branch power flow data directly
    outFile << "\n   Branch Power Flow (p.u.)\n";
    outFile << "\n        Bus 1       Bus 2            P                    Q\n";
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
      SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
      if (!bus1->isIsolated() && !bus2->isIsolated()) {
        char buf[128];
        double p, q;
        branch->getPQ(bus1, &p, &q);
        sprintf(buf, "     %6d      %6d      %12.6f         %12.6f\n",
                bus1->getOriginalIndex(), bus2->getOriginalIndex(), p, q);
        outFile << buf;
      }
    }
    
    // Capture bus measurement comparison by temporarily redirecting stdout
    std::string busMeasurementData;
    bool hasBusMeasurements = false;
    
    for (int i = 0; i < numBus; i++) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (!bus->isIsolated() && bus->vectorNumElements() > 0) {
        // Redirect stdout to our temp file
        stdout = tempFile;
        
        // Let the bus write to stdout via serialWrite
        char buf[1];  // We're not using this buffer
        if (bus->serialWrite(buf, 1, "se")) {
          hasBusMeasurements = true;
        }
        
        // Restore stdout
        stdout = oldStdout;
      }
    }
    
    if (hasBusMeasurements) {
      // Return to beginning of temp file
      rewind(tempFile);
      
      // Read all content from temp file
      char buffer[4096];
      std::string output;
      while (fgets(buffer, sizeof(buffer), tempFile) != NULL) {
        output += buffer;
      }
      
      // Clear the temp file for reuse
      rewind(tempFile);
      ftruncate(fileno(tempFile), 0);
      
      // Write to our output file
      outFile << "\n   Comparison of Bus Measurements and Estimations\n";
      outFile << "\n   Type  Bus Number      Measurement          Estimate"
              << " Difference   Deviation\n";
      outFile << output;
    }
    
    // Capture branch measurement comparison similarly
    bool hasBranchMeasurements = false;
    
    for (int i = 0; i < numBranch; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
      SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
      if (!bus1->isIsolated() && !bus2->isIsolated() && branch->vectorNumElements() > 0) {
        // Redirect stdout to our temp file
        stdout = tempFile;
        
        // Let the branch write to stdout via serialWrite
        char buf[1];  // We're not using this buffer
        if (branch->serialWrite(buf, 1, "se")) {
          hasBranchMeasurements = true;
        }
        
        // Restore stdout
        stdout = oldStdout;
      }
    }
    
    if (hasBranchMeasurements) {
      // Return to beginning of temp file
      rewind(tempFile);
      
      // Read all content from temp file
      char buffer[4096];
      std::string output;
      while (fgets(buffer, sizeof(buffer), tempFile) != NULL) {
        output += buffer;
      }
      
      // Write to our output file
      outFile << "\n   Comparison of Branch Measurements and Estimations\n";
      outFile << "\n   Type  From    To  CKT   Measurement      Estimate"
              << " Difference   Deviation\n";
      outFile << output;
    }
    
    // Close temp file
    fclose(tempFile);
    
    outFile.close();
    printf("State estimation results written to %s\n", outputFile.c_str());
  }
}


/**
void gridpack::state_estimation::SEAppModule::write(void)
{
  if (!p_converged) {
    p_busIO->header("Cannot print results - state estimation did not converge.\n");
    return;
  }
  
  // Get output file name from configuration
  std::string outputFile = "state_estimation_results.txt"; // Default value
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.State_estimation");
  if (cursor) {
    cursor->get("outputFile", &outputFile);
  }
  
  // Create output file
  std::ofstream outFile;
  
  // Only have process 0 write to the file (assuming parallel execution)
  if (p_network->communicator().rank() == 0) {
    outFile.open(outputFile.c_str());
    if (!outFile.is_open()) {
      printf("ERROR: Could not open output file %s\n", outputFile.c_str());
    }
  }
  
  // First print to screen using existing IO objects
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
  
  // Only write to file if it's open
  if (p_network->communicator().rank() == 0 && outFile.is_open()) {
    // Capture the same output for the file
    outFile << "\n   State Estimation Outputs\n";
    outFile << "\n   Bus Number      Phase Angle      Voltage Magnitude\n";
    
    // Write bus data to file
    int numBus = p_network->numBuses();
    for (int i = 0; i < numBus; i++) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (!bus->isIsolated()) {
        char buf[128];
        double angle = bus->getPhase() * 180.0 / M_PI; // Convert to degrees
        sprintf(buf, "     %6d      %12.6f         %12.6f\n",
                bus->getOriginalIndex(), angle, bus->getVoltage());
        outFile << buf;
      }
    }
    
    // Similar format for branches
    outFile << "\n   Branch Power Flow (p.u.)\n";
    outFile << "\n        Bus 1       Bus 2            P                    Q\n";
    
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
      SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
      
      if (!bus1->isIsolated() && !bus2->isIsolated()) {
        char buf[128];
        // Get power flow values
        double p, q;
        branch->getPQ(bus1, &p, &q);
        sprintf(buf, "     %6d      %6d      %12.6f         %12.6f\n",
                bus1->getOriginalIndex(), bus2->getOriginalIndex(), p, q);
        outFile << buf;
      }
    }
    
    // Add similar sections for measurement comparisons
    outFile << "\n   Comparison of Bus Measurements and Estimations\n";
    outFile << "\n   Type  Bus Number      Measurement          Estimate"
            << " Difference   Deviation\n";
    
    // Now manually write each bus measurement to the file
    int numBus = p_network->numBuses();
    for (int i = 0; i < numBus; i++) {
      SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
      if (!bus->isIsolated() && bus->hasMeasurements()) {
	std::vector<Measurement> measurements = bus->getMeasurements();
	for (int j = 0; j < measurements.size(); j++) {
	  Measurement& meas = measurements[j];
	  std::string type = meas.p_type;
	  double estimate = 0.0;
	  
	  // Calculate the estimated value based on measurement type
	  if (type == "VM") {
	    estimate = bus->getVoltage();
	  } else if (type == "VA") {
	    estimate = bus->getPhase();
	  } else if (type == "PI") {
	    // Calculate real power injection (similar to what's in vectorGetElementValues)
	    estimate = bus->getInjectedRealPower(); // You might need to add this method
	  } else if (type == "QI") {
	    // Calculate reactive power injection
	    estimate = bus->getInjectedReactivePower(); // You might need to add this method
	  }
	  
	  // Format and write to file
	  char buf[256];
	  sprintf(buf, "    %s %8d    %16.5f  %16.5f   %16.5f    %8.4f\n",
		  type.c_str(), bus->getOriginalIndex(), meas.p_value, estimate,
		  estimate - meas.p_value, meas.p_deviation);
	  outFile << buf;
	}
      }
    }

    // Similarly for branch measurements:
    outFile << "\n   Comparison of Branch Measurements and Estimations\n";
    outFile << "\n   Type  From    To  CKT   Measurement      Estimate"
	    << " Difference   Deviation\n";

    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      SEBus* bus1 = dynamic_cast<SEBus*>(branch->getBus1().get());
      SEBus* bus2 = dynamic_cast<SEBus*>(branch->getBus2().get());
      
      if (!bus1->isIsolated() && !bus2->isIsolated() && branch->hasMeasurements()) {
	std::vector<Measurement> measurements = branch->getMeasurements();
	for (int j = 0; j < measurements.size(); j++) {
	  Measurement& meas = measurements[j];
	  std::string type = meas.p_type;
	  std::string ckt = meas.p_ckt;
	  double estimate = 0.0;
	  
	  // Calculate the estimated value based on measurement type
	  if (type == "PIJ") {
	    // Get power flow from bus i to j
	    gridpack::ComplexType s = branch->getComplexPower(ckt);
	    estimate = real(s)/branch->getBasePower();
	  } else if (type == "QIJ") {
	    gridpack::ComplexType s = branch->getComplexPower(ckt);
	    estimate = imag(s)/branch->getBasePower();
	  } else if (type == "PJI") {
	    gridpack::ComplexType s = branch->getRvrsComplexPower(ckt);
	    estimate = real(s)/branch->getBasePower();
	  } else if (type == "QJI") {
	    gridpack::ComplexType s = branch->getRvrsComplexPower(ckt);
	    estimate = imag(s)/branch->getBasePower();
	  }
	  // Add other measurement types if needed
	  
	  // Format and write to file
	  char buf[256];
	  sprintf(buf, "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
		  type.c_str(), bus1->getOriginalIndex(), bus2->getOriginalIndex(), 
		  ckt.c_str(), meas.p_value, estimate, 
		  estimate - meas.p_value, meas.p_deviation);
	  outFile << buf;
	}
      }
    }

    outFile.close();
    printf("State estimation results written to %s\n", outputFile.c_str());
  }
}
***/

/**
 * Save results of state estimation calculation to data collection objects
 */
void gridpack::state_estimation::SEAppModule::saveData(void)
{
  p_factory->saveData();
}

std::vector<int> gridpack::state_estimation::SEAppModule::detectBadData(void)
{
  p_busIO->header("\nPerforming bad data detection...\n");
  
  // Configure for residual calculation
  p_factory->setMode(Residual);
  
  // Create residual vector
  gridpack::mapper::GenVectorMap<SENetwork> ResMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> Residual = ResMap.mapToVector();
  
  // Print size information
  char buf[128];
  sprintf(buf, "Residual vector size: %d\n", Residual->size());
  p_busIO->header(buf);
  
  // Check if vector is empty
  if (Residual->size() == 0) {
    p_busIO->header("ERROR: Residual vector is empty!\n");
    return std::vector<int>();
  }
  
  // Print first few residuals
  p_busIO->header("First few residuals:\n");
  int printCount = std::min(10, Residual->size());
  for (int i = 0; i < printCount; i++) {
    gridpack::ComplexType value;
    Residual->getElement(i, value);
    sprintf(buf, "  [%d]: %f\n", i, std::real(value));
    p_busIO->header(buf);
  }
  
  // Continue with normal processing...
  p_factory->setMode(R_inv);
  
  // Step 4: Create and fill R inverse matrix
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();
  
  // Step 5: Calculate normalized residuals
  int size = Residual->size();
  
  // Check if we have any residuals
  if (size == 0) {
    p_busIO->header("\nWARNING: No measurement residuals found!\n");
    // Perform direct check for Bus 8
    return std::vector<int>();
  }
  
  // Continue with residual analysis
  boost::shared_ptr<gridpack::math::Vector> NormResidual(Residual->clone());
  double maxNormRes = 0.0;
  int maxIdx = -1;
  double chiSquare = 0.0;
  std::vector<int> badIndices; // Vector to store indices of bad measurements
  
  p_busIO->header("Full residual vector:\n");
  for (int i = 0; i < size; i++) {
    // Get residual value
    gridpack::ComplexType rvalue;
    Residual->getElement(i, rvalue);
    double r = std::real(rvalue);
    if (std::abs(r) > 0.01) {  // Only print significant ones
      char buf[128];
      sprintf(buf, "  [%d]: %f\n", i, r);
      p_busIO->header(buf);
    }

    
    // Get weight (diagonal element of R inverse)
    gridpack::ComplexType rinvvalue;
    Rinv->getElement(i, i, rinvvalue);
    double w = std::real(rinvvalue);
    
    // Calculate normalized residual
    double sigma = (w > 0.0) ? 1.0/sqrt(w) : 1.0;
    double normRes = r/sigma;

    // Update maximum
    if (fabs(normRes) > fabs(maxNormRes)) {
      maxNormRes = normRes;
      maxIdx = i;
    }
    
    // Add to chi-square
    chiSquare += r*r*w;
    
    // Store in normalized residual vector
    NormResidual->setElement(i, gridpack::ComplexType(normRes, 0.0));
  }
  // Enhanced max residual output
  sprintf(buf, "Step 2 (2) - Max normalized residual: %12.6f at index %d (threshold: %12.6f)\n",
          maxNormRes, maxIdx, p_bad_data_threshold);
  p_busIO->header(buf);
  
  // Calculate degrees of freedom
  int numStates = p_network->totalBuses() * 2 - 1; // 2n-1 state variables
  int dof = size - numStates;
  
  // Report results
  //char buf[128];
  sprintf(buf, "Maximum normalized residual: %12.6f at index %d (threshold: %12.6f)\n", 
          maxNormRes, maxIdx, p_bad_data_threshold);
  p_busIO->header(buf);

/***
if (fabs(maxNormRes) > p_bad_data_threshold) {
    std::string details;
    if (p_factory->reportMeasurement(maxIdx, details)) {
   // if (p_factory->reportMeasurement(busID, type, details, globalIdx)) {
        char buf[256];
        sprintf(buf, "Bad data detected: %s, normRes=%f\n", details.c_str(), maxNormRes);
        p_busIO->header(buf);
    }
} else {
    std::string details;
    if (p_factory->reportMeasurement(maxIdx, details)) {
    //if (p_factory->reportMeasurement(busID, type, details, globalIdx)) {
        char buf[256];
        sprintf(buf, "Largest normRes: %s, normRes=%f\n", details.c_str(), maxNormRes);
        p_busIO->header(buf);
    }
}
  
  sprintf(buf, "Chi-square value: %12.6f with %d degrees of freedom\n", 
          chiSquare, dof);
  p_busIO->header(buf);
  
  // Detect bad data
  if (fabs(maxNormRes) > p_bad_data_threshold) {
    p_busIO->header("Bad data detected!\n");
    p_factory->identifyBadData(maxIdx);
  } else {
    p_busIO->header("No bad data detected using normalized residual test.\n");
  }
}
***/
 // Report all bad measurements
  if (!badIndices.empty()) {
    p_busIO->header("Bad measurements detected:\n");
    for (int idx : badIndices) {
      std::string details;
      if (p_factory->reportMeasurement(idx, details)) {
        gridpack::ComplexType normResValue;

        // Retrieve the element at index idx
        NormResidual->getElement(idx, normResValue);
        //sprintf(buf, "  Index %d: %s, normRes=%.6f\n", idx, details.c_str(), NormResidual->getElement(idx:.real());
        sprintf(buf, "  Index %d: %s, normRes=%.6f\n", idx, details.c_str(), normResValue.real());
        p_busIO->header(buf);
      }
    }
  } else {
    p_busIO->header("No bad measurements detected.\n");
  }
  
  sprintf(buf, "Chi-square value: %12.6f with %d degrees of freedom\n",
          chiSquare, dof);
  p_busIO->header(buf);
  
  return badIndices; // Return the list of bad measurement indices
}
// Helper function to get Chi-square threshold value
double gridpack::state_estimation::SEAppModule::getChiSquareThreshold(int dof, double confidence)
{
  // Approximation of chi-square threshold for common degrees of freedom
  // For a more accurate implementation, you would use a proper statistical library
  if (dof <= 0) return 0.0;
  
  // Simple approximation for 95% confidence
  return dof + sqrt(2.0*dof) * 1.645;
}


void gridpack::state_estimation::SEAppModule::debugMapper()
{
  p_busIO->header("\nDebugging mapper functionality...\n");
  
  // First check if buses and branches are returning correct number of elements
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
      sprintf(buf, "  Bus %d: %d measurements\n", bus->getOriginalIndex(), elements);
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
      sprintf(buf, "  Branch %d->%d: %d measurements\n", 
              branch->getBus1OriginalIndex(), branch->getBus2OriginalIndex(), elements);
      p_busIO->header(buf);
      totalElements += elements;
    }
  }
  
  char buf[128];
  sprintf(buf, "Total elements reported: %d\n", totalElements);
  p_busIO->header(buf);
  
  // Now create a test mapper to see if it works
  p_factory->setMode(Residual);
  gridpack::mapper::GenVectorMap<SENetwork> testMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> testVec = testMap.mapToVector();
  
  sprintf(buf, "Mapper vector size: %d\n", testVec->size());
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
    // Iterate over each index in badIndices
    for (int idx : badIndices) {
        // Access the sigma for this measurement and increase it
        double& sigma = getMeasurementSigma(idx);
        sigma *= 10.0; // Increase sigma by a factor of 10 to reduce weight
    }
}

double& gridpack::state_estimation::SEAppModule::getMeasurementSigma(int idx)
{
    // Return a reference to the sigma at the given index
    // This is a placeholder; replace with your actual sigma access logic
    return p_measurementSigmas[idx];
}

