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
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/configuration/configuration.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/gen_matrix_map.hpp"
#include "gridpack/mapper/gen_vector_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "se_app_module.hpp"

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
        ret.push_back(measurement);
      } 
    }
  }
  return ret;
}


enum Parser{PTI23, PTI33};

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
    } else {
      printf("No network configuration file specified\n");
      return;
    }
  }

  // Convergence and iteration parameters
  p_tolerance = secursor->get("tolerance",1.0e-3);
  p_max_iteration = secursor->get("maxIteration",20);

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
}

/**
 * Solve the state estimation problem
 */
void gridpack::state_estimation::SEAppModule::solve(void)
{
  // set YBus components so that you can create Y matrix  
  p_factory->setYBus();

  // set some state estimation parameters
  p_factory->configureSE();

  p_factory->setMode(YBus);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();

  // Create mapper to push voltage data back onto buses
  p_factory->setMode(Voltage);
  gridpack::mapper::BusVectorMap<SENetwork> VMap(p_network);

  // Create initial version of  H Jacobian and estimation vector
  p_factory->setMode(Jacobian_H);
  gridpack::mapper::GenMatrixMap<SENetwork> HJacMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> HJac = HJacMap.mapToMatrix();

  gridpack::mapper::GenVectorMap<SENetwork> EzMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> Ez = EzMap.mapToVector();

  // Convergence and iteration parameters
  ComplexType tol;
  tol = 2.0*p_tolerance;
  int iter = 0;

  p_factory->setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();

  // Start N-R loop
  while (real(tol) > p_tolerance && iter < p_max_iteration) {

    
    // Form estimation vector
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
}

/**
 * Write final results of state estimation calculation to standard
 * output
 */
void gridpack::state_estimation::SEAppModule::write(void)
{

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

  // Output 
}

/**
 * Save results of state estimation calculation to data collection objects
 */
void gridpack::state_estimation::SEAppModule::saveData(void)
{
  p_factory->saveData();
}
