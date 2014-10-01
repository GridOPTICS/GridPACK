/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app.cpp
 * @author Yousu Chen 
 * @date   2014-09-18 12:27:18 d3g096
 * Last updated: 8/5/2014 
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "se_app.hpp"
#include "gridpack/serial_io/serial_io.hpp"

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::state_estimation::SEApp::SEApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::state_estimation::SEApp::~SEApp(void)
{
}

/**
 * Get list of measurements from external file
 * @param cursor pointer to contingencies in input deck
 * @return vector of measurements
 */
std::vector<gridpack::state_estimation::Measurement>
  gridpack::state_estimation::SEApp::getMeasurements(
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

/**
 * Execute application
 */
void gridpack::state_estimation::SEApp::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator p_comm;
  boost::shared_ptr<SENetwork> network(new SENetwork(p_comm));

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,p_comm);
  } else {
    config->open("input.xml",p_comm);
  }
  gridpack::utility::Configuration::CursorPtr cursor, secursor;
  secursor = config->getCursor("Configuration.State_estimation");
  std::string filename;
  if (!secursor->get("networkConfiguration",&filename)) {
     printf("No network configuration specified\n");
     return;
  }


  // load input file
  gridpack::parser::PTI23_parser<SENetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();


  ///////////////////////////////////////////////////////////////////
  // Read in measurement file
  std::string measurementfile;
  if (!secursor->get("measurementList", &measurementfile)) {
    measurementfile = "IEEE14_meas.xml";
  }
  bool ok = config->open(measurementfile, p_comm);

  // get a list of measurements
  cursor = config->getCursor("Measurements");
  gridpack::utility::Configuration::ChildCursors measurements;
  if (cursor) cursor->children(measurements);
  std::vector<gridpack::state_estimation::Measurement>
    meas = getMeasurements(measurements);

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
  ///////////////////////////////////////////////////////////////////


  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<SENetwork> busIO(128, network);
  gridpack::serial_io::SerialBranchIO<SENetwork> branchIO(256, network);
  char ioBuf[128];

  // create factory
  gridpack::state_estimation::SEFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  factory.setExchange();

  // Add measurements to buses and branches
  factory.setMeasurements(meas);

  // Create bus data exchange
  network->initBusUpdate();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  // set some state estimation parameters
  factory.configureSE();

  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  branchIO.header("\nybus:\n");
  ybus->print();

  // Create mapper to push voltage data back onto buses
  factory.setMode(Voltage);
  gridpack::mapper::BusVectorMap<SENetwork> VMap(network);

  // Create initial version of  H Jacobian and estimation vector
  factory.setMode(Jacobian_H);
  gridpack::mapper::GenMatrixMap<SENetwork> HJacMap(network);
  boost::shared_ptr<gridpack::math::Matrix> HJac = HJacMap.mapToMatrix();
  HJac->print();

  gridpack::mapper::GenVectorMap<SENetwork> EzMap(network);
  boost::shared_ptr<gridpack::math::Vector> Ez = EzMap.mapToVector();

  // Convergence and iteration parameters
  double tolerance = cursor->get("tolerance",1.0e-6);
  int max_iteration = cursor->get("maxIteration",20);
  ComplexType tol;

  tol = 2.0*tolerance;
  int iter = 0;

  factory.setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();
//  Rinv->print();

  // Start N-R loop
  while (real(tol) > tolerance && iter < max_iteration) {

    
    // Form estimation vector
    factory.setMode(Jacobian_H);
//    printf("Got to HJac\n");
    HJacMap.mapToMatrix(HJac);
//  HJac->print();
//  printf("Got to H'\n");

    // Form H'
    boost::shared_ptr<gridpack::math::Matrix> trans_HJac(transpose(*HJac));
//  trans_HJac->print();
//  printf("Got to Ez\n");

    // Build measurement equation
    EzMap.mapToVector(Ez);
//  Ez->print();
//  printf("Got to Gain\n");

    // Form Gain matrix
    boost::shared_ptr<gridpack::math::Matrix> Gain1(multiply(*trans_HJac, *Rinv));
    boost::shared_ptr<gridpack::math::Matrix> Gain(multiply(*Gain1, *HJac));
//  Gain->print();
//  printf("Got to H'*Rinv\n");

    // Form right hand side vector
    boost::shared_ptr<gridpack::math::Matrix> HTR(multiply(*trans_HJac, *Rinv));
//  HTR->print();
//  printf("Got to RHS\n");

//  printf("HTR iDim: %d jDim: %d Ez len: %d\n",HTR->rows(),HTR->cols(),Ez->size());
    boost::shared_ptr<gridpack::math::Vector> RHS(multiply(*HTR, *Ez));
//  printf("Create Solver\n");
//  RHS->print();
//  printf("Got to Solver\n");

    // create a linear solver
    gridpack::math::LinearSolver solver(*Gain);
    solver.configure(secursor);
//  printf("Got to DeltaX\n");

    // Solve linear equation
    boost::shared_ptr<gridpack::math::Vector> X(RHS->clone()); 
//    printf("Got to Solve\n");
    X->zero(); //might not need to do this
    solver.solve(*RHS, *X);
//  X->print();
//  printf("Got to updateBus\n");
//    boost::shared_ptr<gridpack::math::Vector> X(solver.solve(*RHS)); 
    tol = X->normInfinity();
    printf("Iteration %d Tol: %12.6e\n",iter+1,real(tol));

     // Push solution back onto bus variables
    factory.setMode(Voltage);
     //
    VMap.mapToBus(X);
  
    // update values
    network->updateBuses();
//  printf("Last sentence\n");

    iter++;

  // End N-R loop
  }

//  gridpack::serial_io::SerialBranchIO<SENetwork> branchIO(512,network);
//  branchIO.header("\n   Branch Power Flow\n");
//  branchIO.header("\n        Bus 1       Bus 2   CKT         P"
//                  "                    Q\n");
//  branchIO.write();


  busIO.header("\n   State Estimation Outputs\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();
  branchIO.header("\n   Branch Power Flow (p.u.)\n");
  branchIO.header("\n        Bus 1       Bus 2            P                    Q\n");
  branchIO.write();

  busIO.header("\n   Comparison of Bus Measurements and Estimations\n");
  busIO.header("\n   Type  Bus Number       Measurement       Estimate"
                 "           Difference          Deviation\n");
  busIO.write("se");

  branchIO.header("\n   Comparison of Branch Measurements and Estimations\n");
  branchIO.header("\n   Type     FromBus    ToBus       Measurement       Estimate"
                 "           Difference          Deviation\n");
  branchIO.write("se");

  // Output 

}

