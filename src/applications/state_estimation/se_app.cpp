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
/*        printf("%s %d %d %s %f %f\n", ckt.c_str(), measurement.p_fbusid,
          measurement.p_tbusid, meas_type.c_str(), measurement.p_value,
          measurement.p_deviation);*/
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

  gridpack::utility::CoarseTimer *timer =
  gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Total Application");
  timer->start(t_total);
 
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

  int t_pti = timer->createCategory("PTI Parser");
  timer->start(t_pti);
  // load input file
  gridpack::parser::PTI23_parser<SENetwork> parser(network);
  parser.parse(filename.c_str());
  timer->stop(t_pti);

  // partition network
  int t_part = timer->createCategory("Partition");
  timer->start(t_part);
  network->partition();
  timer->stop(t_part);

  ///////////////////////////////////////////////////////////////////
  // Read in measurement file
  std::string measurementfile;
  if (!secursor->get("measurementList", &measurementfile)) {
    measurementfile = "IEEE14_meas.xml";
  }
  bool ok = config->open(measurementfile, p_comm);

  // get a list of measurements
  int t_rmea = timer->createCategory("Read Meas");
  timer->start(t_rmea);
  cursor = config->getCursor("Measurements");
  gridpack::utility::Configuration::ChildCursors measurements;
  if (cursor) cursor->children(measurements);
  std::vector<gridpack::state_estimation::Measurement>
    meas = getMeasurements(measurements);
  timer->stop(t_rmea);
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
  ///////////////////////////////////////////////////////////////////

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<SENetwork> busIO(512, network);
  gridpack::serial_io::SerialBranchIO<SENetwork> branchIO(512, network);
  char ioBuf[128];

  // create factory
  gridpack::state_estimation::SEFactory factory(network);
  int t_load = timer->createCategory("Factory: Load");
  timer->start(t_load);
  factory.load();
  timer->stop(t_load);


  // set network components using factory
  int t_setc = timer->createCategory("Factory: Set Components");
  timer->start(t_setc);
  factory.setComponents();
  timer->stop(t_setc);

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  int t_setx = timer->createCategory("Factory: Set Exchange");
  timer->start(t_setx);
  factory.setExchange();
  timer->stop(t_setx);

  // Add measurements to buses and branches
  int t_setm = timer->createCategory("Factory: Set Measurement");
  timer->start(t_setm);
  factory.setMeasurements(meas);
  timer->stop(t_setm);

  // Create bus data exchange
  int t_updt = timer->createCategory("Bus Update");
  timer->start(t_updt);
  network->initBusUpdate();

  if (factory.checkLoneBus(busIO.getStream().get())) {
  //  sprintf(ioBuf,"\nIsolated bus found for contingency %s\n",
  //      contingency.p_name.c_str());
  //  busIO.header(ioBuf);
    // Exchange bus status between processors
    network->updateBuses();
  }
  timer->stop(t_updt);

  // set YBus components so that you can create Y matrix  
  int t_fact = timer->createCategory("Factory setYBus");
  timer->start(t_fact);
  factory.setYBus();
  timer->stop(t_fact);

  // set some state estimation parameters
  int t_factcfgSE = timer->createCategory("Factory configureSE");
  timer->start(t_factcfgSE);
  factory.configureSE();
  timer->stop(t_factcfgSE);

  int t_mapper = timer->createCategory("Mapper");
  timer->start(t_mapper);
  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  branchIO.header("\nybus:\n");
 // ybus->print();

  // Create mapper to push voltage data back onto buses
  factory.setMode(Voltage);
  gridpack::mapper::BusVectorMap<SENetwork> VMap(network);

  // Create initial version of  H Jacobian and estimation vector
  factory.setMode(Jacobian_H);
  gridpack::mapper::GenMatrixMap<SENetwork> HJacMap(network);
  boost::shared_ptr<gridpack::math::Matrix> HJac = HJacMap.mapToMatrix();
  //printf(" ------------HJacMap.mapToMatrix------\n");
//  branchIO.header("\nHJac:\n");
//  HJac->print();

  gridpack::mapper::GenVectorMap<SENetwork> EzMap(network);
  boost::shared_ptr<gridpack::math::Vector> Ez = EzMap.mapToVector();
  timer->stop(t_mapper);

  // Convergence and iteration parameters
  double tolerance = cursor->get("tolerance",1.0e-5);
  int max_iteration = cursor->get("maxIteration",20);
  ComplexType tol;

  tol = 2.0*tolerance;
  int iter = 0;

  int t_rinv = timer->createCategory("R_Inv");
  timer->start(t_rinv);
  factory.setMode(R_inv);
  gridpack::mapper::GenMatrixMap<SENetwork> RinvMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Rinv = RinvMap.mapToMatrix();
  timer->stop(t_rinv);
//  busIO.header("\n Print Rinv vector outside\n");
//  Rinv->print();

  // Start N-R loop
  int t_solv = timer->createCategory("Solver NRloop");
  timer->start(t_solv);

  while (real(tol) > tolerance && iter < max_iteration) {

    
    // Form estimation vector
    int t_H = timer->createCategory("Jacobian_H");
    timer->start(t_H);
    factory.setMode(Jacobian_H);
    printf("Got to HJac\n");
    HJacMap.mapToMatrix(HJac);
    busIO.header("\n Print HJac matrix\n");
    timer->stop(t_H);
//    HJac->print();
//  printf("Got to H'\n");

    // Form H'
    int t_HT = timer->createCategory("Jacobian_HT");
    timer->start(t_HT);
    boost::shared_ptr<gridpack::math::Matrix> trans_HJac(transpose(*HJac));
    busIO.header("\n Print HJac' matrix\n");
    timer->stop(t_HT);
//    trans_HJac->print();
//  printf("Got to Ez\n");

    // Build measurement equation
    EzMap.mapToVector(Ez);
    busIO.header("\n Print Ez Vector\n");
//    Ez->print();
//  printf("Got to Gain\n");

    // Form Gain matrix
    int t_gain = timer->createCategory("Gain");
    timer->start(t_gain);
    boost::shared_ptr<gridpack::math::Matrix> Gain1(multiply(*trans_HJac, *Rinv));
    boost::shared_ptr<gridpack::math::Matrix> Gain(multiply(*Gain1, *HJac));
    timer->stop(t_gain);
//  Gain->print();
//  printf("Got to H'*Rinv\n");

    // Form right hand side vector
    int t_rhs = timer->createCategory("RHS");
    timer->start(t_rhs);
    boost::shared_ptr<gridpack::math::Matrix> HTR(multiply(*trans_HJac, *Rinv));
//    busIO.header("\n Print HTR Vector\n");
//    HTR->print();
//  printf("Got to RHS\n");

//    printf("HTR iDim: %d jDim: %d Ez len: %d\n",HTR->rows(),HTR->cols(),Ez->size());
    boost::shared_ptr<gridpack::math::Vector> RHS(multiply(*HTR, *Ez));
    timer->stop(t_rhs);
//  printf("Create Solver\n");
//  RHS->print();
//  printf("Got to Solver\n");

    // create a linear solver
    int t_sol = timer->createCategory("SE Solver");
    timer->start(t_sol);
    gridpack::math::LinearSolver solver(*Gain);

    solver.configure(secursor);
    
//    busIO.header("\n Print Gain matrix\n");
//    Gain->print();
//  printf("Got to DeltaX\n");

    // Solve linear equation
    boost::shared_ptr<gridpack::math::Vector> X(RHS->clone()); 
//    printf("Got to Solve\n");
    busIO.header("\n Print RHS vector\n");
//    RHS->print();
    X->zero(); //might not need to do this
    int t_soleq = timer->createCategory("Pure SE Solver");
    timer->start(t_soleq);
    solver.solve(*RHS, *X);
    timer->stop(t_sol);
    timer->stop(t_soleq);
//  X->print();
//  printf("Got to updateBus\n");
//    boost::shared_ptr<gridpack::math::Vector> X(solver.solve(*RHS)); 
    tol = X->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    busIO.header(ioBuf);

    int t_up = timer->createCategory("SE Update");
    timer->start(t_up);
     // Push solution back onto bus variables
    factory.setMode(Voltage);
     //
    VMap.mapToBus(X);
  
    // update values
    network->updateBuses();
    timer->stop(t_up);
//  printf("Last sentence\n");

    iter++;

  // End N-R loop
  }
  timer->stop(t_solv);

//  gridpack::serial_io::SerialBranchIO<SENetwork> branchIO(512,network);
//  branchIO.header("\n   Branch Power Flow\n");
//  branchIO.header("\n        Bus 1       Bus 2   CKT         P"
//                  "                    Q\n");
//  branchIO.write();

  int t_io = timer->createCategory("IO");
  timer->start(t_io);

  busIO.header("\n   State Estimation Outputs\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();
  branchIO.header("\n   Branch Power Flow (p.u.)\n");
  branchIO.header("\n        Bus 1       Bus 2            P                    Q\n");
  branchIO.write();

  busIO.header("\n   Comparison of Bus Measurements and Estimations\n");
  busIO.header("\n   Type  Bus Number      Measurement          Estimate"
                 "         Difference   Deviation\n");
  busIO.write("se");

  branchIO.header("\n   Comparison of Branch Measurements and Estimations\n");
  branchIO.header("\n   Type      From        To  CKT      Measurement          Estimate"
                 "         Difference   Deviation\n");
  timer->stop(t_io);
  branchIO.write("se");

  // Output 
  timer->stop(t_total);
  timer->dump();

}

