/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

//#define USE_TIMESTAMP

#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/parser/PSSE_seq_parser.hpp"
//#include "gridpack/mapper/full_map.hpp"
//#include "gridpack/mapper/bus_vector_map.hpp"
//#include "gridpack/math/math.hpp"
#include "gridpack/parallel/global_vector.hpp"
#include "dsf_app_module.hpp"
//#include "gridpack/component/base_component.hpp"
//#include "hadrec_app_module.hpp"
#include <iostream>
#include <string>
#include <vector>
#include "gridpack/utilities/string_utils.hpp"

using namespace std;

#ifdef USE_HELICS
//#include "helics/ValueFederates.hpp"
//#include <helics/shared_api_library/ValueFederate.h>
#include <helics/helics.hpp>
#endif

//#define MAP_PROFILE

// Calling program for dynamic simulation application

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DSFullApp::DSFullApp(void)
{
  p_internal_watch_file_name = false; 
  p_generatorWatch = false;
  p_loadWatch = false;
  p_generators_read_in = false;
  p_save_time_series = false;
  p_monitorGenerators = false;
  p_bDynSimuDone = false;
  p_suppress_watch_files = false;
  Simu_Current_Step = 0;

  bapplyLineTripAction = false;
  bapplyLoadChangeP = false;
  bapplyLoadChangeQ = false;
  p_report_dummy_obs = false;
  p_biterative_solve_network = false;
  p_iterative_network_debug = false;
  p_generator_observationpower_systembase = true;
  ITER_TOL = 1.0e-7;
  MAX_ITR_NO = 8;

  p_current_time = 0.0;
  p_time_step = 0.005;
  
}

/**
 * Basic constructor with commmunicator argument
 * @param comm communicator that application object is restricted to
 */
gridpack::dynamic_simulation::DSFullApp::DSFullApp(gridpack::parallel::Communicator comm)
  : p_comm(comm)
{
  p_internal_watch_file_name = false; 
  p_generatorWatch = false;
  p_loadWatch = false;
  p_generators_read_in = false;
  p_save_time_series = false;
  p_monitorGenerators = false;
  p_bDynSimuDone = false;
  p_suppress_watch_files = false;
  Simu_Current_Step = 0;
  
  bapplyLineTripAction = false;
  bapplyLoadChangeP = false;
  bapplyLoadChangeQ = false;
  p_report_dummy_obs = false;
  p_biterative_solve_network = false;
  p_iterative_network_debug = false;
  ITER_TOL = 1.0e-7;
  MAX_ITR_NO = 8;
  
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSFullApp::~DSFullApp(void)
{
}

enum Format{PTI23, PTI33, PTI34, PTI35};
/**
 * Read in and partition the dynamic simulation network. The input file is read
 * directly from the Dynamic_simulation block in the configuration file so no
 * external file names or parameters need to be passed to this routine
 * @param network pointer to a DSFullNetwork object. This should not have any
 * buses or branches defined on it.
 * @param config pointer to open configuration file
 * @param otherfile name of network configuration file if different from the
 * one in the input deck
 */
void gridpack::dynamic_simulation::DSFullApp::readNetwork(
    boost::shared_ptr<DSFullNetwork> &network,
    gridpack::utility::Configuration *config,
    const char *otherfile)
{
  p_comm = network->communicator();
  p_network = network;
  p_config = config;

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  std::string filename;
  int filetype = PTI23;
  if (otherfile == NULL) {
    if (cursor->get("networkConfiguration",&filename)) {
      filetype = PTI23;
    } else if (cursor->get("networkConfiguration_v33",&filename)) {
      filetype = PTI33;
    } else if (cursor->get("networkConfiguration_v34",&filename)) {
      filetype = PTI34;
    } else if (cursor->get("networkConfiguration_v35",&filename)) {
      filetype = PTI35;
    } else {
      printf("No network configuration specified\n");
    }
  } else {
    filetype = PTI23;
    filename = otherfile;
  }

  p_sim_time = cursor->get("simulationTime",0.0);
  if (p_sim_time == 0.0) {
    // TODO: some kind of error
  }
  p_time_step = cursor->get("timeStep",0.0);
  if (p_time_step == 0.0) {
    // TODO: some kind of error
  }
  
  //--------------whether iteratively compute network current-----------
  p_biterative_solve_network = cursor->get("iterativeNetworkInterface",false);
  p_iterative_network_debug = cursor->get("iterativeNetworkInterfaceDebugPrint",false);
  
  ITER_TOL = cursor->get("iterativeNetworkInterfaceTol", 1.0e-7);
  MAX_ITR_NO = cursor->get("iterativeNetworkInterfaceMaxItrNo", 8);
  
  //printf ("-----rk debug in gridpack::dynamic_simulation::DSFullApp::readNetwork( ): ITER_TOL: %15.12f, MAX_ITR_NO: %d \n\n", ITER_TOL, MAX_ITR_NO);

  // Monitor generators for frequency violations
  p_monitorGenerators = cursor->get("monitorGenerators",false);
  p_report_dummy_obs = cursor->get("reportNonExistingElements",false);
  p_maximumFrequency = cursor->get("frequencyMaximum",61.8);

  // load input file
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else if (filetype == PTI34) {
    gridpack::parser::PTI34_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else if (filetype == PTI35) {
    gridpack::parser::PTI35_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else {
    printf("Unknown filetype\n");
  }
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  filename = cursor->get("generatorParameters","");

  // partition network
  network->partition();
  p_analytics.reset(new gridpack::analysis::NetworkAnalytics<DSFullNetwork>(network));

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(512, network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<DSFullNetwork>(128, network));
}

/**
 * Assume that DSFullNetwork already exists and just cache an internal pointer
 * to it. This routine does not call the partition function. Also read in
 * simulation parameters from configuration file
 * @param network pointer to a complete DSFullNetwork object.
 * @param config pointer to open configuration file
 */
void gridpack::dynamic_simulation::DSFullApp::setNetwork(
    boost::shared_ptr<DSFullNetwork> &network,
    gridpack::utility::Configuration *config)
{
  p_comm = network->communicator();
  p_network = network;
  p_config = config;

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  /*std::string filename;
  if (!cursor->get("networkConfiguration",&filename)) {
    printf("No network configuration specified\n");
  }*/
  std::string filename = cursor->get("networkConfiguration", 
      "No network configuration specified");
  p_sim_time = cursor->get("simulationTime",0.0);
  if (p_sim_time == 0.0) {
    // TODO: some kind of error
  }
  p_time_step = cursor->get("timeStep",0.0);
  if (p_time_step == 0.0) {
    // TODO: some kind of error
  }
  
  //--------------whether iteratively compute network current-----------
  p_biterative_solve_network = cursor->get("iterativeNetworkInterface",false);
  p_iterative_network_debug = cursor->get("iterativeNetworkInterfaceDebugPrint",false);
  p_generator_observationpower_systembase = cursor->get("generatorObservationPowerSystemBase",true);
  
  ITER_TOL = cursor->get("iterativeNetworkInterfaceTol", 1.0e-7);
  MAX_ITR_NO =  cursor->get("iterativeNetworkInterfaceMaxItrNo", 8);
  
  //printf ("-----rk debug in gridpack::dynamic_simulation::DSFullApp::setNetwork( ): ITER_TOL: %15.12f, MAX_ITR_NO: %d \n\n", ITER_TOL, MAX_ITR_NO);

  // Monitor generators for frequency violations
  p_monitorGenerators = cursor->get("monitorGenerators",false);
  p_report_dummy_obs = cursor->get("reportNonExistingElements",false);
  p_maximumFrequency = cursor->get("frequencyMaximum",61.8);
  p_analytics.reset(new gridpack::analysis::NetworkAnalytics<DSFullNetwork>(network));

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(512, network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<DSFullNetwork>(128, network));
}

/**
 * Read generator parameters. These will come from a separate file (most
 * likely). The name of this file comes from the input configuration file.
 * @param ds_idx index of dyr file if a list of dyr files are provided.
 */
void gridpack::dynamic_simulation::DSFullApp::readGenerators(int ds_idx)
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  std::string filename;
  gridpack::parser::PTI23_parser<DSFullNetwork> parser(p_network);
  if (ds_idx == -1) {
    cursor = p_config->getCursor("Configuration.Dynamic_simulation");
    filename = cursor->get("generatorParameters","");
  } else if (ds_idx >= 0) {
    gridpack::utility::Configuration::CursorPtr dyr_cursor;
    dyr_cursor = p_config->getCursor(
        "Configuration.Dynamic_simulation.generatorFiles");
    gridpack::utility::Configuration::ChildCursors files;
    if (dyr_cursor) dyr_cursor->children(files);
    if (ds_idx < files.size()) {
      if (!files[ds_idx]->get("generatorParams",&filename)) {
        printf("Unknown generator parameter file specified\n");
        return;
      }
    }
  }
  //printf("p[%d] generatorParameters: %s\n",p_comm.rank(),filename.c_str());
  if (filename.size() > 0) parser.externalParse(filename.c_str());
  //printf("p[%d] finished Generator parameters\n",p_comm.rank());
}

/**
 * Read sequence data from a file.
 */
void gridpack::dynamic_simulation::DSFullApp::readSequenceData()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  std::string filename = cursor->get("sequenceDataFile","");
  if (filename != "") {
    gridpack::parser::PSSE_seq_parser<DSFullNetwork> parser(p_network);
    parser.parse(filename);
  }
}

/**
 * Check to see if system is secure
 */
int gridpack::dynamic_simulation::DSFullApp::isSecure()
{
  return p_insecureAt;
}

/**
 * Set up exchange buffers and other internal parameters and
 * initialize
 * network components using data from data collection
 */
void gridpack::dynamic_simulation::DSFullApp::initialize()
{
  // create factory
  p_factory.reset(new gridpack::dynamic_simulation::DSFullFactory(p_network));
  // p_factory->dumpData();
  p_factory->load();

  // set network components using factory
  p_factory->setComponents();
  
  // set voltages for the extended buses from composite load model
  p_factory->setExtendedCmplBusVoltage();
  
  // load parameters for the extended buses from composite load model
  p_factory->LoadExtendedCmplBus();

  // set YBus components so that you can create Y matrix  
  p_factory->setYBus();

  if (!p_factory->checkGen()) {
    p_busIO->header("Missing generators on at least one processor\n");
    return;
  }
}

/**
 * Reinitialize calculation from data collections
 */
void gridpack::dynamic_simulation::DSFullApp::reload()
{
  orgYbus.reset();
  ybusyl.reset();
  ybuspg.reset();
  ybus_jxd.reset();
  ybus.reset();
  ybus_fy.reset();
  ybus_posfy.reset();

  ybusMap_sptr.reset();
  ngenMap_sptr.reset();
  nbusMap_sptr.reset();

  volt.reset();
  INorton_full.reset();
  INorton_full_chk.reset();
  volt_full.reset();

  solver_sptr.reset();
  solver_fy_sptr.reset();
  solver_posfy_sptr.reset();

  p_factory->load();
  p_factory->setYBus();
}

/**
 * Execute the time integration portion of the application
 */
void gridpack::dynamic_simulation::DSFullApp::solve(
    gridpack::dynamic_simulation::Event fault)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

  t_solve = timer->createCategory("DS Solve: Total");
  t_misc = timer->createCategory("DS Solve: Miscellaneous");
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif
  timer->start(t_solve);
  timer->start(t_misc);

  // Get cursor for setting solver options
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  timer->stop(t_misc);

  t_mode = timer->createCategory("DS Solve: Set Mode");
  timer->start(t_mode);
  p_factory->setMode(YBUS);
  timer->stop(t_mode);
  t_ybus = timer->createCategory("DS Solve: Make YBus");
  timer->start(t_ybus);
  
  ybusMap_sptr.reset(new gridpack::mapper::FullMatrixMap<DSFullNetwork> (p_network));
  orgYbus = ybusMap_sptr->mapToMatrix();
  
  //printf("\n=== org ybus: ============\n");
  //orgYbus->print();
  //orgYbus->save("ybus_GridPACK_org.m");
  //exit(0);

  //p_factory->addLoadAdmittance();

  // Form constant impedance load admittance yl for all buses and add it to
  // system Y matrix: ybus = ybus + yl
  p_factory->setMode(YL);
  ybusyl = ybusMap_sptr->mapToMatrix();
  timer->stop(t_ybus);
  //branchIO.header("\n=== ybus after added yl: ============\n");
  //printf("\n=== ybus after added yl: ============\n");
  //ybusyl->print();
  //ybusyl->save("ybus_GridPACK_yl.m");
  //exit(0);

  p_factory->setMode(PG);
  ybuspg = ybusMap_sptr->mapToMatrix();
  //printf("\n=== ybus after added pg: ============\n");
  //ybuspg->print();
  //exit(0);

  //printf("# of buses in the network: %d\n", p_network->totalBuses());

  // Add j*Xd' to system Y matrix:
  // Extract appropriate xdprime and xdpprime from machine data
  timer->start(t_mode);
  p_factory->setMode(jxd);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybus_jxd = ybusMap_sptr->mapToMatrix();
  //branchIO.header("\n=== ybusyl after added j*Xd': =============\n");
  //printf("\n=== ybusyl after added j*Xd': =============\n");
  //ybus_jxd->print();
  //ybus_jxd->save("ybus_GridPACK_jxd.m");
  
  
  // Add dynamic load impedance to system Y matrix:
  timer->start(t_mode);
  p_factory->setMode(YDYNLOAD);
  timer->stop(t_mode);
  ybus = ybusMap_sptr->mapToMatrix();
  //branchIO.header("\n=== ybus_jxd after added dynamic load impedance': =============\n");
  //printf("\n=== ybus_dynload after added dynamic load impedance': =============\n");
  //ybus->print();
  //ybus->save("ybus_GridPACK_dynload.m");
  
  //exit(0);

  // Compute ybus_fy for fault on stage
  ybus_fy.reset(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_misc);
  p_factory->setEvent(fault);
  timer->stop(t_misc);
  timer->start(t_mode);
  p_factory->setMode(onFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap_sptr->overwriteMatrix(ybus_fy);
  //branchIO.header("\n=== ybus_fy: ============\n");
  //printf("\n=== ybus_fy: ============\n");
  //ybus_fy->print();
  //ybus_fy->save("ybus_fy_GridPACK_jxd.m");

  // Compute ybus_posfy for fault clear stage
  ybus_posfy.reset(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_mode);
  p_factory->setMode(posFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap_sptr->incrementMatrix(ybus_posfy);
  //branchIO.header("\n=== ybus_posfy: ============\n");
  //printf("\n=== ybus_posfy: ============\n");
  //ybus_posfy->print();
  //ybus_posfy->save("ybus_posfy_GridPACK_jxd.m");
  timer->stop(t_ybus);

  // Simulation related variables
  t_init = timer->createCategory("DS Solve: Initialization");
  timer->start(t_init);
  
  int t_step[20];
  double t_width[20];

  //const double sysFreq = 60.0;
  //double pi = 4.0*atan(1.0);
  //const double basrad = 2.0 * pi * sysFreq;
  //gridpack::ComplexType jay(0.0, 1.0);

  // switch info is set up here
  int nswtch = 4;
  static double sw1[4];
  static double sw7[4];
  sw1[0] = 0.0;
  sw1[1] = fault.start;
  sw1[2] = fault.end;
  sw1[3] = p_sim_time;
  sw7[0] = p_time_step;
  sw7[1] = p_time_step;
  sw7[2] = p_time_step;
  sw7[3] = p_time_step;
  simu_total_steps = 0;
  for (int i = 0; i < nswtch-1; i++) {
    t_step[i] = (int) ((sw1[i+1] -sw1[i]) / sw7[i]);
    t_width[i] = (sw1[i+1] - sw1[i]) / t_step[i];
    simu_total_steps += t_step[i];
  }
  simu_total_steps++;
  
  // Initialize vectors for integration 
  p_factory->initDSVect(p_time_step);
  
  p_factory->setGeneratorObPowerBaseFlag(p_generator_observationpower_systembase);
  //exit(0);

  ngenMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork> (p_network));
  
  // Map to create vector volt
  volt = ngenMap_sptr->mapToVector();
  //p_busIO->header("\n=== volt: ===\n");
  //volt->print();

  solver_sptr.reset(new gridpack::math::LinearSolver (*ybus));
  solver_sptr->configure(cursor);
  
  //gridpack::math::LinearSolver solver_fy(*ybus_fy);
  solver_fy_sptr.reset(new gridpack::math::LinearSolver (*ybus_fy));
  solver_fy_sptr->configure(cursor);
  
  //gridpack::math::LinearSolver solver_posfy(*ybus_posfy);
  //gridpack::math::LinearSolver solver_posfy(*ybus); 
  solver_posfy_sptr.reset(new gridpack::math::LinearSolver (*ybus));
  solver_posfy_sptr->configure(cursor);

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagP = 0;
  flagC = 0;
  S_Steps = 1;
  last_S_Steps = -1;

  p_insecureAt = -1;

  p_factory->setMode(make_INorton_full);
  //gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
  nbusMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork>(p_network));
  INorton_full = nbusMap_sptr->mapToVector();
  INorton_full_chk = nbusMap_sptr->mapToVector();
  max_INorton_full = 0.0;
  volt_full.reset(INorton_full->clone());

  timer->stop(t_init);
  if (!p_suppress_watch_files) {
#ifdef USE_TIMESTAMP
    if (p_generatorWatch) p_generatorIO->header("t, t_stamp");//bus_id,ckt,x1d_1,x2w_1,x3Eqp_1,x4Psidp_1,x5Psiqpp_1");
    //#  if (p_generatorWatch) p_generatorIO->header("t, t_stamp,bus_id,ckt,x1d_1,x2w_1,x3Eqp_1,x4Psidp_1,x5Psiqpp_1");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    if (p_generatorWatch) p_generatorIO->header("\n");

    if (p_loadWatch) p_loadIO->header("t, t_stamp");
    if (p_loadWatch) p_loadIO->write("load_watch_header");
    if (p_loadWatch) p_loadIO->header("\n");
#else
    if (p_generatorWatch) p_generatorIO->header("t");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    if (p_generatorWatch) p_generatorIO->header("\n");

    if (p_loadWatch) p_loadIO->header("t");
    if (p_loadWatch) p_loadIO->write("load_watch_header");
    if (p_loadWatch) p_loadIO->header("\n");
#endif
  }
#ifdef USEX_GOSS
  if (p_generatorWatch) p_generatorIO->dumpChannel();
  if (p_loadWatch) p_loadIO->dumpChannel();
#endif
  p_frequencyOK = true;
  // Save initial time step
  //saveTimeStep();
 
	
 
#ifdef USE_HELICS
	//std::cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << std::endl;
	cout << "-------------!!!helics test: HELICS Version: " << helics::versionString << endl;
	string configFile = "/home/huan495/gridpack-dev/src/build/applications/dynamic_simulation_full_y/testcase/helics_39bus_3.json";
    helics::ValueFederate fed(configFile);
	helics::Publication pub;
	helics::Input sub;
	double helics_requestTime = 0.0;
	
	//to get publication definitions
    int pubCount = fed.getPublicationCount();
	
	printf("-------------helics test: num of pub: %d \n", pubCount);
    for(int i = 0; i < pubCount; i++) {
        pub = fed.getPublication(i);
        string pubInfo = pub.getInfo();
        // do stuff to tie pub to GridPACK object property
    }
    
	//to get subscription definitions
    int subCount = fed.getInputCount();
	printf("-------------helics test: num of sub: %d \n", subCount);
	
    for(int i = 0; i < subCount; i++) {
        sub = fed.getInput(i);
        string subInfo = sub.getInfo();
        // do stuff to tie pub to GridPACK object property
    }
         
	//let helics broker know you are ready to start simulation 
	fed.enterExecutingMode();	

#endif  //end if of HELICS


  for (Simu_Current_Step = 0; Simu_Current_Step < simu_total_steps - 1; Simu_Current_Step++) {
  //for (Simu_Current_Step = 0; Simu_Current_Step < 200; Simu_Current_Step++) {
    //char step_str[128];
    //sprintf(step_str,"\nIter %d\n", Simu_Current_Step);
    //p_busIO->header(step_str);
    timer->start(t_misc);
    printf("Step %d\ttime %5.3f sec: \n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
    //printf("\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
    ///char step_str[128];
    ///sprintf(step_str, "\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
     ///p_busIO->header(step_str);
    S_Steps = Simu_Current_Step;

    if (Simu_Current_Step < steps1) {
      flagP = 0;
      flagC = 0;
    } else if (Simu_Current_Step == steps1) {
      flagP = 0;
      //flagC = 1;
      flagC = 0;
    } else if ((Simu_Current_Step > steps1) && (Simu_Current_Step < steps2)) {
      flagP = 1;
      flagC = 1;
    } else if (Simu_Current_Step == steps2) {
      flagP = 1;
      //flagC = 2;
      flagC = 1;
    } else if (Simu_Current_Step > steps2) {
      flagP = 2;
      flagC = 2;
    }
    timer->stop(t_misc);
    
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor_currentInjection(false);
    } else {
      p_factory->predictor_currentInjection(true);
    }

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    t_mIf = timer->createCategory("DS Solve: Modified Euler Predictor: Make INorton");
    timer->start(t_mIf);
	p_factory->setMode(make_INorton_full);
    nbusMap_sptr->mapToVector(INorton_full);
    ///gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
    ///boost::shared_ptr<gridpack::math::Vector> INorton_full = nbusMap_sptr->mapToVector();
    //p_busIO->header("\n=== [Predictor] INorton_full: ===\n");
    //printf("renke test \n=== [Predictor] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_mIf);
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif
 
    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2) 
    // to calculate terminal volt: ----------
    t_psolve = timer->createCategory("DS Solve: Modified Euler Predictor: Linear Solver");
    timer->start(t_psolve);
    //boost::shared_ptr<gridpack::math::Vector> volt_full(INorton_full->clone());
    volt_full->zero();
#if 0
    bool flag_chk = true;
    while (flag_chk == true ) {
		
			volt_full->zero();
			
			if (flagP == 0) {
				solver_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 1) {
				solver_fy_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 2) {
				solver_posfy_sptr->solve(*INorton_full, *volt_full);
			}
			

			printf("1: itr test:----previous predictor_INorton_full:\n");
			INorton_full->print();

			INorton_full_chk->equate(*INorton_full);
			printf("2: itr test:----predictor_INorton_full_chk:\n");
			INorton_full_chk->print();

			nbusMap_sptr->mapToBus(volt_full);
			p_factory->setVolt(false);
			
			if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
				p_factory->predictor_currentInjection(false);
			} else {
				p_factory->predictor_currentInjection(true);
			}
			
# if 0	
			printf("3: itr test:----previous predictor_INorton_full:\n");
			INorton_full->print();

			INorton_full_chk->equate(*INorton_full);
			printf("4: itr test:----predictor_INorton_full_chk:\n");
			INorton_full_chk->print();
# endif			
			p_factory->setMode(make_INorton_full);
			nbusMap_sptr->mapToVector(INorton_full);
			
			printf("5: itr test:----predictor_INorton_full:\n");
			INorton_full->print();
			
			//multiply(*ybus_fy, *volt_full, *INorton_full_chk);
			INorton_full_chk->add(*INorton_full, -1.0);
			max_INorton_full=abs(INorton_full_chk->normInfinity());
			
			if (max_INorton_full <1.0e-8) {
				flag_chk = false;
			} else {
				
				printf("max_INorton_full = %8.4f \n", max_INorton_full);
				//printf("-----INorton_full : \n");
				//INorton_full->print();
				//printf("-----INorton_full_chk - INorton_full : \n");
				//INorton_full_chk->print();
			}
    }
#else
    if (flagP == 0) {
      solver_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy_sptr->solve(*INorton_full, *volt_full);
    }
#endif
    timer->stop(t_psolve);

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    //p_busIO->header("\n=== [Predictor] volt_full: ===\n");
    //volt_full->print();
    //if (Simu_Current_Step==4){
    //	 exit(0);
   //	}

    t_vmap= timer->createCategory("DS Solve: Map Volt to Bus");
    timer->start(t_vmap);
	
	//printf("after first volt sovle, before first volt map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap_sptr->mapToBus(volt_full);
	
	//printf("after first volt sovle, after first volt map: \n");
	
	if ( Simu_Current_Step==0 ) {
		//printf("enter the initial update oldbusvoltage, Timestep: %d \n", Simu_Current_Step);
		p_factory->updateoldbusvoltage(); //renke add, first timestep, copy volt_full to volt_full_old
	}
    timer->stop(t_vmap);

    t_volt= timer->createCategory("DS Solve: Set Volt");
    timer->start(t_volt);
    p_factory->setVolt(false);
	p_factory->updateBusFreq(h_sol1);
	
	
	std::vector <double> vwideareafreqs;
	vwideareafreqs = p_factory->grabWideAreaFreq();
	//printf("-----!!renke debug dsf_app_module.cpp: grabWideAreaFreq: bus 30: %12.6f, bus 30: %12.6f, delta_freq bus34-bus30: %12.6f \n", 
	//		vwideareafreqs[0], vwideareafreqs[1], vwideareafreqs[2]);
	int tmp = vwideareafreqs.size();
	double widearea_deltafreq = vwideareafreqs[tmp-1];

#ifdef USE_HELICS
	 
	 //pub.publish(widearea_deltafreq);
	 for(int i = 0; i < pubCount; i++) {
            pub = fed.getPublication(i);
            string pubInfo = pub.getInfo();
			//std::cout << "-------------!!!helics test: HELICS pub info: " << pubInfo << std::endl;
			pub.publish(vwideareafreqs[i]);
            // do stuff to tie pub to GridPACK object property
          }

	 helics_requestTime =       double (Simu_Current_Step*h_sol1);
	 //printf("-------------!!!Helics request time: %12.6f \n", helics_requestTime); 
	 double helics_grantime;
	 helics_grantime = fed.requestTime(helics_requestTime);
	 //printf("-------------!!!Helics grant time: %12.6f \n", helics_grantime); 
	 
	 double subvalue = 0.0;
	 
	 for(int i = 0; i < subCount; i++) {
        sub = fed.getInput(i);
		//printf("-------------!!!helics debug entering  sub loop\n"); 
		//if(sub.isUpdated()) {
            //auto value = sub.getValue();
			subvalue = fed.getDouble(sub);
			//printf("-------------!!!Helics sub value: %12.6f \n", subvalue);
                             //update GridPACK object property with value
        //}

	 }
	 
	//printf("-------------!!!Outside Helics def sub value: %12.6f \n", subvalue);
	 
	p_factory->setWideAreaFreqforPSS(subvalue);

	//p_factory->setWideAreaFreqforPSS(widearea_deltafreq);
	 
#else	 
	 
	p_factory->setWideAreaFreqforPSS(widearea_deltafreq);
	 
#endif
		
    timer->stop(t_volt);
	
	//printf("before update relay, after first volt solv: \n");
	//p_factory->printallbusvoltage();
    //renke add, compute bus freq if necessary
    //printf("Timestep, %d \n", Simu_Current_Step);
    bool flagBus = p_factory->updateBusRelay(false, h_sol1);
    bool flagBranch = p_factory->updateBranchRelay(false, h_sol1);
	
	// update dynamic load internal relay functions here
	p_factory->dynamicload_post_process(h_sol1, false);
    
	// if bus relay trips, modify the corresponding Ymatrix, renke modified
    if (flagBus) {
        printf("DSFull_APP::Solve: updatebusrelay return trigger siganl: TURE!!! \n");
		
        //please update the bus contribution to the Y bus matrix here. //Shuangshuang tbd
	if (flagP == 0) { 
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",Simu_Current_Step );
			 
	     ybus_fy->save(sybus);
	 
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
	     ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

	     ybus->save(sybus);
        
             printf("DSFull_APP::Solve: bus relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

			 
        } else if (flagP == 2) {
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
             printf("DSFull_APP::Solve: bus relay trip after fault, ybus changed:\n");
	     ybus->print();
	     char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

	     ybus->save(sybus);

             printf("DSFull_APP::Solve: bus relay trip after fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

        }
    }
	
	// if branch relay trips, modify the corresponding Ymatrix, renke modified
	if (flagBranch) {
        
        printf("DSFull_APP::Solve: updatebranchrelay return trigger siganl: TURE!!! \n");

        //please update the bus contribution to the Y bus matrix here. //Shuangshuang tbd
	if (flagP == 0) { 
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",Simu_Current_Step );
			 
	     ybus_fy->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
             ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

             ybus->save(sybus);

			 
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_posfy);
	     ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

	     ybus_posfy->save(sybus);
			 
        } else if (flagP == 2) {
             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
             ybus->print();
             char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

             ybus->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

        }
    }
	
    //renke add, update old busvoltage first
    p_factory->updateoldbusvoltage(); //renke add
	
	//printf("after updateoldbus voltage: \n");
	//p_factory->printallbusvoltage();
	
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif

    t_predictor = timer->createCategory("DS Solve: Modified Euler Predictor");
    //printf("Test: predictor begins: \n");
    timer->start(t_predictor);
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor(h_sol1, false);
    } else { 
      p_factory->predictor(h_sol1, true);
    }
    timer->stop(t_predictor);

    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->corrector_currentInjection(false);
    } else {
      p_factory->corrector_currentInjection(true);
    }

    //INorton_full = nbusMap_sptr->mapToVector();
    int t_cmIf = timer->createCategory("DS Solve: Modified Euler Corrector: Make INorton");
    timer->start(t_cmIf);
    p_factory->setMode(make_INorton_full);
    nbusMap_sptr->mapToVector(INorton_full);
    //p_busIO->header("\n=== [Corrector] INorton_full: ===\n");
    //printf("\nrelaytest=== [Corrector] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_cmIf);

    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2)
    // to calculate terminal volt: ----------
    t_csolve = timer->createCategory("DS Solve: Modified Euler Corrector: Linear Solver");
    timer->start(t_csolve);
    volt_full->zero();

#if 0
    flag_chk = true;
    while (flag_chk == true ) {
		
			volt_full->zero();
			
			if (flagP == 0) {
				solver_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 1) {
				solver_fy_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 2) {
				solver_posfy_sptr->solve(*INorton_full, *volt_full);
			}
			nbusMap_sptr->mapToBus(volt_full);
			p_factory->setVolt(false);
			
			if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
				p_factory->corrector_currentInjection(false);
			} else {
				p_factory->corrector_currentInjection(true);
			}
			
			INorton_full_chk->equate(*INorton_full);
			printf("itr test:----corrector_INorton_full_chk:\n");
			INorton_full_chk->print();
			
			p_factory->setMode(make_INorton_full);
			nbusMap_sptr->mapToVector(INorton_full);
			
			printf("itr test:----corrector_INorton_full:\n");
			INorton_full->print();
			
			//multiply(*ybus_fy, *volt_full, *INorton_full_chk);
			INorton_full_chk->add(*INorton_full, -1.0);
			max_INorton_full=abs(INorton_full_chk->normInfinity());
			
			if (max_INorton_full <1.0e-8) {
				flag_chk = false;
			} else {
				printf("max_INorton_full = %8.4f \n", max_INorton_full);
				//printf("-----INorton_full : \n");
				//INorton_full->print();
				//printf("-----INorton_full_chk - INorton_full : \n");
				//INorton_full_chk->print();
			}
    }
#else
    if (flagP == 0) {
      solver_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy_sptr->solve(*INorton_full, *volt_full);
    }
#endif

    timer->stop(t_csolve);

    //p_busIO->header("\n=== [Corrector] volt_full: ===\n");
    //printf("relaytest \n=== [Corrector] volt_full: ===\n");
    //volt_full->print();
    timer->start(t_vmap);
	
	//printf("after second solve, before second map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap_sptr->mapToBus(volt_full);
	
	//printf("after second solve, after second map: \n");
	//p_factory->printallbusvoltage();
	
    timer->stop(t_vmap);

    timer->start(t_volt);
    p_factory->setVolt(false);
	p_factory->updateBusFreq(h_sol1);
    timer->stop(t_volt);

    t_corrector = timer->createCategory("DS Solve: Modified Euler Corrector");
    timer->start(t_corrector);
    //printf("Test: corrector begins: \n");
    if (last_S_Steps != S_Steps) {
      p_factory->corrector(h_sol2, false);
    } else {
      p_factory->corrector(h_sol2, true);
    }
    timer->stop(t_corrector);

    //if (Simu_Current_Step == simu_total_steps - 1) 
      //p_busIO->write();

    if (Simu_Current_Step == steps1) {
      solver_fy_sptr->solve(*INorton_full, *volt_full);
//      printf("\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
//      printf("\n=== [Corrector] volt_full: ===\n");
//      volt_full->print();
      nbusMap_sptr->mapToBus(volt_full);
      p_factory->setVolt(false);
	  p_factory->updateBusFreq(h_sol1);
    } else if (Simu_Current_Step == steps2) {
      solver_posfy_sptr->solve(*INorton_full, *volt_full);
//      printf("\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
//      printf("\n=== [Corrector] volt_full: ===\n");
//      volt_full->print();
      nbusMap_sptr->mapToBus(volt_full);
      p_factory->setVolt(true);
	  p_factory->updateBusFreq(h_sol1);
    }
    if (Simu_Current_Step == 1) {
//      printf("\n Dynamic Step 1 [Corrector] volt_full: ===\n");
//      volt_full->print();
//      printf("\n Dynamic Step 1 [Corrector] Norton_full: ===\n");
//      INorton_full->print();
    }
    t_secure = timer->createCategory("DS Solve: Check Security");
    timer->start(t_secure);
    if (p_generatorWatch && Simu_Current_Step%p_generatorWatchFrequency == 0) {
      char tbuf[32];
      if (!p_suppress_watch_files) {
#ifdef USE_TIMESTAMP
        sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(Simu_Current_Step)*p_time_step,
            timer->currentTime());
        if (p_generatorWatch) p_generatorIO->header(tbuf);
        if (p_generatorWatch) p_generatorIO->write("watch");
        if (p_generatorWatch) p_generatorIO->header("\n");

        //      if (p_generatorWatch) p_generatorIO->write("watch");

        //      sprintf(tbuf,"%8.4f, %20.4f",mac_ang_s0, mac_spd_s0);
        //      if (p_generatorWatch) p_generatorIO->header(tbuf);
        //      if (p_generatorWatch) p_generatorIO->write("watch");
        //      if (p_generatorWatch) p_generatorIO->header("\n");
#else
        sprintf(tbuf,"%8.4f",static_cast<double>(Simu_Current_Step)*p_time_step);
        if (p_generatorWatch) p_generatorIO->header(tbuf);
        if (p_generatorWatch) p_generatorIO->write("watch");
        if (p_generatorWatch) p_generatorIO->header("\n");
#endif
      }
#ifdef USEX_GOSS
      if (p_generatorWatch) p_generatorIO->dumpChannel();
#endif
    }
    if (p_loadWatch && Simu_Current_Step%p_loadWatchFrequency == 0) {
      char tbuf[32];
#ifdef USE_TIMESTAMP
      sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(Simu_Current_Step)*p_time_step,
          timer->currentTime());
      if (p_loadWatch) p_loadIO->header(tbuf);
      if (p_loadWatch) p_loadIO->write("load_watch");
      if (p_loadWatch) p_loadIO->header("\n");
#else
      sprintf(tbuf,"%8.4f",static_cast<double>(Simu_Current_Step)*p_time_step);
      if (p_loadWatch) p_loadIO->header(tbuf);
      if (p_loadWatch) p_loadIO->write("load_watch");
      if (p_loadWatch) p_loadIO->header("\n");
#endif
#ifdef USEX_GOSS
      if (p_loadWatch) p_loadIO->dumpChannel();
#endif
    }
    saveTimeStep();
    //    if ((!p_factory->securityCheck()) && p_insecureAt == -1)  
    //       p_insecureAt = Simu_Current_Step;

    
/*    // Print to screen
    if (last_S_Steps != S_Steps) {
      //sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps);
      //p_busIO->header(ioBuf);
      printf("\n==============S_Steps = %d==============\n", S_Steps);
      mac_ang_s0->print();
      mac_spd_s0->print();
      //pmech->print();
      //pelect->print();
      //sprintf(ioBuf, "========================End of S_Steps = %d=========================\n\n", S_Steps);
      //p_busIO->header(ioBuf);
    }
    if (Simu_Current_Step == simu_total_steps) {
      printf("\n==============S_Steps = %d==============\n", S_Steps);
      mac_ang_s1->print();
      mac_spd_s1->print();
      p_factory->setMode(init_mac_ang);
      ngenMap_sptr->mapToBus(mac_ang_s1);
      p_factory->setMode(init_mac_spd);
      ngenMap_sptr->mapToBus(mac_spd_s1);
      p_factory->setMode(init_pmech);
      ngenMap_sptr->mapToBus(pmech);
      p_factory->setMode(init_pelect);
      ngenMap_sptr->mapToBus(pelect);
      sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps+1);
      p_busIO->header(ioBuf);
      sprintf(ioBuf, "\n         Bus ID     Generator ID"
          "    mac_ang         mac_spd         mech            elect\n\n");
      p_busIO->header(ioBuf);
      mac_ang_s1->print();
      mac_spd_s1->print();
      pmech->print();
      pelect->print();
      p_busIO->write();
      sprintf(ioBuf, "\n========================End of S_Steps = %d=========================\n\n", S_Steps+1);
      p_busIO->header(ioBuf);
    } // End of Print to screen

*/    //exit(0);
    last_S_Steps = S_Steps;
    timer->stop(t_secure);
    if (p_monitorGenerators) {
      double presentTime = static_cast<double>(Simu_Current_Step)*p_time_step;
      p_frequencyOK = p_frequencyOK && checkFrequency(0.5,presentTime);
      if (!p_frequencyOK) Simu_Current_Step = simu_total_steps;
    }
  }
  
#if 0
  printf("\n=== ybus after simu: ============\n");
  ybus->print();
  ybus->save("ybus_aftersimu.m");
  
  printf("\n=== ybus_fy after simu:============\n");
  ybus_fy->print();
  ybus_fy->save("ybus_fy_aftersimu.m");
  
  printf("\n=== ybus_posfy after simu: ============\n");
  ybus_posfy->print();
  ybus_posfy->save("ybus_posfy_aftersimu.m");
  
#endif

  //char msg[128];
  //if (p_insecureAt == -1) sprintf(msg, "\nThe system is secure!\n");
  //else sprintf(msg, "\nThe system is insecure from step %d!\n", p_insecureAt);

  char secureBuf[128];
  if (p_insecureAt == -1) {
    char *ptr;
    sprintf(secureBuf,"\nThe system is secure");
    ptr = secureBuf + strlen(secureBuf);
    if (fault.isBus) {
      sprintf(ptr," for fault at bus %d\n", fault.bus_idx);
    } else if (fault.isLine) {
      sprintf(ptr," for fault at line %s from bus %d to bus %d\n",fault.tag.c_str(),
          fault.from_idx,fault.to_idx);
    } else {
      sprintf(ptr,"!\n");
    }
  } else { 
    char *ptr;
    sprintf(secureBuf,"\nThe system is insecure from step %d", p_insecureAt);
    ptr = secureBuf + strlen(secureBuf);
    if (fault.isBus) {
      sprintf(ptr," for fault on bus %d\n",fault.bus_idx);
    } else if (fault.isLine) {
      sprintf(ptr," for fault at line %s from bus %d to bus %d\n",fault.tag.c_str(),
          fault.from_idx,fault.to_idx);
    } else {
      sprintf(ptr,"!\n");
    }
  }
  p_busIO->header(secureBuf);

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
  timer->stop(t_solve);
  //timer->dump();
  
#ifdef USE_HELICS

	fed.finalize();
	
#endif
 
  
}

/**
 * Write out final results of dynamic simulation calculation to
 * standard output
 */
void gridpack::dynamic_simulation::DSFullApp::write(const char* signal)
{
  p_busIO->write(signal);
}

/**
 * Read in generators that should be monitored during simulation
   with a cursor ptr given
*/
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch(gridpack::utility::Configuration::CursorPtr cursor)
{
  bool noprint = gridpack::NoPrint::instance()->status();
  gridpack::utility::Configuration::ChildCursors generators;
  if (cursor) cursor->children(generators);
  int i, j, idx, id, len;
  int ncnt = generators.size();
  std::string generator, tag, clean_tag;
  gridpack::dynamic_simulation::DSFullBus *bus;
  if (!noprint) {
	if (ncnt > 0) p_busIO->header("Monitoring generators:\n");
  }
  std::vector<int> buses;
  std::vector<std::string> tags;
  for (i=0; i<ncnt; i++) {
    // Parse contents of "generator" field to get bus ID and generator tag
    generators[i]->get("busID",&id);
    generators[i]->get("generatorID",&tag);
    gridpack::utility::StringUtils util;
    clean_tag = util.clean2Char(tag);
    buses.push_back(id);
    tags.push_back(clean_tag);
  }
  setGeneratorWatch(buses,tags,true);
}

/**
 * Read in generators that should be monitored during simulation
 */
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  if (!cursor->get("generatorWatchFrequency",&p_generatorWatchFrequency)) {
    p_generatorWatchFrequency = 1;
  }
  cursor = p_config->getCursor("Configuration.Dynamic_simulation.generatorWatch");
  setGeneratorWatch(cursor);  
}

/**
 * Read in generators that should be monitored during simulation
 * @param filename set filename from calling program instead of input
 *        deck
 */
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch(const char *filename)
{
  p_gen_watch_file = filename;
  p_internal_watch_file_name = true;
  setGeneratorWatch();
}

/**
 * Read in generators that should be monitored during simulation
 * @param filename set filename from calling program instead of input
 *        deck
 */
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch(const char *filename,gridpack::utility::Configuration::CursorPtr cursor)
{
  p_gen_watch_file = filename;
  p_internal_watch_file_name = true;
  setGeneratorWatch(cursor);
}


/**
 * Read in generators that should be monitored during simulation
 * @param buses IDs of buses containing generators
 * @param tags generator IDs for watched generators
 * @param writeFile true if external file is to be written
 */
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch(
    std::vector<int> &buses, std::vector<std::string> &tags, bool writeFile)
{
  int ncnt = buses.size();
  bool noprint = gridpack::NoPrint::instance()->status();														 
  if (ncnt != tags.size()) {
    printf("setGeneratorWatch: size mismatch between buses: and tags: vectors\n",
        (int)buses.size(),(int)tags.size());
    // TODO: some kind of error
  }
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  if (!cursor->get("generatorWatchFrequency",&p_generatorWatchFrequency)) {
    p_generatorWatchFrequency = 1;
  }
  std::string generator, tag;
  char buf[128];
  gridpack::dynamic_simulation::DSFullBus *bus;
  p_watch_bus_ids.clear();
  p_watch_gen_ids.clear();
  p_gen_buses.clear();
  p_gen_ids.clear();
  int i, j, id;
  for (i=0; i<ncnt; i++) {
    id = buses[i];
    tag = tags[i];
    std::pair<int,std::string> list_item = std::pair<int,std::string>(id,tag);
    p_watch_list.insert(std::pair<std::pair<int,std::string>, int>(list_item,i));
    p_watch_bus_ids.push_back(id);
    p_watch_gen_ids.push_back(tag);
    // Find local bus indices for generator. If generator is not on this
    // processor then local_ids will have zero length.
    std::vector<int> local_ids = p_network->getLocalBusIndices(id);
    for (j=0; j<local_ids.size(); j++) {
      bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(local_ids[j]).get());
      bus->setWatch(tag,true);
      if (p_network->getActiveBus(local_ids[j])) {
        p_gen_buses.push_back(local_ids[j]);
        p_gen_ids.push_back(tag);
      }
    }
	if (!noprint) {				
		sprintf(buf,"  Bus: %8d Generator ID: %2s\n",id,tag.c_str());
		p_busIO->header(buf);
	}
    if (ncnt > 0) {
      p_generators_read_in = true;
      p_generatorWatch = true;
	  if (!noprint) {				  
		sprintf(buf,"Generator Watch Frequency: %d\n",p_generatorWatchFrequency);
		p_busIO->header(buf);
	  }
    }
  }

  // If storing time series data, set up vector to hold results
  if (writeFile) {
    openGeneratorWatchFile();
    p_monitorGenerators = false;
  } else {
    p_generatorWatch = false;
    p_monitorGenerators = true;
  }
  if (p_save_time_series) {
    p_time_series.clear();
    for (i=0; i<p_gen_buses.size(); i++) {
      std::vector<double> vec0;
      p_time_series.push_back(vec0);
      p_time_series.push_back(vec0);
      p_time_series.push_back(vec0);
      p_time_series.push_back(vec0);
    }
  }
}

/**
 * Check to see if frequency variations on monitored generators are okay
 * @param start time at which to start monitoring
 * @param time current value of time
 * @return true if all watched generators are within acceptable bounds
 */
bool gridpack::dynamic_simulation::DSFullApp::checkFrequency(
    double start, double time)
{
  int nbus = p_network->numBuses();
  int i;
  bool ret = true;
  bool ok = true;
  p_violations.clear();
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(i)) {
      ok = p_network->getBus(i)->checkFrequency(start,time);
      if (!ok) {
        p_violations.push_back(p_network->getBus(i)->getOriginalIndex());
      }
      ret = ret && ok;
    }
  }
  return p_factory->checkTrue(ret);
}

/**
 * Set parameters for monitoring frequency
 * @param flag true if frequency monitoring is turned on
 * @param maxFreq maximum allowable frequency deviation
 */
void gridpack::dynamic_simulation::DSFullApp::setFrequencyMonitoring(
    bool flag, double maxFreq)
{
  p_monitorGenerators = flag;
  p_maximumFrequency = maxFreq;
}

/**
 * Check to see if frequency variations on monitored generators are okay
 * @param limit maximum upper limit on frequency deviation
 * @return true if all watched generators are within acceptable bounds
 */
bool gridpack::dynamic_simulation::DSFullApp::checkFrequency(double limit)
{
  int nbus = p_network->numBuses();
  int i;
  bool ret = true;
  bool ok = true;
  p_violations.clear();
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(i)) {
      ok = p_network->getBus(i)->checkFrequency(limit);
      if (!ok) {
        p_violations.push_back(p_network->getBus(i)->getOriginalIndex());
      }
      ret = ret && ok;
    }
  }
  return p_factory->checkTrue(ret);
}

/**
 * @return true if no frequency violations occured on monitored generators
 */
bool gridpack::dynamic_simulation::DSFullApp::frequencyOK()
{
  return p_frequencyOK;
}

/**
 * Scale generator real power. If zone less than 1 then scale all
 * generators in the area.
 * @param scale factor to scale real power generation
 * @param area index of area for scaling generation
 * @param zone index of zone for scaling generation
 */
void gridpack::dynamic_simulation::DSFullApp::scaleGeneratorRealPower(
    double scale, int area, int zone)
{
  p_factory->scaleGeneratorRealPower(scale,area,zone);
}

/**
 * Scale load power. If zone less than 1 then scale all
 * loads in the area.
 * @param scale factor to scale load real power
 * @param area index of area for scaling load
 * @param zone index of zone for scaling load
 */
void gridpack::dynamic_simulation::DSFullApp::scaleLoadPower(
    double scale, int area, int zone)
{
  return p_factory->scaleLoadPower(scale,area,zone);
}

/**
 * Return the total real power load for all loads in the zone. If zone
 * less than 1, then return the total load for the area
 * @param area index of area
 * @param zone index of zone
 * @return total load
 */
double gridpack::dynamic_simulation::DSFullApp::getTotalLoadRealPower(int area,
    int zone)
{
  return p_factory->getTotalLoadRealPower(area,zone);
}

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
void gridpack::dynamic_simulation::DSFullApp::getGeneratorMargins(int area,
    int zone, double *total, double *pmin, double *pmax)
{
  p_factory->getGeneratorMargins(area,zone,total,pmin,pmax);
}

/**
 * Reset power of loads and generators to original values
 */
void gridpack::dynamic_simulation::DSFullApp::resetPower()
{
  return p_factory->resetPower();
}

/**
 * Read in loads that should be monitored during simulation
 */
void gridpack::dynamic_simulation::DSFullApp::setLoadWatch()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  if (!cursor->get("loadWatchFrequency",&p_loadWatchFrequency)) {
    p_loadWatchFrequency = 1;
  }
  char buf[128];
  cursor = p_config->getCursor("Configuration.Dynamic_simulation.loadWatch");
  gridpack::utility::Configuration::ChildCursors loads;
  if (cursor) cursor->children(loads);
  int i, j, idx, id, len;
  int ncnt = loads.size();
  std::string load, tag, clean_tag;
  gridpack::dynamic_simulation::DSFullBus *bus;
  if (ncnt > 0) p_busIO->header("Monitoring loads:\n");
  for (i=0; i<ncnt; i++) {
    // Parse contents of "load" field to get bus ID and load tag
    loads[i]->get("busID",&id);
    loads[i]->get("loadID",&tag);
    gridpack::utility::StringUtils util;
    clean_tag = util.clean2Char(tag);
    // Find local bus indices for load
    std::vector<int> local_ids = p_network->getLocalBusIndices(id);
    for (j=0; j<local_ids.size(); j++) {
      bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(local_ids[j]).get());
      bus->setWatch(clean_tag,true);
    }
    sprintf(buf,"  Bus: %8d Load ID: %2s\n",id,clean_tag.c_str());
    p_busIO->header(buf);
  }
  if (ncnt > 0) {
    p_loadWatch = true;
    sprintf(buf,"Load Watch Frequency: %d\n",p_loadWatchFrequency);
    p_busIO->header(buf);
    openLoadWatchFile();
  }
}

/**
 * Save watch series to an internal data vector
 * @param flag if true, save time series data
 */
void gridpack::dynamic_simulation::DSFullApp::saveTimeSeries(bool flag)
{
  p_save_time_series = flag;
}

/**
 * Save time series data for watched generators
 */
void gridpack::dynamic_simulation::DSFullApp::saveTimeStep()
{
  if (!p_save_time_series) return;
  int nbus = p_gen_buses.size();
  int i, j;
  int icnt = 0;
  gridpack::dynamic_simulation::DSFullBus *bus;
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(p_gen_buses[i])) {
      bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(p_gen_buses[i]).get());
      std::vector<double> vals = bus->getWatchedValues(p_gen_ids[i]);
      for (j=0; j<vals.size(); j++) {
        p_time_series[icnt].push_back(vals[j]);
        icnt++;
      }
    }
  }
}

/**
 * Return global map of time series values
 * @return map of time series indices (local to global)
 */
std::vector<int> gridpack::dynamic_simulation::DSFullApp::getTimeSeriesMap()
{
  std::vector<int> ret;
  if (p_save_time_series) {
    std::vector<int> orig_idx;
    std::vector<std::string> tags;
    int nbus = p_network->numBuses();
    int i, j;
    int icnt = 0;
    gridpack::dynamic_simulation::DSFullBus *bus;
    for (i=0; i<nbus; i++) {
      if (p_network->getActiveBus(i)) {
        bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
          (p_network->getBus(i).get());
        std::vector<std::string> watched = bus->getWatchedGenerators();
        for (j=0; j<watched.size(); j++) {
          std::pair<int,std::string> gen
            = std::pair<int,std::string>(bus->getOriginalIndex(),watched[j]);
          std::map<std::pair<int,std::string>,int>::iterator it;
          it = p_watch_list.find(gen);
          if (it != p_watch_list.end()) {
            ret.push_back(4*(it->second));
            ret.push_back(4*(it->second)+1);
            ret.push_back(4*(it->second)+2);
            ret.push_back(4*(it->second)+3);
          } else {
            printf("Could not find generator %s on bus %d\n",
                watched[j].c_str(),bus->getOriginalIndex());
          }
          icnt++;
        }
      }
    }
  }
  return ret;
}

/**
 * Return a list of original bus IDs and tags for all monitored
 * generators
 * @param bus_ids list of original bus indices for monitored generators
 * @param gen_ids list of tags for monitored generators
 */
void gridpack::dynamic_simulation::DSFullApp::getListWatchedGenerators(
    std::vector<int> &bus_ids, std::vector<std::string> &gen_ids)
{
  bus_ids.clear();
  gen_ids.clear();
  int nsize = p_watch_bus_ids.size();
  int i;
  for (i=0; i<nsize; i++) {
    bus_ids.push_back(p_watch_bus_ids[i]);
    gen_ids.push_back(p_watch_gen_ids[i]);
  }
}

/**
 * Return vector of time series data for watched generators
 * @return vector of time series for generators on this
 processor
 */
std::vector<std::vector<double> > gridpack::dynamic_simulation::DSFullApp::getGeneratorTimeSeries()
{
  std::vector<std::vector<double> > ret;
  if (p_save_time_series) {
    int ngen = p_time_series.size();
    int i, j;
    for (i=0; i<ngen; i++) {
      std::vector<double> series;
      int nsteps = (p_time_series[i]).size();
      for (j=0; j<nsteps; j++) {
        series.push_back((p_time_series[i])[j]);
      }
      ret.push_back(series);
    }
  }
  return ret;
}

/**
 * Redirect output from standard out
 * @param filename name of file to write results to
 */
void gridpack::dynamic_simulation::DSFullApp::open(const char *filename)
{
  printf("open busIO (%s)\n",filename);
  p_busIO->open(filename);
  printf("open branchIO\n");
  p_branchIO->setStream(p_busIO->getStream());
  printf("finished open\n");
}

void gridpack::dynamic_simulation::DSFullApp::close()
{
  printf("close busIO\n");
  p_busIO->close();
  printf("close branchIO\n");
  p_branchIO->setStream(p_busIO->getStream());
  printf("finished close\n");
}

/**
 * Print string. This can be used to direct output to the file opened using
 * the open command
 * @param buf string to be printed
 */
void gridpack::dynamic_simulation::DSFullApp::print(const char *buf)
{
    p_busIO->header(buf);
}

/**
 * Open file containing generator watch results
 */
void gridpack::dynamic_simulation::DSFullApp::openGeneratorWatchFile()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
#ifndef USEX_GOSS
  std::string filename;
  std::string flag;
  cursor->get("suppressWatchFiles", &flag);
  gridpack::utility::StringUtils util;
  p_suppress_watch_files = util.getBool(flag.c_str());
  if (!p_internal_watch_file_name) {
    if (!p_suppress_watch_files) {
      if (cursor->get("generatorWatchFileName",&filename)) {
        p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(512,
              p_network));
        p_generatorIO->open(filename.c_str());
      } else {
		  bool noprint = gridpack::NoPrint::instance()->status();
		  if (!noprint) {														   	   
			p_busIO->header("No Generator Watch File Name Found\n");
		  }
        p_generatorWatch = false;
      }
    }
  } else {
    if (!p_suppress_watch_files) {
      p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(512,
            p_network));
      p_generatorIO->open(p_gen_watch_file.c_str());
    }
  }
#else
  std::string topic, URI, username, passwd;
  bool ok = true;
  ok = ok && cursor->get("channelTopic",&topic);
  ok = ok && cursor->get("channelURI",&URI);
  ok = ok && cursor->get("username",&username);
  ok = ok && cursor->get("password",&passwd);
  if (p_internal_watch_file_name) {
    topic = p_gen_watch_file;
  }
  printf("channeltopic %s \n", topic.c_str());
  printf("channelURI %s \n", URI.c_str());
  printf("username %s \n", username.c_str());
  printf("password %s \n", passwd.c_str());
  if (ok) {
    p_generatorIO.reset(new
        gridpack::serial_io::SerialBusIO<DSFullNetwork>(512,
          p_network));
    p_generatorIO->openChannel(topic.c_str());
  } else {
    p_busIO->header("Unable to open channel\n");
    p_generatorWatch = false;
  }
#endif
}

/**
 * Close file containing generator watch results
 */
void gridpack::dynamic_simulation::DSFullApp::closeGeneratorWatchFile()
{
  if (p_generatorWatch) {
#ifndef USEX_GOSS
    if (!p_suppress_watch_files) {
      p_generatorIO->close();
    }
#else
    p_generatorIO->closeChannel();
#endif
  }
}

/**
 * Open file containing load watch results
 */
void gridpack::dynamic_simulation::DSFullApp::openLoadWatchFile()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
#ifndef USEX_GOSS
  std::string filename;
  std::string flag;
  cursor->get("suppressWatchFiles", &flag);
  gridpack::utility::StringUtils util;
  p_suppress_watch_files = util.getBool(flag.c_str());
  if (!p_suppress_watch_files) {
    if (cursor->get("loadWatchFileName",&filename)) {
      p_loadIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(128,
            p_network));
      p_loadIO->open(filename.c_str());
    } else {
      p_busIO->header("No Load Watch File Name Found\n");
      p_loadWatch = false;
    }
  }
#else
  std::string topic, URI, username, passwd;
  bool ok = true;
  ok = ok && cursor->get("channelTopic",&topic);
  ok = ok && cursor->get("channelURI",&URI);
  ok = ok && cursor->get("username",&username);
  ok = ok && cursor->get("password",&passwd);
  printf("channeltopic %s \n", topic.c_str());
  printf("channelURI %s \n", URI.c_str());
  printf("username %s \n", username.c_str());
  printf("password %s \n", passwd.c_str());
  if (ok) {
    p_loadIO.reset(new
        gridpack::serial_io::SerialBusIO<DSFullNetwork>(512,
          p_network));
    p_loadIO->openChannel(topic.c_str());
  } else {
    p_busIO->header("Unable to open channel\n");
    p_loadWatch = false;
  }
#endif
}

/**
 * Close file contain load watch results
 */
void gridpack::dynamic_simulation::DSFullApp::closeLoadWatchFile()
{
  if (p_loadWatch) {
#ifndef USEX_GOSS
    if (!p_suppress_watch_files) {
      p_loadIO->close();
    }
#else
    p_loadIO->closeChannel();
#endif
  }
}

/**
 * Get observations and store them internally
 * @param cursor configuration pointer to observation block
 */
void gridpack::dynamic_simulation::DSFullApp::setObservations(
    gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("observations");
  gridpack::utility::Configuration::ChildCursors observations;
  p_obs_genBus.clear();
  p_obs_genIDs.clear();
  p_obs_vBus.clear();
  p_obs_vBusfreq.clear();
  gridpack::utility::StringUtils util;
  // Parser observation block
  std::vector<int> foundGen;
  std::vector<int> foundBus;
  std::vector<int> foundBusfreq;
  std::vector<int> foundLoad;
  int idx;
  if (list) {
    list->children(observations);
    int size = observations.size();
    // Find out which observations actually correspond to elements that exist
    for (idx=0; idx<size; idx++) {
      std::string type;
      if (!observations[idx]->get("type",&type)) continue;
      util.trim(type);
      util.toLower(type);
      if (type == "generator") {
        int bus;
        std::string genID, tID;
        if (observations[idx]->get("busID",&bus) &&
            observations[idx]->get("generatorID",&genID)) {
          tID = util.clean2Char(genID);
          p_obs_genBus.push_back(bus);
          p_obs_genIDs.push_back(tID);
          foundGen.push_back(0);
        }
      } else if (type == "bus") {
        int bus;
        if (observations[idx]->get("busID",&bus)) {
          p_obs_vBus.push_back(bus);
          foundBus.push_back(0);
        }
      }else if (type == "busfrequency") {
        int bus;
        if (observations[idx]->get("busID",&bus)) {
          p_obs_vBusfreq.push_back(bus);
          foundBusfreq.push_back(0);
        }
      }
	  else if (type == "load") {
        int bus;
        std::string loadID, tID;
        if (observations[idx]->get("busID",&bus) &&
            observations[idx]->get("loadID",&loadID)) {
          tID = util.clean2Char(loadID);
          p_obs_loadBus.push_back(bus);
          p_obs_loadIDs.push_back(tID);
          foundLoad.push_back(0);
        }
      } else {
        printf("Unknown observation type: %s\n",type.c_str());
      }
    }
  }
  // create global vectors to store values of observations. Check to see if
  // any observations are on this processor
  p_obs_vMag.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_vAng.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_vBusfreqVal.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_rSpd.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_rAng.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_genP.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_genQ.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  p_obs_fOnline.reset(new gridpack::parallel::GlobalVector<double>(p_comm));
  if (p_obs_genBus.size() > 0) {
    int nbus = p_obs_genBus.size();
    p_obs_gActive.resize(nbus);
    p_obs_lGenBus.clear();
    p_obs_lGenIDs.clear();
    p_obs_lGenIdx.clear();
    p_obs_gUse.clear();
    std::vector<double> zeros;
    std::vector<double> ones;
    int i, j, k, lidx;
    for (i = 0; i<nbus; i++) {
      std::vector<int> localIndices;
      localIndices = p_network->getLocalBusIndices(p_obs_genBus[i]);
      bool isLocal = false;
      p_obs_gActive[i] = 0;
      // Check to see if generator host is active on this processor
      for (j=0; j<localIndices.size(); j++) {
        if (p_network->getActiveBus(localIndices[j])) {
          // Check to see if generator is on this bus
          std::vector<std::string> tags
            = p_network->getBus(localIndices[j])->getGenerators();
          for (k = 0; k<tags.size(); k++) {
            if (tags[k] == p_obs_genIDs[i]) {
              lidx = localIndices[j];
              isLocal = true;
              p_obs_gActive[i] = 1;
              break;
            }
          }
          if (isLocal) break;
        }
      }
      if (isLocal) {
        foundGen[i] = 1;
        p_obs_lGenIdx.push_back(i);
        p_obs_GenIdx.push_back(lidx);
        p_obs_lGenBus.push_back(p_obs_genBus[i]);
        p_obs_lGenIDs.push_back(p_obs_genIDs[i]);
        p_obs_gUse.push_back(0);
        zeros.push_back(0.0);
        ones.push_back(1.0);
      }
    }
    // Sum foundGen vector over all processors to find out if any observations
    // do not correspond to existing elements
    p_network->communicator().sum(&foundGen[0],foundGen.size());
    // Add observations not found on any processors to observations on process 0
    if (p_comm.rank() == 0 && p_report_dummy_obs) {
      for (i = 0; i<foundGen.size(); i++) {
        if (!foundGen[i]) {
          p_obs_lGenIdx.push_back(i);
          p_obs_GenIdx.push_back(-1);
          p_obs_lGenBus.push_back(p_obs_genBus[i]);
          p_obs_lGenIDs.push_back(p_obs_genIDs[i]);
          p_obs_gUse.push_back(1);
          zeros.push_back(0.0);
          ones.push_back(1.0);
        }
      }
    }
    p_obs_rSpd->addElements(p_obs_lGenIdx, zeros);
    p_obs_rSpd->upload();
    p_obs_rAng->addElements(p_obs_lGenIdx, zeros);
    p_obs_rAng->upload();
    p_obs_genP->addElements(p_obs_lGenIdx, zeros);
    p_obs_genP->upload();
    p_obs_genQ->addElements(p_obs_lGenIdx, zeros);
    p_obs_genQ->upload();
    p_comm.sum(&p_obs_gActive[0],p_comm.size());
  }
  if (p_obs_loadBus.size() > 0) {
    int nbus = p_obs_loadBus.size();
    p_obs_lActive.resize(nbus);
    p_obs_lLoadBus.clear();
    p_obs_lLoadIDs.clear();
    p_obs_lLoadIdx.clear();
    p_obs_lUse.clear();
    std::vector<double> ones;
    int i, j, k, lidx;
    for (i = 0; i<nbus; i++) {
      std::vector<int> localIndices;
      localIndices = p_network->getLocalBusIndices(p_obs_loadBus[i]);
      bool isLocal = false;
      p_obs_lActive[i] = 0;
      // Check to see if load host is active on this processor
      for (j=0; j<localIndices.size(); j++) {
        if (p_network->getActiveBus(localIndices[j])) {
          // Check to see if regular load is on this bus
          std::vector<std::string> tags
            = p_network->getBus(localIndices[j])->getLoads();
          for (k = 0; k<tags.size(); k++) {
            if (tags[k] == p_obs_loadIDs[i]) {
              lidx = localIndices[j];
              isLocal = true;
              p_obs_lActive[i] = 1;
              break;
            }
          }
          // Check to see if dynamic load is on this bus
          tags = p_network->getBus(localIndices[j])->getDynamicLoads();
          for (k = 0; k<tags.size(); k++) {
            if (tags[k] == p_obs_loadIDs[i]) {
              lidx = localIndices[j];
              isLocal = true;
              p_obs_lActive[i] = 1;
              break;
            }
          }
        }
        if (isLocal) break;
      }
      if (isLocal) {
        foundLoad[i] = 1;
        p_obs_lLoadIdx.push_back(i);
        p_obs_LoadIdx.push_back(lidx);
        p_obs_lLoadBus.push_back(p_obs_loadBus[i]);
        p_obs_lLoadIDs.push_back(p_obs_loadIDs[i]);
        p_obs_lUse.push_back(0);
        ones.push_back(1.0);
      }
    }
    // Sum foundLoad vector over all processors to find out if any observations
    // do not correspond to existing elements
    p_network->communicator().sum(&foundLoad[0],foundLoad.size());
    // Add observations not found on any processors to observations on process 0
    if (p_comm.rank() == 0 && p_report_dummy_obs) {
      for (i = 0; i<foundLoad.size(); i++) {
        if (!foundLoad[i]) {
          p_obs_lLoadIdx.push_back(i);
          p_obs_LoadIdx.push_back(-1);
          p_obs_lLoadBus.push_back(p_obs_loadBus[i]);
          p_obs_lLoadIDs.push_back(p_obs_loadIDs[i]);
          p_obs_lUse.push_back(1);
          ones.push_back(1.0);
        }
      }
    }
    p_obs_fOnline->addElements(p_obs_lLoadIdx, ones);
    p_obs_fOnline->upload();
    p_comm.sum(&p_obs_lActive[0],p_comm.size());
  }
  if (p_obs_vBus.size() > 0) {
    int nbus = p_obs_vBus.size();
    p_obs_vActive.resize(nbus);
    p_obs_lVBus.clear();
    p_obs_lVIdx.clear();
    p_obs_vUse.clear();
    std::vector<double> zeros;
    std::vector<double> ones;
    int i, j, lidx;
    for (i = 0; i<nbus; i++) {
      std::vector<int> localIndices;
      localIndices = p_network->getLocalBusIndices(p_obs_vBus[i]);
      bool isLocal = false;
      p_obs_vActive[i] = 0;
      // Check to see if bus is active on this processor
      for (j=0; j<localIndices.size(); j++) {
        if (p_network->getActiveBus(localIndices[j])) {
          lidx = localIndices[j];
          p_obs_vActive[i] = 1;
          isLocal = true;
        }
      }
      if (isLocal) {
        foundBus[i] = 1;
        p_obs_lVIdx.push_back(i);
        p_obs_VIdx.push_back(lidx);
        p_obs_lVBus.push_back(p_obs_vBus[i]);
        p_obs_vUse.push_back(0);
        zeros.push_back(0.0);
        ones.push_back(1.0);
      }
    }
    // Sum foundBus vector over all processors to find out if any observations
    // do not correspond to existing elements
    p_network->communicator().sum(&foundBus[0],foundBus.size());
    // Add observations not found on any processors to observations on process 0
    if (p_comm.rank() == 0 && p_report_dummy_obs) {
      for (i = 0; i<foundBus.size(); i++) {
        if (!foundBus[i]) {
          p_obs_lVIdx.push_back(i);
          p_obs_VIdx.push_back(-1);
          p_obs_lVBus.push_back(p_obs_vBus[i]);
          p_obs_vUse.push_back(1);
          zeros.push_back(0.0);
          ones.push_back(1.0);
        }
      }
    }
    p_obs_vMag->addElements(p_obs_lVIdx, ones);
    p_obs_vMag->upload();
    p_obs_vAng->addElements(p_obs_lVIdx, zeros);
    p_obs_vAng->upload();
    p_comm.sum(&p_obs_vActive[0],p_comm.size());
  }
  if (p_obs_vBusfreq.size() > 0) {  // this is for bus frequency ob
    int nbusfreq = p_obs_vBusfreq.size();
    p_obs_vActivefreq.resize(nbusfreq);
    p_obs_lVBusfreq.clear();
    p_obs_lVIdxfreq.clear();
    std::vector<double> ones;
    int i, j, lidx;
    for (i = 0; i<nbusfreq; i++) {
      std::vector<int> localIndicesfreq;
      localIndicesfreq = p_network->getLocalBusIndices(p_obs_vBusfreq[i]);
      bool isLocalfreq = false;
      p_obs_vActivefreq[i] = 0;
      // Check to see if bus is active on this processor
      for (j=0; j<localIndicesfreq.size(); j++) {
        if (p_network->getActiveBus(localIndicesfreq[j])) {
          lidx = localIndicesfreq[j];
          p_obs_vActivefreq[i] = 1;
          isLocalfreq = true;	  		  
        }
      }
      if (isLocalfreq) {
        foundBusfreq[i] = 1;
        p_obs_lVIdxfreq.push_back(i);
        p_obs_VIdxfreq.push_back(lidx);
        p_obs_lVBusfreq.push_back(p_obs_vBus[i]);
        p_obs_vUsefreq.push_back(0);
        ones.push_back(1.0);

        //need to setup Busfreq computation flag for this specific bus  here!!!!!!!!!!!!!!
        p_network->getBus(lidx)->setBusVolFrequencyFlag(true);
		
      }
    }
    // Sum foundBusfreq vector over all processors to find out if any observations
    // do not correspond to existing elements
    p_network->communicator().sum(&foundBusfreq[0],foundBusfreq.size());
    // Add observations not found on any processors to observations on process 0
    if (p_comm.rank() == 0 && p_report_dummy_obs) {
      for (i = 0; i<foundBusfreq.size(); i++) {
        if (!foundBusfreq[i]) {
          p_obs_lVIdxfreq.push_back(i);
          p_obs_VIdxfreq.push_back(-1);
          p_obs_lVBusfreq.push_back(p_obs_vBus[i]);
          p_obs_vUsefreq.push_back(1);
          ones.push_back(1.0);
        }
      }
    }
    p_obs_vBusfreqVal->addElements(p_obs_lVIdxfreq, ones);
    p_obs_vBusfreqVal->upload();
    p_comm.sum(&p_obs_vActivefreq[0],p_comm.size());
  } //bus frequency ob ends here
    
}

/**
 * Get bus and generator IDs for all observations
 * @param genBuses host IDs for all observed generators
 * @param genIDs character identifiers for all observed generators
 * @param loadBuses host IDs for all observed dynamic loads
 * @param loadIDs character identifiers for all observed dynamic loads
 * @param busIDs bus IDs for all observed buses
 */
void gridpack::dynamic_simulation::DSFullApp::getObservationLists(
    std::vector<int> &genBuses, std::vector<std::string> &genIDs,
    std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
    std::vector<int> &busIDs)
{
  genBuses.clear();
  genIDs.clear();
  busIDs.clear();
  int i;
  int nbus = p_obs_genBus.size();
  bool use;
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_gActive[i]);
    }
    if (use) {
      genBuses.push_back(p_obs_genBus[i]);
      genIDs.push_back(p_obs_genIDs[i]);
    }
  }
  nbus = p_obs_loadBus.size();
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_lActive[i]);
    }
    if (use) {
      loadBuses.push_back(p_obs_loadBus[i]);
      loadIDs.push_back(p_obs_loadIDs[i]);
    }
  }
  nbus = p_obs_vBus.size();
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_vActive[i]);
    }
    if (use) {
      busIDs.push_back(p_obs_vBus[i]);
    }
  }
}

/**
 * Get bus and generator IDs for all observations including bus frequency ob
 * @param genBuses host IDs for all observed generators
 * @param genIDs character identifiers for all observed generators
 * @param loadBuses host IDs for all observed dynamic loads
 * @param loadIDs character identifiers for all observed dynamic loads
 * @param busIDs bus IDs for all observed buses
 * @param busfreqIDs bus IDs for all observed buses for bus frequency
 */
void gridpack::dynamic_simulation::DSFullApp::getObservationLists_withBusFreq(
    std::vector<int> &genBuses, std::vector<std::string> &genIDs,
    std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
    std::vector<int> &busIDs, std::vector<int> &busfreqIDs)
{
  genBuses.clear();
  genIDs.clear();
  busIDs.clear();
  busfreqIDs.clear();
  int i;
  int nbus = p_obs_genBus.size();
  bool use;
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_gActive[i]);
    }
    if (use) {
      genBuses.push_back(p_obs_genBus[i]);
      genIDs.push_back(p_obs_genIDs[i]);
    }
  }
  nbus = p_obs_loadBus.size();
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_lActive[i]);
    }
    if (use) {
      loadBuses.push_back(p_obs_loadBus[i]);
      loadIDs.push_back(p_obs_loadIDs[i]);
    }
  }
  nbus = p_obs_vBus.size();
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_vActive[i]);
    }
    if (use) {
      busIDs.push_back(p_obs_vBus[i]);
    }
  }
  
  nbus = p_obs_vBusfreq.size();
  for (i=0; i<nbus; i++) {
    if (p_report_dummy_obs) {
      use = true;
    } else {
      use = static_cast<bool>(p_obs_vActivefreq[i]);
    }
    if (use) {
      busfreqIDs.push_back(p_obs_vBusfreq[i]);
    }
  }
}

/**
 * Get current values of observations
 * @param vMag voltage magnitude for observed buses
 * @param vAng voltage angle for observed buses
 * @param rSpd rotor speed on observed generators
 * @param rAng rotor angle on observed generators
 * @param genP real power on observed generators
 * @param genQ reactive power on observed generators
 * @param fOnline fraction of load shed
 */
void gridpack::dynamic_simulation::DSFullApp::getObservations(
    std::vector<double> &vMag, std::vector<double> &vAng,
    std::vector<double> &rSpd, std::vector<double> &rAng,
	std::vector<double> &genP, std::vector<double> &genQ,
    std::vector<double> &fOnline)
{
  vMag.clear(); 
  vAng.clear(); 
  rSpd.clear(); 
  rAng.clear(); 
  genP.clear(); 
  genQ.clear();
  fOnline.clear(); 
  std::vector<double> tvMag;
  std::vector<double> tvAng;
  std::vector<double> trSpd;
  std::vector<double> trAng;
  std::vector<double> tgenP;
  std::vector<double> tgenQ;
  std::vector<double> tfOnline;
  bool use;
  if (p_obs_genBus.size()) {
    int i, j;
    int nbus =  p_obs_lGenBus.size();
    for (i=0; i<nbus; i++) {
      if (!static_cast<bool>(p_obs_gUse[i])) {
        std::vector<std::string> tags
          = p_network->getBus(p_obs_GenIdx[i])->getGenerators();
        for (j=0; j<tags.size(); j++) {
          if (tags[j] == p_obs_lGenIDs[i]) {
            double speed, angle, gPtmp, gQtmp;
            p_network->getBus(p_obs_GenIdx[i])->getWatchedValues(j,&speed,&angle, &gPtmp, &gQtmp);
            trSpd.push_back(speed);
            trAng.push_back(angle);
            tgenP.push_back(gPtmp);
            tgenQ.push_back(gQtmp);
            break;
          }
        }
      } else {
        trSpd.push_back(0.0);
        trAng.push_back(0.0);
        tgenP.push_back(0.0);
        tgenQ.push_back(0.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lGenIdx.size() != trSpd.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_rSpd->resetElements(p_obs_lGenIdx, trSpd);
    p_obs_rSpd->reload();
    p_obs_rSpd->getAllData(trSpd);
    p_obs_rAng->resetElements(p_obs_lGenIdx, trAng);
    p_obs_rAng->reload();
    p_obs_rAng->getAllData(trAng);
	
	p_obs_genP->resetElements(p_obs_lGenIdx, tgenP);
    p_obs_genP->reload();
    p_obs_genP->getAllData(tgenP);
	p_obs_genQ->resetElements(p_obs_lGenIdx, tgenQ);
    p_obs_genQ->reload();
    p_obs_genQ->getAllData(tgenQ);
	
    nbus = trSpd.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_gActive[i]);
      }
      if (use) {
        rSpd.push_back(trSpd[i]);
        rAng.push_back(trAng[i]);
        genP.push_back(tgenP[i]);
        genQ.push_back(tgenQ[i]);

      }
    }
  }
  if (p_obs_loadBus.size()) {
    int i, j;
    int nbus =  p_obs_lLoadBus.size();
    for (i=0; i<nbus; i++) {
      bool found = false;
      if (!static_cast<bool>(p_obs_lUse[i])) {
        std::vector<std::string> tags
          = p_network->getBus(p_obs_LoadIdx[i])->getDynamicLoads();
        for (j=0; j<tags.size(); j++) {
          if (tags[j] == p_obs_lLoadIDs[i]) {
            double frac;
            frac = p_network->getBus(p_obs_LoadIdx[i])->getOnlineLoadFraction(j);
            tfOnline.push_back(frac);
            found = true;
            break;
          }
        }
      } else {
        tfOnline.push_back(1.0);
        found = true;
      }
      if (!found) {
        tfOnline.push_back(1.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lLoadIdx.size() != tfOnline.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_fOnline->resetElements(p_obs_lLoadIdx, tfOnline);
    p_obs_fOnline->reload();
    p_obs_fOnline->getAllData(tfOnline);
    nbus = tfOnline.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_lActive[i]);
      }
      if (use) {
        fOnline.push_back(tfOnline[i]);
      }
    }
  }
  if (p_obs_vBus.size() > 0) {
    int i, j;
    int nbus =  p_obs_lVBus.size();
    for (i=0; i<nbus; i++) {
      if (!static_cast<bool>(p_obs_vUse[i])) {
        gridpack::ComplexType voltage =
          p_network->getBus(p_obs_VIdx[i])->getComplexVoltage();
        double rV = real(voltage);
        double iV = imag(voltage);
        double V = sqrt(rV*rV+iV*iV);
        double Ang = acos(rV/V);
        if (iV < 0) {
          Ang = -Ang;
        }
        tvMag.push_back(V);
        tvAng.push_back(Ang);
      } else {
        tvMag.push_back(1.0);
        tvAng.push_back(0.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lVIdx.size() != tvMag.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_vMag->resetElements(p_obs_lVIdx, tvMag);
    p_obs_vMag->reload();
    p_obs_vMag->getAllData(tvMag);
    p_obs_vAng->resetElements(p_obs_lVIdx, tvAng);
    p_obs_vAng->reload();
    p_obs_vAng->getAllData(tvAng);
    nbus = tvMag.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_vActive[i]);
      }
      if (use) {
        vMag.push_back(tvMag[i]);
        vAng.push_back(tvAng[i]);
      }
    }
  }
}

/**
 * Get current values of observations including bus frequency ob
 * @param vMag voltage magnitude for observed buses
 * @param vAng voltage angle for observed buses
 * @param rSpd rotor speed on observed generators
 * @param rAng rotor angle on observed generators
 * @param genP real power on observed generators
 * @param genQ reactive power on observed generators
 * @param fOnline fraction of load shed
 * @param busfreq frequency of the buses in ob list
 */
void gridpack::dynamic_simulation::DSFullApp::getObservations_withBusFreq(
    std::vector<double> &vMag, std::vector<double> &vAng,
    std::vector<double> &rSpd, std::vector<double> &rAng,
	std::vector<double> &genP, std::vector<double> &genQ,
    std::vector<double> &fOnline, std::vector<double> &busfreq)
{
  vMag.clear(); 
  vAng.clear(); 
  rSpd.clear(); 
  rAng.clear(); 
  genP.clear(); 
  genQ.clear();
  fOnline.clear(); 
  busfreq.clear();
  std::vector<double> tvMag;
  std::vector<double> tvAng;
  std::vector<double> trSpd;
  std::vector<double> trAng;
  std::vector<double> tgenP;
  std::vector<double> tgenQ;
  std::vector<double> tfOnline;
  std::vector<double> tbusfreq;
  bool use;
  if (p_obs_genBus.size()) {
    int i, j;
    int nbus =  p_obs_lGenBus.size();
    for (i=0; i<nbus; i++) {
      if (!static_cast<bool>(p_obs_gUse[i])) {
        std::vector<std::string> tags
          = p_network->getBus(p_obs_GenIdx[i])->getGenerators();
        for (j=0; j<tags.size(); j++) {
          if (tags[j] == p_obs_lGenIDs[i]) {
            double speed, angle, gPtmp, gQtmp;
            p_network->getBus(p_obs_GenIdx[i])->getWatchedValues(j,&speed,&angle, &gPtmp, &gQtmp);
            trSpd.push_back(speed);
            trAng.push_back(angle);
            tgenP.push_back(gPtmp);
            tgenQ.push_back(gQtmp);
            break;
          }
        }
      } else {
        trSpd.push_back(0.0);
        trAng.push_back(0.0);
        tgenP.push_back(0.0);
        tgenQ.push_back(0.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lGenIdx.size() != trSpd.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_rSpd->resetElements(p_obs_lGenIdx, trSpd);
    p_obs_rSpd->reload();
    p_obs_rSpd->getAllData(trSpd);
    p_obs_rAng->resetElements(p_obs_lGenIdx, trAng);
    p_obs_rAng->reload();
    p_obs_rAng->getAllData(trAng);
	
	p_obs_genP->resetElements(p_obs_lGenIdx, tgenP);
    p_obs_genP->reload();
    p_obs_genP->getAllData(tgenP);
	p_obs_genQ->resetElements(p_obs_lGenIdx, tgenQ);
    p_obs_genQ->reload();
    p_obs_genQ->getAllData(tgenQ);
	
    nbus = trSpd.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_gActive[i]);
      }
      if (use) {
        rSpd.push_back(trSpd[i]);
        rAng.push_back(trAng[i]);
		genP.push_back(tgenP[i]);
		genQ.push_back(tgenQ[i]);
		
      }
    }
  }
  if (p_obs_loadBus.size()) {
    int i, j;
    int nbus =  p_obs_lLoadBus.size();
    for (i=0; i<nbus; i++) {
      bool found = false;
      if (!static_cast<bool>(p_obs_lUse[i])) {
        std::vector<std::string> tags
          = p_network->getBus(p_obs_LoadIdx[i])->getDynamicLoads();
        for (j=0; j<tags.size(); j++) {
          if (tags[j] == p_obs_lLoadIDs[i]) {
            double frac;
            frac = p_network->getBus(p_obs_LoadIdx[i])->getOnlineLoadFraction(j);
            tfOnline.push_back(frac);
            found = true;
            break;
          }
        }
      } else {
        tfOnline.push_back(1.0);
        found = true;
      }
      if (!found) {
        tfOnline.push_back(1.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lLoadIdx.size() != tfOnline.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_fOnline->resetElements(p_obs_lLoadIdx, tfOnline);
    p_obs_fOnline->reload();
    p_obs_fOnline->getAllData(tfOnline);
    nbus = tfOnline.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_lActive[i]);
      }
      if (use) {
        fOnline.push_back(tfOnline[i]);
      }
    }
  }
  
  if (p_obs_vBus.size() > 0) {
    int i, j;
    int nbus =  p_obs_lVBus.size();
    for (i=0; i<nbus; i++) {
      if (!static_cast<bool>(p_obs_vUse[i])) {
        gridpack::ComplexType voltage =
          p_network->getBus(p_obs_VIdx[i])->getComplexVoltage();
        double rV = real(voltage);
        double iV = imag(voltage);
        double V = sqrt(rV*rV+iV*iV);
        double Ang = acos(rV/V);
        if (iV < 0) {
          Ang = -Ang;
        }
        tvMag.push_back(V);
        tvAng.push_back(Ang);
      } else {
        tvMag.push_back(1.0);
        tvAng.push_back(0.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lVIdx.size() != tvMag.size()) {
      printf("Mismatch in vector sizes when resetting global vectors\n");
    }
    p_obs_vMag->resetElements(p_obs_lVIdx, tvMag);
    p_obs_vMag->reload();
    p_obs_vMag->getAllData(tvMag);
    p_obs_vAng->resetElements(p_obs_lVIdx, tvAng);
    p_obs_vAng->reload();
    p_obs_vAng->getAllData(tvAng);
    nbus = tvMag.size();
    for (i=0; i<nbus; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_vActive[i]);
      }
      if (use) {
        vMag.push_back(tvMag[i]);
        vAng.push_back(tvAng[i]);
      }
    }
  }
  
  if (p_obs_vBusfreq.size() > 0) { // bus frequency ob get
    int i, j;
    int nbusfreq =  p_obs_lVBusfreq.size();
	double dbusfreqtmp;
    for (i=0; i<nbusfreq; i++) {
      if (!static_cast<bool>(p_obs_vUsefreq[i])) {
        dbusfreqtmp = p_network->getBus(p_obs_VIdxfreq[i])->getBusVolFrequency();
        tbusfreq.push_back(dbusfreqtmp);
      } else {
        tbusfreq.push_back(0.0);
      }
    }
    // Check to make sure that local vectors still match
    if (p_obs_lVIdxfreq.size() != tbusfreq.size()) {
      printf("Mismatch in vector sizes when resetting global vectors for bus frequency observations\n");
    }
    p_obs_vBusfreqVal->resetElements(p_obs_lVIdxfreq, tbusfreq);
    p_obs_vBusfreqVal->reload();
    p_obs_vBusfreqVal->getAllData(tbusfreq);
	
    nbusfreq = tbusfreq.size();
    for (i=0; i<nbusfreq; i++) {
      if (p_report_dummy_obs) {
        use = true;
      } else {
        use = static_cast<bool>(p_obs_vActivefreq[i]);
      }
      if (use) {
        busfreq.push_back(tbusfreq[i]);
      }
    }
  }
  
  
}

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
void gridpack::dynamic_simulation::DSFullApp::writeRTPRDiagnostics(
    int src_area, int src_zone, int load_area,
    int load_zone, double gen_scale, double load_scale, const char *file)
{
  p_factory->setRTPRParams(src_area,src_zone,load_area,load_zone,
      gen_scale,load_scale);
  p_busIO->open(file);
  double gtotal, ltotal, pmin, pmax, scaled;
  p_factory->getGeneratorMargins(src_area, src_zone,&gtotal,&pmin,&pmax);
  ltotal = p_factory->getTotalLoadRealPower(load_area,load_zone);
  if (gen_scale > 0.0) {
    scaled = gtotal + gen_scale*(pmax-gtotal);
  } else {
    scaled = gtotal + gen_scale*(gtotal-pmin);
  }

  char sbuf[128];
  sprintf(sbuf,"Total Generation:         %16.4f\n",gtotal);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Minimum Generation:     %16.4f\n",pmin);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Maximum Generation:     %16.4f\n",pmax);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Generator Scale Factor: %16.4f\n",gen_scale);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Scaled Generation:      %16.4f\n",scaled);
  p_busIO->header(sbuf);
  p_busIO->header("\nIndividual Scaled Generators\n");
  sprintf(sbuf,"\n     Bus ID   Status Area Zone     Real Power   Scaled Power"
      "           Pmin           Pmax\n\n");
  p_busIO->header(sbuf);
  p_busIO->write("src_gen");
  sprintf(sbuf,"\nTotal Load:               %16.4f\n",ltotal);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Load Scale Factor:      %16.4f\n",load_scale);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Scaled Load:            %16.4f\n",load_scale*ltotal);
  p_busIO->header(sbuf);
  p_busIO->header("\nIndividual Scaled Loads\n");
  sprintf(sbuf,"\n     Bus ID   Status Area Zone     Real Power   Scaled Power"
      " Reactive Power   Scaled Power\n\n");
  p_busIO->header(sbuf);
  p_busIO->write("sink_load");
  p_busIO->close();
}

/**
 * Get a list of buses that had frequency violations
 * @return a list of buses that had frequency failures
 */
std::vector<int> gridpack::dynamic_simulation::DSFullApp::getFrequencyFailures()
{
  std::vector<int> ret;
  gridpack::parallel::Communicator comm = p_network->communicator();
  gridpack::parallel::GlobalVector<int> sumVec(comm);
  int nproc = comm.size();
  int me = comm.rank();
  std::vector<int> sizes(nproc);
  int i;
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = p_violations.size();
  comm.sum(&sizes[0],nproc);

  int offset = 0;
  for (i=1; i<me; i++) offset += sizes[i];
  int total = 0;
  for (i=0; i<nproc; i++) total += sizes[i];
  if (total == 0) return ret;
  std::vector<int> idx;
  int last = offset+sizes[me];
  for (i=offset; i<last; i++) idx.push_back(i);
  sumVec.addElements(idx,p_violations);
  sumVec.upload();
  sumVec.getAllData(ret);
  return ret;
}

/**
 * initialization before the time step integration starts 
 */
void gridpack::dynamic_simulation::DSFullApp::solvePreInitialize(
    gridpack::dynamic_simulation::Event fault)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

  t_misc = timer->createCategory("DS Solve: Miscellaneous");
  t_presolve = timer->createCategory("DS App: solvePreInitialize");
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif

  timer->start(t_misc);
  timer->start(t_presolve);

  // Get cursor for setting solver options
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  timer->stop(t_misc);

  t_mode = timer->createCategory("DS Solve: Set Mode");
  timer->start(t_mode);
  p_factory->setMode(YBUS);
  timer->stop(t_mode);
  t_ybus = timer->createCategory("DS Solve: Make YBus");
  timer->start(t_ybus);
  
  // set the line trip action related flag to be false and clear the vector
  bapplyLineTripAction = false;
  p_vbranches_need_to_trip.clear();
  
  // set the load P and Q change related flag to be false and clear the vector
  bapplyLoadChangeP = false;
  bapplyLoadChangeQ = false;
  p_vbus_need_to_changeP.clear();
  p_vbus_need_to_changeQ.clear();
  
  ybusMap_sptr.reset(new gridpack::mapper::FullMatrixMap<DSFullNetwork> (p_network));
  orgYbus = ybusMap_sptr->mapToMatrix();
  
  //printf("\n=== org ybus: ============\n");
  //orgYbus->print();
  //orgYbus->save("ybus_GridPACK_org.m");
  //exit(0);

  //p_factory->addLoadAdmittance();

  // Form constant impedance load admittance yl for all buses and add it to
  // system Y matrix: ybus = ybus + yl
  p_factory->setMode(YL);
  ybusyl = ybusMap_sptr->mapToMatrix();
  timer->stop(t_ybus);
  //branchIO.header("\n=== ybus after added yl: ============\n");
  //printf("\n=== ybus after added yl: ============\n");
  //ybusyl->print();
  //ybusyl->save("ybus_GridPACK_yl.m");
  //exit(0);

  p_factory->setMode(PG);
  ybuspg = ybusMap_sptr->mapToMatrix();
  //printf("\n=== ybus after added pg: ============\n");
  //ybuspg->print();
  //exit(0);

  //printf("# of buses in the network: %d\n", p_network->totalBuses());

  // Add j*Xd' to system Y matrix:
  // Extract appropriate xdprime and xdpprime from machine data
  timer->start(t_mode);
  p_factory->setMode(jxd);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybus_jxd = ybusMap_sptr->mapToMatrix();
  //branchIO.header("\n=== ybusyl after added j*Xd': =============\n");
  //printf("\n=== ybusyl after added j*Xd': =============\n");
  //ybus_jxd->print();
  //ybus_jxd->save("ybus_GridPACK_jxd.m");
  
  
  // Add dynamic load impedance to system Y matrix:
  timer->start(t_mode);
  p_factory->setMode(YDYNLOAD);
  timer->stop(t_mode);
  ybus = ybusMap_sptr->mapToMatrix();
  //branchIO.header("\n=== ybus_jxd after added dynamic load impedance': =============\n");
  //printf("\n=== ybus_dynload after added dynamic load impedance': =============\n");
  //ybus->print();
  //ybus->save("ybus_GridPACK_dynload.m");
  
  //exit(0);

  // Compute ybus_fy for fault on stage
  ybus_fy.reset(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_misc);
  p_factory->setEvent(fault);
  timer->stop(t_misc);
  timer->start(t_mode);
  p_factory->setMode(onFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap_sptr->overwriteMatrix(ybus_fy);
  //branchIO.header("\n=== ybus_fy: ============\n");
  //printf("\n=== ybus_fy: ============\n");
  //ybus_fy->print();
  //ybus_fy->save("ybus_fy_GridPACK_jxd.m");

  // Compute ybus_posfy for fault clear stage
  ybus_posfy.reset(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_mode);
  p_factory->setMode(posFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap_sptr->incrementMatrix(ybus_posfy);
  //branchIO.header("\n=== ybus_posfy: ============\n");
  //printf("\n=== ybus_posfy: ============\n");
  //ybus_posfy->print();
  //ybus_posfy->save("ybus_posfy_GridPACK_jxd.m");
  timer->stop(t_ybus);

  // Simulation related variables
  t_init = timer->createCategory("DS Solve: Initialization");
  timer->start(t_init);
  
  int t_step[20];
  double t_width[20];

  //const double sysFreq = 60.0;
  //double pi = 4.0*atan(1.0);
  //const double basrad = 2.0 * pi * sysFreq;
  //gridpack::ComplexType jay(0.0, 1.0);

  // switch info is set up here
  int nswtch = 4;
  static double sw1[4];
  static double sw7[4];
  sw1[0] = 0.0;
  sw1[1] = fault.start;
  sw1[2] = fault.end;
  sw1[3] = p_sim_time;
  sw7[0] = p_time_step;
  sw7[1] = p_time_step;
  sw7[2] = p_time_step;
  sw7[3] = p_time_step;
  simu_total_steps = 0;
  for (int i = 0; i < nswtch-1; i++) {
    t_step[i] = (int) ((sw1[i+1] -sw1[i]) / sw7[i]);
    t_width[i] = (sw1[i+1] - sw1[i]) / t_step[i];
    simu_total_steps += t_step[i];
  }
  simu_total_steps++;
  
  // Initialize vectors for integration 
  p_factory->initDSVect(p_time_step);
  
  p_factory->setGeneratorObPowerBaseFlag(p_generator_observationpower_systembase);
  //exit(0);

  ngenMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork> (p_network));
  
  // Map to create vector volt
  volt = ngenMap_sptr->mapToVector();
  //p_busIO->header("\n=== volt: ===\n");
  //volt->print();

  solver_sptr.reset(new gridpack::math::LinearSolver (*ybus));
  solver_sptr->configure(cursor);
  
  //gridpack::math::LinearSolver solver_fy(*ybus_fy);
  solver_fy_sptr.reset(new gridpack::math::LinearSolver (*ybus_fy));
  solver_fy_sptr->configure(cursor);
  
  //gridpack::math::LinearSolver solver_posfy(*ybus_posfy);
  //gridpack::math::LinearSolver solver_posfy(*ybus); 
  solver_posfy_sptr.reset(new gridpack::math::LinearSolver (*ybus));
  solver_posfy_sptr->configure(cursor);

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagP = 0;
  flagC = 0;
  S_Steps = 1;
  last_S_Steps = -1;

  p_insecureAt = -1;

  p_factory->setMode(make_INorton_full);
  //gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
  nbusMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork>(p_network));
  INorton_full = nbusMap_sptr->mapToVector();
  INorton_full_chk = nbusMap_sptr->mapToVector();
  max_INorton_full = 0.0;
  volt_full.reset(INorton_full->clone());

  timer->stop(t_init);
  if (!p_suppress_watch_files) {
#ifdef USE_TIMESTAMP
    if (p_generatorWatch) p_generatorIO->header("t, t_stamp");//bus_id,ckt,x1d_1,x2w_1,x3Eqp_1,x4Psidp_1,x5Psiqpp_1");
    //#  if (p_generatorWatch) p_generatorIO->header("t, t_stamp,bus_id,ckt,x1d_1,x2w_1,x3Eqp_1,x4Psidp_1,x5Psiqpp_1");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    if (p_generatorWatch) p_generatorIO->header("\n");

    if (p_loadWatch) p_loadIO->header("t, t_stamp");
    if (p_loadWatch) p_loadIO->write("load_watch_header");
    if (p_loadWatch) p_loadIO->header("\n");
#else
    if (p_generatorWatch) p_generatorIO->header("t");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    if (p_generatorWatch) p_generatorIO->header("\n");

    if (p_loadWatch) p_loadIO->header("t");
    if (p_loadWatch) p_loadIO->write("load_watch_header");
    if (p_loadWatch) p_loadIO->header("\n");
#endif
  }
#ifdef USE_GOSS
  if (p_generatorWatch) p_generatorIO->dumpChannel();
  if (p_loadWatch) p_loadIO->dumpChannel();
#endif
  p_frequencyOK = true;
  // Save initial time step
  //saveTimeStep();
	
  Simu_Current_Step = 0;
  p_bDynSimuDone = false;
  
  timer->stop(t_presolve);
  //printf (" In function solvePreInitialize end, simu_total_steps: %d \n", simu_total_steps);
  
}

/**
 * Execute only one simulation time step 
 */
void gridpack::dynamic_simulation::DSFullApp::executeOneSimuStep( ){
	
    gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

    t_execute_steps = timer->createCategory("DS Solve: execute steps");
	timer->start(t_execute_steps);
	
  //for (Simu_Current_Step = 0; Simu_Current_Step < simu_total_steps - 1; Simu_Current_Step++) {
  //for (Simu_Current_Step = 0; Simu_Current_Step < 200; Simu_Current_Step++) {
    //char step_str[128];
    //sprintf(step_str,"\nIter %d\n", Simu_Current_Step);
    //p_busIO->header(step_str);
    timer->start(t_misc);
    //printf("Step %d\ttime %5.3f sec: \n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
    //printf("\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
    ///char step_str[128];
    ///sprintf(step_str, "\n===================Step %d\ttime %5.3f sec:================\n", Simu_Current_Step+1, (Simu_Current_Step+1) * p_time_step);
     ///p_busIO->header(step_str);
    S_Steps = Simu_Current_Step;

    if (Simu_Current_Step < steps1) {
      flagP = 0;
      flagC = 0;
    } else if (Simu_Current_Step == steps1) {
      flagP = 0;
      //flagC = 1;
      flagC = 0;
    } else if ((Simu_Current_Step > steps1) && (Simu_Current_Step < steps2)) {
      flagP = 1;
      flagC = 1;
    } else if (Simu_Current_Step == steps2) {
      flagP = 1;
      //flagC = 2;
      flagC = 1;
    } else if (Simu_Current_Step > steps2) {
      flagP = 2;
      flagC = 2;
    }
    timer->stop(t_misc);
	
	// renke add, if a line trip action is detected, modify the post-fault Ymatrix. 
	// Here we assume line trip action will only happen AFTER FAULT!!!!!!!
	if (bapplyLineTripAction){
		
		char sybus[100];
		
		//sprintf(sybus, "ybus_%d_before_linetrip.m",Simu_Current_Step );
		//ybus->save(sybus);
		
		p_factory->setMode(branch_trip_action);
        ybusMap_sptr->incrementMatrix(ybus);  // in the current code, solver_posfy_sptr is linked with ybus, check Bill
		
		//printf ("-----------renke debug, line tripping action detected----------");		
        //ybus->print();       
        //sprintf(sybus, "ybus_%d_linetrip_test.m",Simu_Current_Step );
        //ybus->save(sybus);
		
		// after Y-matrix is modified, we need to clear this line trip action to 
		// avoid next step still apply the same line trip action
		clearLineTripAction();// in this one, need to clear the flag, vector of each branch and set the status of the branches to be 0);  

	}
   bapplyLineTripAction = false;
   
   
   	// renke add, if a load change P action is detected, modify the post-fault Ymatrix. 
	// Here we assume line trip action will only happen AFTER FAULT!!!!!!!
	if (bapplyLoadChangeP){
		
		char sybus[100];
		
		//sprintf(sybus, "ybus_%d_before_linetrip.m",Simu_Current_Step );
		//ybus->save(sybus);
		
		p_factory->setMode(bus_Yload_change_P);
        ybusMap_sptr->incrementMatrix(ybus);  // in the current code, solver_posfy_sptr is linked with ybus, check Bill
		
		//printf ("-----------renke debug, line tripping action detected----------");		
        //ybus->print();       
        //sprintf(sybus, "ybus_%d_linetrip_test.m",Simu_Current_Step );
        //ybus->save(sybus);
		
		// after Y-matrix is modified, we need to clear this line trip action to 
		// avoid next step still apply the same line trip action
		clearConstYLoad_Change_P();// in this one, need to clear the flag, vector of each branch and set the status of the branches to be 0);  

	}
   bapplyLoadChangeP = false;
   
      	// renke add, if a load change Q action is detected, modify the post-fault Ymatrix. 
	// Here we assume line trip action will only happen AFTER FAULT!!!!!!!
	if (bapplyLoadChangeQ){
		
		char sybus[100];
		
		//sprintf(sybus, "ybus_%d_before_linetrip.m",Simu_Current_Step );
		//ybus->save(sybus);
		
		p_factory->setMode(bus_Yload_change_Q);
        ybusMap_sptr->incrementMatrix(ybus);  // in the current code, solver_posfy_sptr is linked with ybus, check Bill
		
		//printf ("-----------renke debug, line tripping action detected----------");		
        //ybus->print();       
        //sprintf(sybus, "ybus_%d_linetrip_test.m",Simu_Current_Step );
        //ybus->save(sybus);
		
		// after Y-matrix is modified, we need to clear this line trip action to 
		// avoid next step still apply the same line trip action
		clearConstYLoad_Change_Q();// in this one, need to clear the flag, vector of each branch and set the status of the branches to be 0);  

	}
   bapplyLoadChangeQ = false;
   
   
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor_currentInjection(false);
    } else {
      p_factory->predictor_currentInjection(true);
    }

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    t_mIf = timer->createCategory("DS Solve: Modified Euler Predictor: Make INorton");
    timer->start(t_mIf);
	p_factory->setMode(make_INorton_full);
    nbusMap_sptr->mapToVector(INorton_full);
    ///gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
    ///boost::shared_ptr<gridpack::math::Vector> INorton_full = nbusMap_sptr->mapToVector();
    //p_busIO->header("\n=== [Predictor] INorton_full: ===\n");
    //printf("renke test \n=== [Predictor] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_mIf);
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif
 
    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2) 
    // to calculate terminal volt: ----------
    t_psolve = timer->createCategory("DS Solve: Modified Euler Predictor: Linear Solver");
    timer->start(t_psolve);
    //boost::shared_ptr<gridpack::math::Vector> volt_full(INorton_full->clone());
    volt_full->zero();
	
//!!!!!!!!!!iteratively solve or not, rk!!!!
 //double ITER_TOL = 1.0e-7;
 //int MAX_ITR_NO = 8;
 
 //printf ("-----rk debug in gridpack::dynamic_simulation::DSFullApp::executeOneSimuStep( ): ITER_TOL: %15.12f, MAX_ITR_NO: %d \n\n", ITER_TOL, MAX_ITR_NO);
 bool flag_chk;
 int iter_num_record;
 //bool p_iterative_network_debug = false;
 if (p_biterative_solve_network){
    flag_chk = true;
	iter_num_record = 0;
	if (p_iterative_network_debug) printf ("--------------------------in iterative predictor current injection, Simu_Current_Step: %d------------------- \n", Simu_Current_Step);
    while (flag_chk == true ) {			
      if (flagP == 0) {
	solver_sptr->solve(*INorton_full, *volt_full);
      } else if (flagP == 1) {
	solver_fy_sptr->solve(*INorton_full, *volt_full);
      } else if (flagP == 2) {
	solver_posfy_sptr->solve(*INorton_full, *volt_full);
      }
			
      
      //printf("1: itr test:----previous predictor_INorton_full:\n");
      //INorton_full->print();
      
      INorton_full_chk->equate(*INorton_full);
      //printf("2: itr test:----predictor_INorton_full_chk:\n");
      //INorton_full_chk->print();

      nbusMap_sptr->mapToBus(volt_full);
      p_factory->setVolt(false);
      
      if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
	p_factory->predictor_currentInjection(false);
      } else {
	p_factory->predictor_currentInjection(true);
      }
      
# if 0	
      printf("3: itr test:----previous predictor_INorton_full:\n");
      INorton_full->print();
      
      INorton_full_chk->equate(*INorton_full);
      printf("4: itr test:----predictor_INorton_full_chk:\n");
      INorton_full_chk->print();
# endif			
      p_factory->setMode(make_INorton_full);
      nbusMap_sptr->mapToVector(INorton_full);
      
      //printf("5: itr test:----predictor_INorton_full:\n");
      //INorton_full->print();
      
      //multiply(*ybus_fy, *volt_full, *INorton_full_chk);
      INorton_full_chk->add(*INorton_full, -1.0);
      max_INorton_full=abs(INorton_full_chk->normInfinity());
      
      if (max_INorton_full < ITER_TOL || iter_num_record >= MAX_ITR_NO) {
	flag_chk = false;
      } else {
	iter_num_record += 1;
	if (p_iterative_network_debug) printf("              predictor iter number: %d, max_INorton_full = %13.10f \n", iter_num_record, max_INorton_full);
	//printf("-----INorton_full : \n");
	//INorton_full->print();
	//printf("-----INorton_full_chk - INorton_full : \n");
	//INorton_full_chk->print();
      }
    }// end of while
 }else {// p_biterative_solve_network = false
    if (flagP == 0) {
      solver_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy_sptr->solve(*INorton_full, *volt_full);
    }
  } // end of if (p_biterative_solve_network)
    timer->stop(t_psolve);

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    //p_busIO->header("\n=== [Predictor] volt_full: ===\n");
    //volt_full->print();
    //if (Simu_Current_Step==4){
    //	 exit(0);
   //	}

    t_vmap= timer->createCategory("DS Solve: Map Volt to Bus");
    timer->start(t_vmap);
	
	//printf("after first volt sovle, before first volt map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap_sptr->mapToBus(volt_full);
	
	//printf("after first volt sovle, after first volt map: \n");
	
	if ( Simu_Current_Step==0 ) {
		//printf("enter the initial update oldbusvoltage, Timestep: %d \n", Simu_Current_Step);
		p_factory->updateoldbusvoltage(); //renke add, first timestep, copy volt_full to volt_full_old
	}
    timer->stop(t_vmap);

    t_volt= timer->createCategory("DS Solve: Set Volt");
    timer->start(t_volt);
    p_factory->setVolt(false);
	p_factory->updateBusFreq(h_sol1);
	
	
	std::vector <double> vwideareafreqs;
	vwideareafreqs = p_factory->grabWideAreaFreq();
	//printf("-----!!renke debug dsf_app_module.cpp: grabWideAreaFreq: bus 30: %12.6f, bus 30: %12.6f, delta_freq bus34-bus30: %12.6f \n", 
	//		vwideareafreqs[0], vwideareafreqs[1], vwideareafreqs[2]);
	int tmp = vwideareafreqs.size();
	double widearea_deltafreq = vwideareafreqs[tmp-1];

	//p_factory->setWideAreaFreqforPSS(widearea_deltafreq);
	 		
    timer->stop(t_volt);
	
	//printf("before update relay, after first volt solv: \n");
	//p_factory->printallbusvoltage();
    //renke add, compute bus freq if necessary
    //printf("Timestep, %d \n", Simu_Current_Step);
    bool flagBus = p_factory->updateBusRelay(false, h_sol1);
    bool flagBranch = p_factory->updateBranchRelay(false, h_sol1);
	
	// update dynamic load internal relay functions here
	p_factory->dynamicload_post_process(h_sol1, false);
    
	// if bus relay trips, modify the corresponding Ymatrix, renke modified
    if (flagBus) {
        printf("DSFull_APP::Solve: updatebusrelay return trigger siganl: TURE!!! \n");
		
        //please update the bus contribution to the Y bus matrix here. //Shuangshuang tbd
	if (flagP == 0) { 
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",Simu_Current_Step );
			 
	     ybus_fy->save(sybus);
	 
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
	     ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

	     ybus->save(sybus);
        
             printf("DSFull_APP::Solve: bus relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

			 
        } else if (flagP == 2) {
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus);
             printf("DSFull_APP::Solve: bus relay trip after fault, ybus changed:\n");
	     ybus->print();
	     char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

	     ybus->save(sybus);

             printf("DSFull_APP::Solve: bus relay trip after fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap_sptr->overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

        }
    }
	
	// if branch relay trips, modify the corresponding Ymatrix, renke modified
	if (flagBranch) {
        
        printf("DSFull_APP::Solve: updatebranchrelay return trigger siganl: TURE!!! \n");

        //please update the bus contribution to the Y bus matrix here. //Shuangshuang tbd
	if (flagP == 0) { 
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",Simu_Current_Step );
			 
	     ybus_fy->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
             ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

             ybus->save(sybus);

			 
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_posfy);
	     ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

	     ybus_posfy->save(sybus);
			 
        } else if (flagP == 2) {
             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus);
             ybus->print();
             char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",Simu_Current_Step );

             ybus->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap_sptr->incrementMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",Simu_Current_Step );

             ybus_posfy->save(sybus);

        }
    }
	
    //renke add, update old busvoltage first
    p_factory->updateoldbusvoltage(); //renke add
	
	//printf("after updateoldbus voltage: \n");
	//p_factory->printallbusvoltage();
	
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif

    t_predictor = timer->createCategory("DS Solve: Modified Euler Predictor");
    //printf("Test: predictor begins: \n");
    timer->start(t_predictor);
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor(h_sol1, false);
    } else { 
      p_factory->predictor(h_sol1, true);
    }
    timer->stop(t_predictor);

    
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->corrector_currentInjection(false);
    } else {
      p_factory->corrector_currentInjection(true);
    }
    

    //INorton_full = nbusMap_sptr->mapToVector();
    t_cmIf = timer->createCategory("DS Solve: Modified Euler Corrector: Make INorton");
    timer->start(t_cmIf);
    p_factory->setMode(make_INorton_full);
    nbusMap_sptr->mapToVector(INorton_full);
    //p_busIO->header("\n=== [Corrector] INorton_full: ===\n");
    //printf("\nrelaytest=== [Corrector] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_cmIf);

    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2)
    // to calculate terminal volt: ----------
    t_csolve = timer->createCategory("DS Solve: Modified Euler Corrector: Linear Solver");
    timer->start(t_csolve);
    
    volt_full->zero();

//!!!!!!!!!!iteratively solve or not, rk!!!!
  if (p_biterative_solve_network){
    flag_chk = true;
	iter_num_record = 0;
    while (flag_chk == true ) {
		
			volt_full->zero();
			
			if (flagP == 0) {
				solver_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 1) {
				solver_fy_sptr->solve(*INorton_full, *volt_full);
			} else if (flagP == 2) {
				solver_posfy_sptr->solve(*INorton_full, *volt_full);
			}
			nbusMap_sptr->mapToBus(volt_full);
			p_factory->setVolt(false);
			
			if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
				p_factory->corrector_currentInjection(false);
			} else {
				p_factory->corrector_currentInjection(true);
			}
			
			INorton_full_chk->equate(*INorton_full);
			//printf("itr test:----corrector_INorton_full_chk:\n");
			//INorton_full_chk->print();
			
			p_factory->setMode(make_INorton_full);
			nbusMap_sptr->mapToVector(INorton_full);
			
			//printf("itr test:----corrector_INorton_full:\n");
			//INorton_full->print();
			
			//multiply(*ybus_fy, *volt_full, *INorton_full_chk);
			INorton_full_chk->add(*INorton_full, -1.0);
			max_INorton_full=abs(INorton_full_chk->normInfinity());
			
			if (max_INorton_full < ITER_TOL || iter_num_record>= MAX_ITR_NO ) {
				flag_chk = false;
			} else {
				iter_num_record += 1;
				if (p_iterative_network_debug)  printf("                corrector: iter number: %d, max_INorton_full = %13.10f \n", iter_num_record, max_INorton_full);
				//printf("-----INorton_full : \n");
				//INorton_full->print();
				//printf("-----INorton_full_chk - INorton_full : \n");
				//INorton_full_chk->print();
			}
    }// end of while
  }else{ // p_biterative_solve_network = false
    if (flagP == 0) {
      solver_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy_sptr->solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy_sptr->solve(*INorton_full, *volt_full);
    }
  } //p_biterative_solve_network end here

    timer->stop(t_csolve);

    //p_busIO->header("\n=== [Corrector] volt_full: ===\n");
    //printf("relaytest \n=== [Corrector] volt_full: ===\n");
    //volt_full->print();
    timer->start(t_vmap);
	
	//printf("after second solve, before second map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap_sptr->mapToBus(volt_full);
	
	//printf("after second solve, after second map: \n");
	//p_factory->printallbusvoltage();
	
    timer->stop(t_vmap);

    timer->start(t_volt);
    p_factory->setVolt(false);
	p_factory->updateBusFreq(h_sol1);
    timer->stop(t_volt);

    t_corrector = timer->createCategory("DS Solve: Modified Euler Corrector");
    timer->start(t_corrector);
    //printf("Test: corrector begins: \n");
    if (last_S_Steps != S_Steps) {
      p_factory->corrector(h_sol2, false);
    } else {
      p_factory->corrector(h_sol2, true);
    }
    timer->stop(t_corrector);

    //if (Simu_Current_Step == simu_total_steps - 1) 
      //p_busIO->write();


    //printf("----------!renke debug, after solve INorton_full and map back voltage ----------\n");
    t_secure = timer->createCategory("DS Solve: Check Security");
    timer->start(t_secure);
    if (p_generatorWatch && Simu_Current_Step%p_generatorWatchFrequency == 0) {
      char tbuf[32];
      if (!p_suppress_watch_files) {
#ifdef USE_TIMESTAMP
        sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(Simu_Current_Step)*p_time_step,
            timer->currentTime());
        if (p_generatorWatch) p_generatorIO->header(tbuf);
        if (p_generatorWatch) p_generatorIO->write("watch");
        if (p_generatorWatch) p_generatorIO->header("\n");

        //      if (p_generatorWatch) p_generatorIO->write("watch");

        //      sprintf(tbuf,"%8.4f, %20.4f",mac_ang_s0, mac_spd_s0);
        //      if (p_generatorWatch) p_generatorIO->header(tbuf);
        //      if (p_generatorWatch) p_generatorIO->write("watch");
        //      if (p_generatorWatch) p_generatorIO->header("\n");
#else
        sprintf(tbuf,"%8.4f",static_cast<double>(Simu_Current_Step)*p_time_step);
        if (p_generatorWatch) p_generatorIO->header(tbuf);
        if (p_generatorWatch) p_generatorIO->write("watch");
        if (p_generatorWatch) p_generatorIO->header("\n");
#endif
      }
#ifdef USE_GOSS
      if (p_generatorWatch) p_generatorIO->dumpChannel();
#endif
    }
    //printf("------------------!!!renke debug after the generator watch ------------\n");

    if (p_loadWatch && Simu_Current_Step%p_loadWatchFrequency == 0) {
      char tbuf[32];
      if (!p_suppress_watch_files) {
#ifdef USE_TIMESTAMP
        sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(Simu_Current_Step)*p_time_step,
            timer->currentTime());
        if (p_loadWatch) p_loadIO->header(tbuf);
        if (p_loadWatch) p_loadIO->write("load_watch");
        if (p_loadWatch) p_loadIO->header("\n");
#else
        sprintf(tbuf,"%8.4f",static_cast<double>(Simu_Current_Step)*p_time_step);
        if (p_loadWatch) p_loadIO->header(tbuf);
        if (p_loadWatch) p_loadIO->write("load_watch");
        if (p_loadWatch) p_loadIO->header("\n");
#endif
      }
#ifdef USE_GOSS
      if (p_loadWatch) p_loadIO->dumpChannel();
#endif
    }
    saveTimeStep();
    //    if ((!p_factory->securityCheck()) && p_insecureAt == -1)  
    //       p_insecureAt = Simu_Current_Step;

    
/*    // Print to screen
    if (last_S_Steps != S_Steps) {
      //sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps);
      //p_busIO->header(ioBuf);
      printf("\n==============S_Steps = %d==============\n", S_Steps);
      mac_ang_s0->print();
      mac_spd_s0->print();
      //pmech->print();
      //pelect->print();
      //sprintf(ioBuf, "========================End of S_Steps = %d=========================\n\n", S_Steps);
      //p_busIO->header(ioBuf);
    }
    if (Simu_Current_Step == simu_total_steps) {
      printf("\n==============S_Steps = %d==============\n", S_Steps);
      mac_ang_s1->print();
      mac_spd_s1->print();
      p_factory->setMode(init_mac_ang);
      ngenMap_sptr->mapToBus(mac_ang_s1);
      p_factory->setMode(init_mac_spd);
      ngenMap_sptr->mapToBus(mac_spd_s1);
      p_factory->setMode(init_pmech);
      ngenMap_sptr->mapToBus(pmech);
      p_factory->setMode(init_pelect);
      ngenMap_sptr->mapToBus(pelect);
      sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps+1);
      p_busIO->header(ioBuf);
      sprintf(ioBuf, "\n         Bus ID     Generator ID"
          "    mac_ang         mac_spd         mech            elect\n\n");
      p_busIO->header(ioBuf);
      mac_ang_s1->print();
      mac_spd_s1->print();
      pmech->print();
      pelect->print();
      p_busIO->write();
      sprintf(ioBuf, "\n========================End of S_Steps = %d=========================\n\n", S_Steps+1);
      p_busIO->header(ioBuf);
    } // End of Print to screen

*/    //exit(0);
    last_S_Steps = S_Steps;
    timer->stop(t_secure);
    if (p_monitorGenerators) {
      double presentTime = static_cast<double>(Simu_Current_Step)*p_time_step;
      p_frequencyOK = p_frequencyOK && checkFrequency(0.5,presentTime);
      if (!p_frequencyOK) Simu_Current_Step = simu_total_steps;
    }
  //} // main for loop ends here
  
  Simu_Current_Step++;
  
  if (Simu_Current_Step >= simu_total_steps - 1){
	  p_bDynSimuDone = true;
  }
	  
  
#if 0
  printf("\n=== ybus after simu: ============\n");
  ybus->print();
  ybus->save("ybus_aftersimu.m");
  
  printf("\n=== ybus_fy after simu:============\n");
  ybus_fy->print();
  ybus_fy->save("ybus_fy_aftersimu.m");
  
  printf("\n=== ybus_posfy after simu: ============\n");
  ybus_posfy->print();
  ybus_posfy->save("ybus_posfy_aftersimu.m");
  
#endif

  //char msg[128];
  //if (p_insecureAt == -1) sprintf(msg, "\nThe system is secure!\n");
  //else sprintf(msg, "\nThe system is insecure from step %d!\n", p_insecureAt);
  //char secureBuf[128];
  //if (p_insecureAt == -1) sprintf(secureBuf,"\nThe system is secure!\n");
  //else sprintf(secureBuf,"\nThe system is insecure from step %d!\n", p_insecureAt);
  //p_busIO->header(secureBuf);

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
  timer->stop(t_execute_steps);
  //timer->dump(); 
  
}

/**
 * Check whether the dynamic simulation is done
 */
bool gridpack::dynamic_simulation::DSFullApp::isDynSimuDone( ){
	
	return p_bDynSimuDone;

}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ
 */
 
void gridpack::dynamic_simulation::DSFullApp::scatterInjectionLoad(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
	
	std::vector<int> vec_busintidx;
	int ival, nvals, ibus, nbus, bus_number;
	gridpack::dynamic_simulation::DSFullBus *bus;	

	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		vec_busintidx = p_network->getLocalBusIndices(bus_number);
		nbus = vec_busintidx.size();
		for(ibus=0; ibus<nbus; ibus++){
			bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
			//printf("----renke debug scatterInjectionLoad, in dsf full app, \n");
			bus->scatterInjectionLoad(vloadP[ival], vloadQ[ival]);
		
		}
	}
		
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ - this implemnetation keeps the Y load component of the bus still at the bus, 
 * while only compenstates the difference
 */
 
void gridpack::dynamic_simulation::DSFullApp::scatterInjectionLoadNew_compensateY(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
		
	std::vector<int> vec_busintidx;
	int ival, nvals, ibus, nbus, bus_number;
	gridpack::dynamic_simulation::DSFullBus *bus;
	double orgp, orgq;
	
	// treat the new load p and q as current source
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		vec_busintidx = p_network->getLocalBusIndices(bus_number);
		nbus = vec_busintidx.size();
		for(ibus=0; ibus<nbus; ibus++){
			bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
			//printf("----renke debug scatterInjectionLoad, in dsf full app, \n");
			bus->scatterInjectionLoad_compensateY(vloadP[ival], vloadQ[ival]);
		
		}
	}
		
}


/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
 * and model the entire load change as injection current
 */
 
void gridpack::dynamic_simulation::DSFullApp::scatterInjectionLoadNew(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
		
	std::vector<int> vec_busintidx;
	int ival, nvals, ibus, nbus, bus_number;
	gridpack::dynamic_simulation::DSFullBus *bus;
	double orgp, orgq;
	
	//first modify the original values of the Contant Y load P and Q to zero, 
	// note: only the first time receive the command of scatter InjectionLoadNew needs to do the clear of the original load values!!!!!!
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		setConstYLoadtoZero_P(bus_number);
		setConstYLoadtoZero_Q(bus_number);
	}
	
	// treat the new load p and q as current source
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		vec_busintidx = p_network->getLocalBusIndices(bus_number);
		nbus = vec_busintidx.size();
		for(ibus=0; ibus<nbus; ibus++){
			bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
			//printf("----renke debug scatterInjectionLoad, in dsf full app, \n");
			bus->scatterInjectionLoad(vloadP[ival], vloadQ[ival]);
		
		}
	}
		
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
 * and model the entire load change as injection current, also add a large parallel connecting impedance at the certain bus
*/
void gridpack::dynamic_simulation::DSFullApp::scatterInjectionLoadNew_Norton(const std::vector<int>& vbusNum, 
							const std::vector<double>& vloadP, const std::vector<double>& vloadQ, 
							const std::vector<double>& vimpedanceR, const std::vector<double>& vimpedanceI){
	
	std::vector<int> vec_busintidx;
	int ival, nvals, ibus, nbus, bus_number;
	gridpack::dynamic_simulation::DSFullBus *bus;
	double orgp, orgq, impr, impi, impedancer, impedancei;
	
	//first modify the original values of the Contant Y load P and Q to zero, 
	// note: only the first time receive the command of scatter InjectionLoadNew needs to do the clear of the original load values!!!!!!
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		impedancer = vimpedanceR[ival];
		impedancei = vimpedanceI[ival];
		//setConstYLoadtoZero_P(bus_number);
		//setConstYLoadtoZero_Q(bus_number);
		setConstYLoadImpedance(bus_number, impedancer, impedancei);
	}
	
	// treat the new load p and q as current source
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		vec_busintidx = p_network->getLocalBusIndices(bus_number);
		nbus = vec_busintidx.size();
		for(ibus=0; ibus<nbus; ibus++){
			bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
			//printf("----renke debug scatterInjectionLoad, in dsf full app, \n");
			bus->scatterInjectionLoad(vloadP[ival], vloadQ[ival]);
		
		}
	}	
								
}

/**
 * execute load scattering with constant current load , the values of the STATIC load current at certain buses vbusNum will be changed to the values of 
 * the vector  vCurR and vCurI - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
 * and model the entire load change as injection current
 */
 
void gridpack::dynamic_simulation::DSFullApp::scatterInjectionLoadNewConstCur(const std::vector<int>& vbusNum, const std::vector<double>& vCurR, const std::vector<double>& vCurI){
		
	std::vector<int> vec_busintidx;
	int ival, nvals, ibus, nbus, bus_number;
	gridpack::dynamic_simulation::DSFullBus *bus;
	
	//first modify the original values of the Contant Y load P and Q to zero, 
	// note: only the first time receive the command of scatter InjectionLoadNew needs to do the clear of the original load values!!!!!!
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		setConstYLoadtoZero_P(bus_number);
		setConstYLoadtoZero_Q(bus_number);
	}
	
	// treat the new load p and q as constant current source
	nvals = vbusNum.size();	
	for (ival=0; ival<=nvals; ival++){
		bus_number = vbusNum[ival];
		vec_busintidx = p_network->getLocalBusIndices(bus_number);
		nbus = vec_busintidx.size();
		for(ibus=0; ibus<nbus; ibus++){
			bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
			//printf("----renke debug scatterInjectionLoad, in dsf full app, \n");
			bus->scatterInjectionLoadConstCurrent(vCurR[ival], vCurI[ival]);
		
		}
	}
		
}
/**
 * execute load shedding	 
 */
void gridpack::dynamic_simulation::DSFullApp::applyLoadShedding(int bus_number, std::string loadid, double percentage){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		bus->applyLoadShedding(loadid, percentage);
	
	}
		
}

/**
 * set the wide area control signals of the PSS of a certain generator
 * input bus_number: generator bus number
 * input bus_number: generator gen ID
 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
 */
void gridpack::dynamic_simulation::DSFullApp::setWideAreaControlSignal(int bus_number, std::string genid, double wideAreaControlSignal){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug GFI Par adjustment, in dsf full app, \n");
		bus->setWideAreaControlSignal(genid, wideAreaControlSignal);
	
	}
	
}

/**
 * execute Grid Forming Inverter control parameters adjustment
 * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid
 * input bus_number: GFI bus number
 * input bus_number: GFI gen ID
 * input newParValScaletoOrg:  GFI new parameter scale value to the very intial value at the begining of dynamic simulation
 */
void gridpack::dynamic_simulation::DSFullApp::applyGFIAdjustment(int controlType, int bus_number, std::string genid, double newParValScaletoOrg){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug GFI Par adjustment, in dsf full app, \n");
		bus->applyGFIAdjustment(controlType, genid, newParValScaletoOrg);
	
	}
	
}

/**
 * execute constant Y load shedding at a curtain bus	 
 * bus number
 * percentage: float load shed percentage, for example -0.2 means shed 20%
 */
void gridpack::dynamic_simulation::DSFullApp::applyConstYLoadShedding(int bus_number, double percentage ){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	//printf("----renke debug applyConstYLoadShedding, in dsf full app, bus_number: %d, nbus: %d\n", bus_number, nbus);
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug applyConstYLoadShedding, in dsf full app, percentage: %f \n", percentage);
		bus->applyConstYLoadShedding(percentage);
	}		
}


/**
 * execute constant Y load P change	 
 */
void gridpack::dynamic_simulation::DSFullApp::applyConstYLoad_Change_P(int bus_number, double loadPChangeMW ){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		bus->applyConstYLoad_Change_P(loadPChangeMW);
		bapplyLoadChangeP = true;
		p_vbus_need_to_changeP.push_back(bus);
	}
	
   // Check to see if a load change P is occuring somewhere in the system
   bapplyLoadChangeP = p_factory->checkTrueSomewhere(bapplyLoadChangeP);
		
}

void gridpack::dynamic_simulation::DSFullApp::clearConstYLoad_Change_P()
{
	//clear the flags and the tripping line pointer vector, as well as the flags in the tripping line
	bapplyLoadChangeP = false;
	int ibus, nbus;
	nbus = p_vbus_need_to_changeP.size();
	for (ibus=0 ; ibus<nbus; ibus++){
		p_vbus_need_to_changeP[ibus]->clearConstYLoad_Change_P();
	}
	p_vbus_need_to_changeP.clear();	
}

/**
 * set constant Y load P to zero for a specific bus	 
 */
void gridpack::dynamic_simulation::DSFullApp::setConstYLoadtoZero_P(int bus_number){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	bool ret;
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		ret = bus->setConstYLoadtoZero_P( );
		if (ret){
			bapplyLoadChangeP = true;
			p_vbus_need_to_changeP.push_back(bus);
		}
	}
	
   // Check to see if a load change P is occuring somewhere in the system
   bapplyLoadChangeP = p_factory->checkTrueSomewhere(bapplyLoadChangeP);
		
}

/**
 * set constant Y load to impedancer and impedancei
 */
void gridpack::dynamic_simulation::DSFullApp::setConstYLoadImpedance(int bus_number, double impedancer, double impedancei){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	bool ret;
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		ret = bus->setConstYLoadtoValue(impedancer, impedancei);
		if (ret){
			bapplyLoadChangeP = true;
			p_vbus_need_to_changeP.push_back(bus);
			bapplyLoadChangeQ = true;
			p_vbus_need_to_changeQ.push_back(bus);
		}				
	}
	
   // Check to see if a load change P is occuring somewhere in the system
   bapplyLoadChangeP = p_factory->checkTrueSomewhere(bapplyLoadChangeP);
   bapplyLoadChangeQ = p_factory->checkTrueSomewhere(bapplyLoadChangeQ);
		
}

/**
 * execute constant Y load Q change	 
 */
void gridpack::dynamic_simulation::DSFullApp::applyConstYLoad_Change_Q(int bus_number, double loadPChangeMVAR ){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		bus->applyConstYLoad_Change_Q(loadPChangeMVAR);
		bapplyLoadChangeQ = true;
		p_vbus_need_to_changeQ.push_back(bus);
	}
	
	// Check to see if a load change Q is occuring somewhere in the system
   bapplyLoadChangeQ = p_factory->checkTrueSomewhere(bapplyLoadChangeQ);
		
}

void gridpack::dynamic_simulation::DSFullApp::clearConstYLoad_Change_Q()
{
	//clear the flags and the tripping line pointer vector, as well as the flags in the tripping line
	bapplyLoadChangeQ = false;
	int ibus, nbus;
	nbus = p_vbus_need_to_changeQ.size();
	for (ibus=0 ; ibus<nbus; ibus++){
		p_vbus_need_to_changeQ[ibus]->clearConstYLoad_Change_Q();
	}
	p_vbus_need_to_changeQ.clear();	
}

/**
 * set constant Y load Q to zero for a specific bus	 
 */
void gridpack::dynamic_simulation::DSFullApp::setConstYLoadtoZero_Q(int bus_number){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	bool ret;
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		ret = bus->setConstYLoadtoZero_Q();
		if (ret){
			bapplyLoadChangeQ = true;
			p_vbus_need_to_changeQ.push_back(bus);
		}
	}
	
	// Check to see if a load change Q is occuring somewhere in the system
   bapplyLoadChangeQ = p_factory->checkTrueSomewhere(bapplyLoadChangeQ);
		
}


/**
 * execute generator tripping
 */
void gridpack::dynamic_simulation::DSFullApp::applyGeneratorTripping(int bus_number, std::string genid){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug generator trip, in dsf full app, \n");
		bus->applyGeneratorTripping(genid);
	
	}
		
}

/**
 * set all the necessery flags for the two buses and one branch for the line needs to trip
 * this function is for single branch flags set-up, may need to be called
 * multiple times for multiple line tripping 
 */
void gridpack::dynamic_simulation::DSFullApp::setLineTripAction
(int brch_from_bus_number, int brch_to_bus_number, std::string branch_ckt){
	
	gridpack::utility::StringUtils util;
	std::string clean_brkckt;
	std::vector<int> vec_branchintidx;
	vec_branchintidx = p_network->getLocalBranchIndices(brch_from_bus_number, brch_to_bus_number);
	int ibr, nbr;
	gridpack::dynamic_simulation::DSFullBranch *pbranch;	
	nbr = vec_branchintidx.size();
	//printf("----renke debug load shed, in dsf full app::setLineTripAction, there are %d branches from bus %d to bus %d \n", nbr, brch_from_bus_number, brch_to_bus_number);
	for(ibr=0; ibr<nbr; ibr++){
		pbranch = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>
			(p_network->getBranch(vec_branchintidx[ibr]).get());
		//printf("----renke debug load shed, in dsf full app::setLineTripAction, from bus:, to bus:\n");
		clean_brkckt = util.clean2Char(branch_ckt);
		
		if (pbranch->setBranchTripAction(clean_brkckt)){
			p_vbranches_need_to_trip.push_back(pbranch);
         bapplyLineTripAction = true;
			break;		
		}	
	}			
   // Check to see if a line trip is occuring somewhere in the system
   bapplyLineTripAction = p_factory->checkTrueSomewhere(bapplyLineTripAction);
}

// trip a branch, given a bus number, just find any one of the connected line(not transformer) with the bus, and trip that one
void gridpack::dynamic_simulation::DSFullApp::setLineTripAction(int bus_number){
	
	std::vector<int> vec_busintidx;
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	int ibus, nbus;
	gridpack::dynamic_simulation::DSFullBus *bus;	
	nbus = vec_busintidx.size();
	std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > vec_nghbrs;
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
			(p_network->getBus(vec_busintidx[ibus]).get());
		//printf("----renke debug load shed, in dsf full app, \n");
		
		bus->getNeighborBranches(vec_nghbrs);
		
		int ibr, nbr;
		nbr = vec_nghbrs.size();
		gridpack::dynamic_simulation::DSFullBranch *pbranch;
		for (ibr=0; ibr<nbr; ibr++){
			
			pbranch = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>
			(vec_nghbrs[ibr].get());
			
			if (pbranch->setBranchTripAction()){ // if the branch is a non-xmfr branch
				p_vbranches_need_to_trip.push_back(pbranch);
            bapplyLineTripAction = true;
				break;
				
			}					
		}	
	}				
   // Check to see if a line trip is occuring somewhere in the system
   bapplyLineTripAction = p_factory->checkTrueSomewhere(bapplyLineTripAction);
}

/**
 * clear all the necessery flags for the all buses and branches for the lines needs to trip
 * this function is for all the branches' flags clear-up, just need to be called
 * once to clear all the tripping lines's flag
 */
void gridpack::dynamic_simulation::DSFullApp::clearLineTripAction()
{
	//clear the flags and the tripping line pointer vector, as well as the flags in the tripping line
	bapplyLineTripAction = false;
	int nbr, ibr;
	nbr = p_vbranches_need_to_trip.size();
	for (ibr=0 ; ibr<nbr; ibr++){
		p_vbranches_need_to_trip[ibr]->clearBranchTripAction();
	}
	p_vbranches_need_to_trip.clear();	
}

/**
 * Return values for total active and reactive load power on bus
 * @param bus_id original bus index
 * @param lp active load power
 * @param lq reactive load power
 * @return false if bus is not found on this processor
 */
bool gridpack::dynamic_simulation::DSFullApp::getBusTotalLoadPower(int bus_id,
    double &total_p, double &total_q)
{
  std::vector<int> indices = p_network->getLocalBusIndices(bus_id);
  int i;
  for (i=0; i<indices.size(); i++) {
    if (p_network->getActiveBus(indices[i])) {
      p_network->getBus(indices[i])->getTotalLoadPower(total_p, total_q);
      return true;
    }
  }
  if (p_report_dummy_obs){
	  total_p = 0.0;
	  total_q = 0.0;
	  return true;
	
  }else{
	 return false; 
  }
}

/**
 * Return real and reactive power produced by requested generator
 * @param bus_id original index for bus hosting generator
 * @param gen_id 2-character identifier for generator
 * @param pg active power produced by generator
 * @param qg reactive power produced by generator
 * @return false if generator is not found on this processor
 */
bool gridpack::dynamic_simulation::DSFullApp::getGeneratorPower(int bus_id,
    std::string gen_id, double &pg, double &qg)
{
  std::vector<int> indices = p_network->getLocalBusIndices(bus_id);
  int i;
  gridpack::utility::StringUtils util;
  std::string clean_genid;
  clean_genid = util.clean2Char(gen_id);
  
  for (i=0; i<indices.size(); i++) {
    if (p_network->getActiveBus(indices[i])) {
      if (p_network->getBus(indices[i])->getGeneratorPower(clean_genid, pg, qg)) {
        return true;
      } else {
			if (p_report_dummy_obs){
				pg = 0.0;
				qg = 0.0;
				return true;
			}else{
				return false; 
			} 
        //return false;
      }
    }
  }
  if (p_report_dummy_obs){
	  pg = 0.0;
	  qg = 0.0;
	  return true;
	
  }else{
	 return false; 
  }
  //return false;
}

/**
 * Get a list of unique zones in the system
 * @param zones a complete list of all zones in the network
 */
void gridpack::dynamic_simulation::DSFullApp::getZoneList(
    std::vector<int> &zones) const
{
  int nbus = p_network->numBuses();
  int i, j, idx;
  // find out how many zones are in system and provide a list of zones to each
  // processor. Start by getting list of zones on this processor
  std::map<int,int> zmap;
  idx = 0;
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
      int iz = bus->getZone();
      if (zmap.find(iz) == zmap.end()) {
        zmap.insert(std::pair<int,int>(iz,idx));
        idx++;
      }
    }
  }
  // exchange zones with each processor. Start by distributing number of zones on
  // each processor
  int nsize = zmap.size();
  int nproc = p_network->communicator().size();
  int me = p_network->communicator().rank();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) {
    sizes[i] = 0;
  }
  sizes[me] = nsize;
  // exchange number of zones found on each process 
  p_network->communicator().sum(&sizes[0],nproc);
  std::vector<int> start(nproc);
  start[0] = 0;
  int ntot = sizes[0];
  // calculate offsets
  for (i=1; i<nproc; i++) {
    start[i] = start[i-1]+sizes[i-1];
    ntot += sizes[i];
  }
  std::vector<int> allzones(ntot); 
  for (i=0; i<ntot; i++) {
    allzones[i] = 0;
  }
  std::map<int,int>::iterator it = zmap.begin();
  idx = 0;
  while (it != zmap.end()) {
    allzones[start[me]+idx] = it->first;
    it++;
    idx++;
  }
  // exchange zones
  p_network->communicator().sum(&allzones[0],ntot);
  // create list of unique zones
  zmap.clear();
  idx = 0;
  for (i=0; i<ntot; i++) {
    if (zmap.find(allzones[i]) == zmap.end()) {
      zmap.insert(std::pair<int,int>(allzones[i],idx));
      idx++;
    }
  }
  int nzones = zmap.size();
  zones.resize(nzones);
  it = zmap.begin();
  while (it != zmap.end()) {
    if (me == 0) {
      //printf("NZONES: %d idx: %d zoneID: %d\n",nzones,it->second,it->first);
    }
    zones[it->second] = it->first;
    it++;
  }
  // sort zones in ascending order (use simple, but inefficient, bubble sort.
  // Assume zones list is not long)
  for (i=0; i<nzones; i++) {
    for (j=i+1; j<nzones; j++) {
      if (zones[i] > zones[j]) {
        int itmp = zones[i];
        zones[i] = zones[j];
        zones[j] = itmp;
      }
    }
  }
}


/**
 * Return total active and reactive loads for each zone
 * @param load_p active load for all zones
 * @param load_q reactive load for all zones
 * @param zone_id label for each zone
 */
void gridpack::dynamic_simulation::DSFullApp::getZoneLoads(std::vector<double> &load_p,
    std::vector<double> &load_q, std::vector<int> &zone_id) const
{
  int nbus = p_network->numBuses();
  int i, idx;

  zone_id.clear();
  getZoneList(zone_id);
  // create map of zones
  std::map<int,int> zones;
  int nzones = zone_id.size();
  for (i=0; i<nzones; i++) {
    zones.insert(std::pair<int,int>(zone_id[i],i));
  }

  load_p.clear();
  load_q.clear();
  load_p.resize(nzones);
  load_q.resize(nzones);
  for (i=0; i<nzones; i++) {
    load_p[i] = 0.0;
    load_q[i] = 0.0;
  }
  // get total loads
  std::map<int,int>::iterator it;
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
      int iz = bus->getZone();
      double p, q;
      bus->getTotalLoadPower(p,q);
      it = zones.find(iz);
      if (it != zones.end()) {
        idx = it->second;
        load_p[idx] += p;
        load_q[idx] += q;
      } else {
        printf("Unknown zone on bus %d: %d\n",bus->getOriginalIndex(),iz);
      }
    }
  }
  p_network->communicator().sum(&load_p[0],nzones);
  p_network->communicator().sum(&load_q[0],nzones);
}

/**
 * Return total active and reactive generator power for each zone
 * @param generator_p active generator power for all zones
 * @param generator_q reactive generator power for all zones
 * @param zone_id label for each zone
 */
void gridpack::dynamic_simulation::DSFullApp::getZoneGeneratorPower(
    std::vector<double> &generator_p, std::vector<double> &generator_q,
    std::vector<int> &zone_id) const
{
  int nbus = p_network->numBuses();
  int i, idx;
  int me = p_network->communicator().rank();

  zone_id.clear();
  getZoneList(zone_id);
  // create map of zones
  std::map<int,int> zones;
  int nzones = zone_id.size();
  for (i=0; i<nzones; i++) {
    zones.insert(std::pair<int,int>(zone_id[i],i));
  }

  generator_p.clear();
  generator_q.clear();
  generator_p.resize(nzones);
  generator_q.resize(nzones);
  for (i=0; i<nzones; i++) {
    generator_p[i] = 0.0;
    generator_q[i] = 0.0;
  }

  // get total generator power
  std::map<int,int>::iterator it;
  for (i=0; i<nbus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
      int iz = bus->getZone();
      double p, q;
      bus->getTotalGeneratorPower(p,q);
      it = zones.find(iz);
      if (it != zones.end()) {
        idx = it->second;
        generator_p[idx] += p;
        generator_q[idx] += q;
      } else {
        printf("Invalid zone on bus %d: %d\n",bus->getOriginalIndex(),iz+1);
      }
    }
  }
  p_network->communicator().sum(&generator_p[0],nzones);
  p_network->communicator().sum(&generator_q[0],nzones);
}

/**
 * Modify generator parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param gen_id two character token specifying generator on bus
 * @param genParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, double value)
{
  gridpack::utility::StringUtils util;
  std::string clean_gen_id = util.clean2Char(gen_id);
  return p_modifyDataCollectionGenParam<double>(bus_id,clean_gen_id,genParam,value);
}
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, int value)
{
  gridpack::utility::StringUtils util;
  std::string clean_gen_id = util.clean2Char(gen_id);
  return p_modifyDataCollectionGenParam<int>(bus_id,clean_gen_id,genParam,value);
}

/**
 * Modify load parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param load_id two character token specifying load on bus
 * @param loadParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, double value)
{
  gridpack::utility::StringUtils util;
  std::string clean_load_id = util.clean2Char(load_id);
  return p_modifyDataCollectionLoadParam<double>(bus_id,clean_load_id,loadParam,value);
}
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, int value)
{
  gridpack::utility::StringUtils util;
  std::string clean_load_id = util.clean2Char(load_id);
  return p_modifyDataCollectionLoadParam<int>(bus_id,clean_load_id,loadParam,value);
}

/**
 * Modify parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param busParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, double value)
{
  return p_modifyDataCollectionBusParam<double>(bus_id,busParam,value);
}
bool gridpack::dynamic_simulation::DSFullApp::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, int value)
{
  return p_modifyDataCollectionBusParam<int>(bus_id,busParam,value);
}

/**
 * Set the state of some device on the network
 * @param bus_id bus ID
 * @param dev_id two character identifier of device
 * @param type of device (GENERATOR, EXCITER, OR GOVERNOR)
 * @param name of the state variable (ANGLE,SPEED_DEV)
 * @param value new value of parameter
 * @return false if this device or
 parameter not found
 */
bool gridpack::dynamic_simulation::DSFullApp::setState(int bus_id,
    std::string dev_id, std::string device, std::string name, double value)
{
  gridpack::utility::StringUtils util;
  bool ret = false;
  util.toUpper(device);
  util.trim(device);
  util.trim(dev_id);
  util.trim(name);
  dev_id = util.clean2Char(dev_id);
  /* return false if ther is no buse corresponding to this bus ID */
  std::vector<int> buses = p_network->getLocalBusIndices(bus_id);
  if (buses.size() > 0) {
    ret = true;
  }
  if (!p_factory->checkTrueSomewhere(ret)) {
    return false;
  }
  ret = false;
  /* search for device corresponding to this dev_id and device type. */
  int i;
  for (i=0; i<buses.size(); i++) {
    DSFullBus *bus = p_network->getBus(buses[i]).get();
    if (device == "GENERATOR"  || device == "GOVERNOR" || device == "EXCITER") {
      BaseGeneratorModel *generator = bus->getGenerator(dev_id);
      if (generator != NULL) {
        if (device == "GOVERNOR") {
          BaseGovernorModel *governor = generator->getGovernor().get();
          if (governor != NULL) {
            ret = governor->setState(name, value);
          }
          if (!ret) {
            std::cout<<"No parameter "<<name<<" found on governor"<<std::endl;
          }
        } else if (device == "EXCITER") {
          BaseExciterModel *exciter = generator->getExciter().get();
          if (exciter != NULL) {
            ret = exciter->setState(name, value);
          }
          if (!ret) {
            std::cout<<"No parameter "<<name<<" found on exciter"<<std::endl;
          }
        } else {
          ret = generator->setState(name, value);
          if (!ret) {
            std::cout<<"No parameter "<<name<<" found on governor"<<std::endl;
          }
        }
      } else {
        ret = false;
        if (!ret) {
          std::cout<<"No generator "<<dev_id<<" found on bus "<<bus_id<<std::endl;
        }
      }
    } else if (device == "RELAY") {
      std::cout<<"No setState function implemented for relays "<<std::endl;
    } else if (device == "LOAD") {
      std::cout<<"No setState function implemented for loads "<<std::endl;
    } else if (device == "PLANT") {
      std::cout<<"No setState function implemented for plants "<<std::endl;
    } else if (device == "PSS") {
      std::cout<<"No setState function implemented for PPS "<<std::endl;
    } else {
      std::cout<<"Device not found on bus "<<bus_id<<
        " for device "<<dev_id<<std::endl;
    } 
  }
  return ret;
}

/**
 * Get the state of some device on the network
 * @param  bus_id bus ID
 * @param  dev_id two character identifier of device
 * @param  type of device (GENERATOR, EXCITER, OR GOVERNOR)
 * @param  name of the state variable (ANGLE, SPEED_DEV)
 * @return value current value of parameter
 * @return false if this device or state variable not found
 */
bool gridpack::dynamic_simulation::DSFullApp::getState(int bus_id,
    std::string dev_id, std::string device, std::string name, double *value)
{
  gridpack::utility::StringUtils util;
  bool ret = false;
  util.toUpper(device);
  util.trim(device);
  util.trim(dev_id);
  util.trim(name);
  dev_id = util.clean2Char(dev_id);
  /* return false if ther is no buse corresponding to this bus ID */
  std::vector<int> buses = p_network->getLocalBusIndices(bus_id);
  if (buses.size() > 0) {
    ret = true;
  }
  if (!p_factory->checkTrueSomewhere(ret)) {
    return false;
  }
  ret = false;
  /* search for device corresponding to this dev_id and device type. */
  int i;
  for (i=0; i<buses.size(); i++) {
    if (p_network->getActiveBus(buses[i])) {
      DSFullBus *bus = p_network->getBus(buses[i]).get();
      if (device == "GENERATOR"  || device == "GOVERNOR" || device == "EXCITER") {
        BaseGeneratorModel *generator = bus->getGenerator(dev_id);
        if (generator != NULL) {
          if (device == "GOVERNOR") {
            BaseGovernorModel *governor = generator->getGovernor().get();
            if (governor != NULL) {
              ret = governor->getState(name, value);
            }
            if (!ret) {
              std::cout<<"No parameter "<<name<<" found on governor"<<std::endl;
            }
          } else if (device == "EXCITER") {
            BaseExciterModel *exciter = generator->getExciter().get();
            if (exciter != NULL) {
              ret = exciter->getState(name, value);
            }
            if (!ret) {
              std::cout<<"No parameter "<<name<<" found on exciter"<<std::endl;
            }
          } else {
            ret = generator->getState(name, value);
            if (!ret) {
              std::cout<<"No parameter "<<name<<" found on governor"<<std::endl;
            }
          }
        } else {
          ret = false;
          if (!ret) {
            std::cout<<"No generator "<<dev_id<<" found on bus "<<bus_id<<std::endl;
          }
        }
      } else if (device == "RELAY") {
        std::cout<<"No getState function implemented for relays "<<std::endl;
      } else if (device == "LOAD") {
        std::cout<<"No getState function implemented for loads "<<std::endl;
      } else if (device == "PLANT") {
        std::cout<<"No getState function implemented for plants "<<std::endl;
      } else if (device == "PSS") {
        std::cout<<"No getState function implemented for PPS "<<std::endl;
      } else {
        std::cout<<"Device not found on bus "<<bus_id<<
          " for device "<<dev_id<<std::endl;
      } 
    }
  }
  return ret;
}

/**
 * Transfer data from power flow to dynamic simulation
 * @param pf_network power flow network
 * @param ds_network dynamic simulation network
 */
void gridpack::dynamic_simulation::DSFullApp::transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
    ds_network)
{
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *pfData;
  gridpack::component::DataCollection *dsData;
  double rval;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    dsData = ds_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    dsData->setValue("BUS_PF_VMAG", rval);
    dsData->setValue(BUS_VOLTAGE_MAG,rval);
    ///printf("Step0 bus%d mag = %f\n", i+1, rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue("BUS_PF_VANG",rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);

    pfData->getValue("LOAD_PL", &rval, 0);
    dsData->setValue(LOAD_PL, rval);
    pfData->getValue("LOAD_QL", &rval, 0);
    dsData->setValue(LOAD_QL, rval);
    
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        dsData->setValue(GENERATOR_PG,rval,j);
        //printf("save PGEN: %f\n", rval);
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
        //printf("save QGEN: %f\n", rval);
      }
    }
  }
}
