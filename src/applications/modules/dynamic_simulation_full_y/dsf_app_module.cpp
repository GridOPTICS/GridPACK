/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shuangshuang Jin
 * @Last modified:  May 13, 2015
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------
//
#define USE_TIMESTAMP

#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/parallel/global_vector.hpp"
#include "dsf_app_module.hpp"
#include <iostream>
#include <string>
#include <vector>

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
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSFullApp::~DSFullApp(void)
{
}

enum Format{PTI23, PTI33};
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

  // Monitor generators for frequency violations
  p_monitorGenerators = cursor->get("monitorGenerators",false);
  p_maximumFrequency = cursor->get("frequencyMaximum",61.8);

  // load input file
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<DSFullNetwork> parser(network);
    if (filename.size() > 0) parser.parse(filename.c_str());
  } else {
    printf("Unknown filetype\n");
  }
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  filename = cursor->get("generatorParameters","");

  // partition network
  network->partition();

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

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(512, network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<DSFullNetwork>(128, network));
}

/**
 * Read generator parameters. These will come from a separate file (most
 * likely). The name of this file comes from the input configuration file.
 */
void gridpack::dynamic_simulation::DSFullApp::readGenerators(void)
{
  int rank = p_network->communicator().rank();
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  gridpack::parser::PTI23_parser<DSFullNetwork> parser(p_network);
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  std::string filename = cursor->get("generatorParameters","");
  printf("p[%d] generatorParameters: %s\n",p_comm.rank(),filename.c_str());
  if (filename.size() > 0) parser.externalParse(filename.c_str());
  printf("p[%d] finished Generator parameters\n",p_comm.rank());
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
  int t_solve = timer->createCategory("DS Solve: Total");
  int t_misc = timer->createCategory("DS Solve: Miscellaneous");
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif
  timer->start(t_solve);
  timer->start(t_misc);

  // Get cursor for setting solver options
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  timer->stop(t_misc);

  int t_mode = timer->createCategory("DS Solve: Set Mode");
  timer->start(t_mode);
  p_factory->setMode(YBUS);
  timer->stop(t_mode);
  int t_ybus = timer->createCategory("DS Solve: Make YBus");
  timer->start(t_ybus);
  gridpack::mapper::FullMatrixMap<DSFullNetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  
  //printf("\n=== org ybus: ============\n");
  //orgYbus->print();
  //orgYbus->save("ybus_GridPACK_org.m");
  //exit(0);

  //p_factory->addLoadAdmittance();

  // Form constant impedance load admittance yl for all buses and add it to
  // system Y matrix: ybus = ybus + yl
  p_factory->setMode(YL);
  boost::shared_ptr<gridpack::math::Matrix> ybusyl = ybusMap.mapToMatrix();
  timer->stop(t_ybus);
  //branchIO.header("\n=== ybus after added yl: ============\n");
  //printf("\n=== ybus after added yl: ============\n");
  //ybusyl->print();
  //ybusyl->save("ybus_GridPACK_yl.m");
  //exit(0);

  p_factory->setMode(PG);
  boost::shared_ptr<gridpack::math::Matrix> ybuspg = ybusMap.mapToMatrix();
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
  boost::shared_ptr<gridpack::math::Matrix> ybus_jxd = ybusMap.mapToMatrix();
  //branchIO.header("\n=== ybusyl after added j*Xd': =============\n");
  //printf("\n=== ybusyl after added j*Xd': =============\n");
  //ybus_jxd->print();
  //ybus_jxd->save("ybus_GridPACK_jxd.m");
  
  
  // Add dynamic load impedance to system Y matrix:
  timer->start(t_mode);
  p_factory->setMode(YDYNLOAD);
  timer->stop(t_mode);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  //branchIO.header("\n=== ybus_jxd after added dynamic load impedance': =============\n");
  //printf("\n=== ybus_dynload after added dynamic load impedance': =============\n");
  //ybus->print();
  //ybus->save("ybus_GridPACK_dynload.m");
  
  //exit(0);
 
  // Get fault information from fautlts Event from input.xml
  int sw2_2 = fault.from_idx - 1;
  int sw3_2 = fault.to_idx - 1;

  // Compute ybus_fy for fault on stage
  boost::shared_ptr<gridpack::math::Matrix> ybus_fy(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_misc);
  p_factory->setEvent(fault);
  timer->stop(t_misc);
  timer->start(t_mode);
  p_factory->setMode(onFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap.overwriteMatrix(ybus_fy);
  //branchIO.header("\n=== ybus_fy: ============\n");
  //printf("\n=== ybus_fy: ============\n");
  //ybus_fy->print();
  //ybus_fy->save("ybus_fy_GridPACK_jxd.m");

  // Compute ybus_posfy for fault clear stage
  boost::shared_ptr<gridpack::math::Matrix> ybus_posfy(ybus->clone());
  timer->stop(t_ybus);
  timer->start(t_mode);
  p_factory->setMode(posFY);
  timer->stop(t_mode);
  timer->start(t_ybus);
  ybusMap.incrementMatrix(ybus_posfy);
  //branchIO.header("\n=== ybus_posfy: ============\n");
  //printf("\n=== ybus_posfy: ============\n");
  //ybus_posfy->print();
  //ybus_posfy->save("ybus_posfy_GridPACK_jxd.m");
  timer->stop(t_ybus);

  // Simulation related variables
  int t_init = timer->createCategory("DS Solve: Initialization");
  timer->start(t_init);
  int simu_k;
  int t_step[20];
  double t_width[20];
  int S_Steps;
  int last_S_Steps;
  int steps3, steps2, steps1;
  double h_sol1, h_sol2;
  int flagP, flagC;
  int I_Steps;

  const double sysFreq = 60.0;
  double pi = 4.0*atan(1.0);
  const double basrad = 2.0 * pi * sysFreq;
  gridpack::ComplexType jay(0.0, 1.0);

  // switch info is set up here
  int nswtch = 4;
  static double sw1[4];
  static double sw7[4];
  sw1[0] = 0.0;
  sw1[1] = fault.start;
  sw1[2] = fault.end;
  sw1[3] = p_sim_time;
  sw7[0] = p_time_step;
  sw7[1] = fault.step;
  sw7[2] = p_time_step;
  sw7[3] = p_time_step;
  simu_k = 0;
  for (int i = 0; i < nswtch-1; i++) {
    t_step[i] = (int) ((sw1[i+1] -sw1[i]) / sw7[i]);
    t_width[i] = (sw1[i+1] - sw1[i]) / t_step[i];
    simu_k += t_step[i];
  }
  simu_k++;
  
  // Initialize vectors for integration 
  p_factory->initDSVect(p_time_step);
  //exit(0);

  gridpack::mapper::BusVectorMap<DSFullNetwork> ngenMap(p_network);
  // Map to create vector volt
  boost::shared_ptr<gridpack::math::Vector> volt = ngenMap.mapToVector();
  //p_busIO->header("\n=== volt: ===\n");
  //volt->print();

  gridpack::math::LinearSolver solver(*ybus);
  solver.configure(cursor);
  gridpack::math::LinearSolver solver_fy(*ybus_fy);
  solver_fy.configure(cursor);
  //gridpack::math::LinearSolver solver_posfy(*ybus_posfy);
  gridpack::math::LinearSolver solver_posfy(*ybus);
  solver_posfy.configure(cursor);

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagP = 0;
  flagC = 0;
  S_Steps = 1;
  last_S_Steps = -1;

  bool flag = true, flag_corrector = true;

  p_insecureAt = -1;

  p_factory->setMode(make_INorton_full);
  gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> INorton_full = nbusMap.mapToVector();
  boost::shared_ptr<gridpack::math::Vector> INorton_full_chk = nbusMap.mapToVector();
  double max_INorton_full = 0.0;
  boost::shared_ptr<gridpack::math::Vector> volt_full(INorton_full->clone());

  timer->stop(t_init);
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


  for (I_Steps = 0; I_Steps < simu_k - 1; I_Steps++) {
  //for (I_Steps = 0; I_Steps < 200; I_Steps++) {
    //char step_str[128];
    //sprintf(step_str,"\nIter %d\n", I_Steps);
    //p_busIO->header(step_str);
    timer->start(t_misc);
    //printf("Step %d\ttime %5.3f sec: \n", I_Steps+1, (I_Steps+1) * p_time_step);
    //printf("\n===================Step %d\ttime %5.3f sec:================\n", I_Steps+1, (I_Steps+1) * p_time_step);
    ///char step_str[128];
    ///sprintf(step_str, "\n===================Step %d\ttime %5.3f sec:================\n", I_Steps+1, (I_Steps+1) * p_time_step);
     ///p_busIO->header(step_str);
    S_Steps = I_Steps;

    if (I_Steps < steps1) {
      flagP = 0;
      flagC = 0;
    } else if (I_Steps == steps1) {
      flagP = 0;
      //flagC = 1;
      flagC = 0;
    } else if ((I_Steps > steps1) && (I_Steps < steps2)) {
      flagP = 1;
      flagC = 1;
    } else if (I_Steps == steps2) {
      flagP = 1;
      //flagC = 2;
      flagC = 1;
    } else if (I_Steps > steps2) {
      flagP = 2;
      flagC = 2;
    }
    timer->stop(t_misc);
    
    if (I_Steps !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor_currentInjection(false);
    } else {
      p_factory->predictor_currentInjection(true);
    }

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    int t_mIf = timer->createCategory("DS Solve: Modified Euler Predictor: Make INorton");
    timer->start(t_mIf);
	p_factory->setMode(make_INorton_full);
    nbusMap.mapToVector(INorton_full);
    ///gridpack::mapper::BusVectorMap<DSFullNetwork> nbusMap(p_network);
    ///boost::shared_ptr<gridpack::math::Vector> INorton_full = nbusMap.mapToVector();
    //p_busIO->header("\n=== [Predictor] INorton_full: ===\n");
    //printf("renke test \n=== [Predictor] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_mIf);
#ifdef MAP_PROFILE
  timer->configTimer(false);
#endif
 
    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2) 
    // to calculate terminal volt: ----------
    int t_psolve = timer->createCategory("DS Solve: Modified Euler Predictor: Linear Solver");
    timer->start(t_psolve);
    //boost::shared_ptr<gridpack::math::Vector> volt_full(INorton_full->clone());
    volt_full->zero();
#if 0
    bool flag_chk = true;
    while (flag_chk == true ) {
		
			volt_full->zero();
			
			if (flagP == 0) {
				solver.solve(*INorton_full, *volt_full);
			} else if (flagP == 1) {
				solver_fy.solve(*INorton_full, *volt_full);
			} else if (flagP == 2) {
				solver_posfy.solve(*INorton_full, *volt_full);
			}
			

			printf("1: itr test:----previous predictor_INorton_full:\n");
			INorton_full->print();

			INorton_full_chk->equate(*INorton_full);
			printf("2: itr test:----predictor_INorton_full_chk:\n");
			INorton_full_chk->print();

			nbusMap.mapToBus(volt_full);
			p_factory->setVolt(false);
			
			if (I_Steps !=0 && last_S_Steps != S_Steps) {
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
			nbusMap.mapToVector(INorton_full);
			
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
      solver.solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy.solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy.solve(*INorton_full, *volt_full);
    }
#endif
    timer->stop(t_psolve);

#ifdef MAP_PROFILE
  timer->configTimer(true);
#endif
    //p_busIO->header("\n=== [Predictor] volt_full: ===\n");
    //volt_full->print();
    //if (I_Steps==4){
    //	 exit(0);
   //	}

    int t_vmap= timer->createCategory("DS Solve: Map Volt to Bus");
    timer->start(t_vmap);
	
	//printf("after first volt sovle, before first volt map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap.mapToBus(volt_full);
	
	//printf("after first volt sovle, after first volt map: \n");
	
	if ( I_Steps==0 ) {
		//printf("enter the initial update oldbusvoltage, Timestep: %d \n", I_Steps);
		p_factory->updateoldbusvoltage(); //renke add, first timestep, copy volt_full to volt_full_old
	}
    timer->stop(t_vmap);

    int t_volt= timer->createCategory("DS Solve: Set Volt");
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

	 helics_requestTime =       double (I_Steps*h_sol1);
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
    //printf("Timestep, %d \n", I_Steps);
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
             ybusMap.overwriteMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(bus_relay);
             ybusMap.overwriteMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",I_Steps );
			 
	     ybus_fy->save(sybus);
	 
	     printf("DSFull_APP::Solve: bus relay trip during fault, ybus changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap.overwriteMatrix(ybus);
	     ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",I_Steps );

	     ybus->save(sybus);
        
             printf("DSFull_APP::Solve: bus relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap.overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",I_Steps );

             ybus_posfy->save(sybus);

			 
        } else if (flagP == 2) {
             p_factory->setMode(bus_relay);
             ybusMap.overwriteMatrix(ybus);
             printf("DSFull_APP::Solve: bus relay trip after fault, ybus changed:\n");
	     ybus->print();
	     char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",I_Steps );

	     ybus->save(sybus);

             printf("DSFull_APP::Solve: bus relay trip after fault, ybus_posfy changed too:\n");
             p_factory->setMode(bus_relay);
             ybusMap.overwriteMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",I_Steps );

             ybus_posfy->save(sybus);

        }
    }
	
	// if branch relay trips, modify the corresponding Ymatrix, renke modified
	if (flagBranch) {
        
        printf("DSFull_APP::Solve: updatebranchrelay return trigger siganl: TURE!!! \n");

        //please update the bus contribution to the Y bus matrix here. //Shuangshuang tbd
	if (flagP == 0) { 
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus);
        } else if (flagP == 1) {
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus_fy);
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_fy changed:\n");
	     ybus_fy->print();
	     char sybus[100];
             sprintf(sybus, "ybus_fy_%d_relay.m",I_Steps );
			 
	     ybus_fy->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus);
             ybus->print();
             sprintf(sybus, "ybus_%d_relay.m",I_Steps );

             ybus->save(sybus);

			 
	     printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus_posfy);
	     ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",I_Steps );

	     ybus_posfy->save(sybus);
			 
        } else if (flagP == 2) {
             printf("DSFull_APP::Solve: branch relay trip during fault, ybus changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus);
             ybus->print();
             char sybus[100];
             sprintf(sybus, "ybus_%d_relay.m",I_Steps );

             ybus->save(sybus);

             printf("DSFull_APP::Solve: branch relay trip during fault, ybus_posfy changed too:\n");
             p_factory->setMode(branch_relay);
             ybusMap.incrementMatrix(ybus_posfy);
             ybus_posfy->print();
             sprintf(sybus, "ybus_posfy_%d_relay.m",I_Steps );

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

    int t_predictor = timer->createCategory("DS Solve: Modified Euler Predictor");
    //printf("Test: predictor begins: \n");
    timer->start(t_predictor);
    if (I_Steps !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor(h_sol1, false);
    } else { 
      p_factory->predictor(h_sol1, true);
    }
    timer->stop(t_predictor);

    if (I_Steps !=0 && last_S_Steps != S_Steps) {
      p_factory->corrector_currentInjection(false);
    } else {
      p_factory->corrector_currentInjection(true);
    }

    //INorton_full = nbusMap.mapToVector();
    int t_cmIf = timer->createCategory("DS Solve: Modified Euler Corrector: Make INorton");
    timer->start(t_cmIf);
    p_factory->setMode(make_INorton_full);
    nbusMap.mapToVector(INorton_full);
    //p_busIO->header("\n=== [Corrector] INorton_full: ===\n");
    //printf("\nrelaytest=== [Corrector] INorton_full: ===\n");
    //INorton_full->print();
    timer->stop(t_cmIf);

    // ---------- CALL ssnetwork_cal_volt(S_Steps+1, flagF2)
    // to calculate terminal volt: ----------
    int t_csolve = timer->createCategory("DS Solve: Modified Euler Corrector: Linear Solver");
    timer->start(t_csolve);
    volt_full->zero();

#if 0
    flag_chk = true;
    while (flag_chk == true ) {
		
			volt_full->zero();
			
			if (flagP == 0) {
				solver.solve(*INorton_full, *volt_full);
			} else if (flagP == 1) {
				solver_fy.solve(*INorton_full, *volt_full);
			} else if (flagP == 2) {
				solver_posfy.solve(*INorton_full, *volt_full);
			}
			nbusMap.mapToBus(volt_full);
			p_factory->setVolt(false);
			
			if (I_Steps !=0 && last_S_Steps != S_Steps) {
				p_factory->corrector_currentInjection(false);
			} else {
				p_factory->corrector_currentInjection(true);
			}
			
			INorton_full_chk->equate(*INorton_full);
			printf("itr test:----corrector_INorton_full_chk:\n");
			INorton_full_chk->print();
			
			p_factory->setMode(make_INorton_full);
			nbusMap.mapToVector(INorton_full);
			
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
      solver.solve(*INorton_full, *volt_full);
    } else if (flagP == 1) {
      solver_fy.solve(*INorton_full, *volt_full);
    } else if (flagP == 2) {
      solver_posfy.solve(*INorton_full, *volt_full);
    }
#endif

    timer->stop(t_csolve);

    //p_busIO->header("\n=== [Corrector] volt_full: ===\n");
    //printf("relaytest \n=== [Corrector] volt_full: ===\n");
    //volt_full->print();
    timer->start(t_vmap);
	
	//printf("after second solve, before second map: \n");
	//p_factory->printallbusvoltage();
	
    nbusMap.mapToBus(volt_full);
	
	//printf("after second solve, after second map: \n");
	//p_factory->printallbusvoltage();
	
    timer->stop(t_vmap);

    timer->start(t_volt);
    p_factory->setVolt(false);
	p_factory->updateBusFreq(h_sol1);
    timer->stop(t_volt);

    int t_corrector = timer->createCategory("DS Solve: Modified Euler Corrector");
    timer->start(t_corrector);
    //printf("Test: corrector begins: \n");
    if (last_S_Steps != S_Steps) {
      p_factory->corrector(h_sol2, false);
    } else {
      p_factory->corrector(h_sol2, true);
    }
    timer->stop(t_corrector);

    //if (I_Steps == simu_k - 1) 
      //p_busIO->write();

    if (I_Steps == steps1) {
      solver_fy.solve(*INorton_full, *volt_full);
//      printf("\n===================Step %d\ttime %5.3f sec:================\n", I_Steps+1, (I_Steps+1) * p_time_step);
//      printf("\n=== [Corrector] volt_full: ===\n");
//      volt_full->print();
      nbusMap.mapToBus(volt_full);
      p_factory->setVolt(false);
	  p_factory->updateBusFreq(h_sol1);
    } else if (I_Steps == steps2) {
      solver_posfy.solve(*INorton_full, *volt_full);
//      printf("\n===================Step %d\ttime %5.3f sec:================\n", I_Steps+1, (I_Steps+1) * p_time_step);
//      printf("\n=== [Corrector] volt_full: ===\n");
//      volt_full->print();
      nbusMap.mapToBus(volt_full);
      p_factory->setVolt(true);
	  p_factory->updateBusFreq(h_sol1);
    }
    if (I_Steps == 1) {
//      printf("\n Dynamic Step 1 [Corrector] volt_full: ===\n");
//      volt_full->print();
//      printf("\n Dynamic Step 1 [Corrector] Norton_full: ===\n");
//      INorton_full->print();
    }
    int t_secure = timer->createCategory("DS Solve: Check Security");
    timer->start(t_secure);
    if (p_generatorWatch && I_Steps%p_generatorWatchFrequency == 0) {
      char tbuf[32];
#ifdef USE_TIMESTAMP
      sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(I_Steps)*p_time_step,
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
      sprintf(tbuf,"%8.4f",static_cast<double>(I_Steps)*p_time_step);
      if (p_generatorWatch) p_generatorIO->header(tbuf);
      if (p_generatorWatch) p_generatorIO->write("watch");
      if (p_generatorWatch) p_generatorIO->header("\n");
#endif
#ifdef USEX_GOSS
      if (p_generatorWatch) p_generatorIO->dumpChannel();
#endif
    }
    if (p_loadWatch && I_Steps%p_loadWatchFrequency == 0) {
      char tbuf[32];
#ifdef USE_TIMESTAMP
      sprintf(tbuf,"%8.4f, %20.4f",static_cast<double>(I_Steps)*p_time_step,
          timer->currentTime());
      if (p_loadWatch) p_loadIO->header(tbuf);
      if (p_loadWatch) p_loadIO->write("load_watch");
      if (p_loadWatch) p_loadIO->header("\n");
#else
      sprintf(tbuf,"%8.4f",static_cast<double>(I_Steps)*p_time_step);
      if (p_loadWatch) p_loadIO->header(tbuf);
      if (p_loadWatch) p_loadIO->write("load_watch");
      if (p_loadWatch) p_loadIO->header("\n");
#endif
#ifdef USEX_GOSS
      if (p_loadWatch) p_loadIO->dumpChannel();
#endif
    }
    saveTimeStep();
    if ((!p_factory->securityCheck()) && p_insecureAt == -1)  
       p_insecureAt = I_Steps;
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
    if (I_Steps == simu_k) {
      printf("\n==============S_Steps = %d==============\n", S_Steps);
      mac_ang_s1->print();
      mac_spd_s1->print();
      p_factory->setMode(init_mac_ang);
      ngenMap.mapToBus(mac_ang_s1);
      p_factory->setMode(init_mac_spd);
      ngenMap.mapToBus(mac_spd_s1);
      p_factory->setMode(init_pmech);
      ngenMap.mapToBus(pmech);
      p_factory->setMode(init_pelect);
      ngenMap.mapToBus(pelect);
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
      p_frequencyOK = p_frequencyOK && checkFrequency(p_maximumFrequency);
      if (!p_frequencyOK) I_Steps = simu_k;
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
    if (fault.isGenerator) {
      sprintf(ptr," for fault at generator %s on bus %d\n",fault.tag,fault.bus_idx);
    } else if (fault.isLine) {
      sprintf(ptr," for fault at line %s from bus %d to bus %d\n",fault.tag,
          fault.from_idx,fault.to_idx);
    } else {
      sprintf(ptr,"!\n");
    }
  } else { 
    char *ptr;
    sprintf(secureBuf,"\nThe system is insecure from step %d", p_insecureAt);
    ptr = secureBuf + strlen(secureBuf);
    if (fault.isGenerator) {
      sprintf(ptr," for fault at generator %s on bus %d\n",fault.tag,fault.bus_idx);
    } else if (fault.isLine) {
      sprintf(ptr," for fault at line %s from bus %d to bus %d\n",fault.tag,
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
 * Read faults from external file and form a list of faults
 * @param cursor pointer to open file contain fault or faults
 * @return a list of fault events
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::dynamic_simulation::DSFullApp::
getFaults(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  std::vector<gridpack::dynamic_simulation::Event> ret;
  if (list) {
    list->children(events);
    int size = events.size();
    int idx;
    // Parse fault events
    for (idx=0; idx<size; idx++) {
      gridpack::dynamic_simulation::Event event;
      event.start = events[idx]->get("beginFault",0.0);
      event.end = events[idx]->get("endFault",0.0);
      std::string indices = events[idx]->get("faultBranch","0 0");
      //Parse indices to get from and to indices of branch
      int ntok1 = indices.find_first_not_of(' ',0);
      int ntok2 = indices.find(' ',ntok1);
      if (ntok2 - ntok1 > 0 && ntok1 != std::string::npos && ntok2 !=
          std::string::npos) {
        event.from_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
        ntok1 = indices.find_first_not_of(' ',ntok2);
        ntok2 = indices.find(' ',ntok1);
        if (ntok1 != std::string::npos && ntok1 < indices.length()) {
          if (ntok2 == std::string::npos) {
            ntok2 = indices.length();
          }
          event.to_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
        } else {
          event.from_idx = 0;
          event.to_idx = 0;
        }
        event.isGenerator = false;
        event.isLine = true;
      } else {
        event.from_idx = 0;
        event.to_idx = 0;
      }
      event.step = events[idx]->get("timeStep",0.0);
      if (event.step != 0.0 && event.end != 0.0 && event.from_idx != event.to_idx) {
        ret.push_back(event);
      }
    }
  }
  return ret;
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
  gridpack::utility::Configuration::ChildCursors generators;
  if (cursor) cursor->children(generators);
  int i, j, idx, id, len;
  int ncnt = generators.size();
  std::string generator, tag, clean_tag;
  gridpack::dynamic_simulation::DSFullBus *bus;
  if (ncnt > 0) p_busIO->header("Monitoring generators:\n");
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
 * @param buses IDs of buses containing generators
 * @param tags generator IDs for watched generators
 * @param writeFile true if external file is to be written
 */
void gridpack::dynamic_simulation::DSFullApp::setGeneratorWatch(
    std::vector<int> &buses, std::vector<std::string> &tags, bool writeFile)
{
  int ncnt = buses.size();
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
    sprintf(buf,"  Bus: %8d Generator ID: %2s\n",id,tag.c_str());
    p_busIO->header(buf);
    if (ncnt > 0) {
      p_generators_read_in = true;
      p_generatorWatch = true;
      sprintf(buf,"Generator Watch Frequency: %d\n",p_generatorWatchFrequency);
      p_busIO->header(buf);
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
    printf("p_gen_buses: %d\n",(int)p_gen_buses.size());
    for (i=0; i<p_gen_buses.size(); i++) {
      std::vector<double> vec0;
      p_time_series.push_back(vec0);
      std::vector<double> vec1;
      p_time_series.push_back(vec1);
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
      std::vector<double> vals = bus->getWatchedValues();
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
            ret.push_back(2*(it->second));
            ret.push_back(2*(it->second)+1);
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
  if (!p_internal_watch_file_name) {
    if (cursor->get("generatorWatchFileName",&filename)) {
      p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(128,
            p_network));
      p_generatorIO->open(filename.c_str());
    } else {
      p_busIO->header("No Generator Watch File Name Found\n");
      p_generatorWatch = false;
    }
  } else {
    p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(128,
          p_network));
    p_generatorIO->open(p_gen_watch_file.c_str());
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
    p_generatorIO->close();
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
  if (cursor->get("loadWatchFileName",&filename)) {
    p_loadIO.reset(new gridpack::serial_io::SerialBusIO<DSFullNetwork>(128,
          p_network));
    p_loadIO->open(filename.c_str());
  } else {
    p_busIO->header("No Load Watch File Name Found\n");
    p_loadWatch = false;
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
    p_loadIO->close();
#else
    p_loadIO->closeChannel();
#endif
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
