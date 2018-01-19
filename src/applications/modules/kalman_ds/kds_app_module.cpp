/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   kds_app_module.cpp
 * @author Da Meng and Yousu Chen 
 * @date   1/06/2015
 *
 * @brief
 *
 * @Modyfied by Xinya Li on 8/20/2015 to include the dynamic simulation
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "kds_app_module.hpp"

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::kalman_filter::KalmanApp::KalmanApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::kalman_filter::KalmanApp::~KalmanApp(void)
{
}

/**
 * Get the time series measurements from a file
 * @param filename name of external file with time series data
 * @param series a vector of pointers to arrays of time series data
 * @param nsteps the number of timesteps in each of the time series
 * @param keys list of bus indices that own the time series data
 */
void gridpack::kalman_filter::KalmanApp::getTimeSeries(std::string filename,
    std::vector<double*> &series, double &delta_t, int &nsteps,
    std::vector<int> &keys)
{
  int me = p_comm.rank();
  int i;
  for (i=0; i<series.size(); i++) {
    delete [] series[i];
  }
  series.clear();
  keys.clear();
  if (me==0) {
    std:: ifstream input;
    input.open(filename.c_str());
    if (!input.is_open()) {
      char buf[512];
      sprintf(buf,"Failed to open time series data file: %s\n\n",
          filename.c_str());
      throw gridpack::Exception(buf);
    }
    std::string line;
    std::getline(input, line);
    std::vector<std::string> tokens;
    boost::split(tokens, line, boost::algorithm::is_any_of(","),
        boost::token_compress_on);
    int ndata = tokens.size() - 1;
    int idx = tokens[0].find_last_of('-');
    std::string token;
    int len;
    // find data after "-" mark
    if (idx != std::string::npos) {
      len = tokens[0].length()-idx;
      token = tokens[0].substr(idx+1,len);
    } else {
      char buf[512];
      sprintf(buf,"Failed to find total number of timesteps: %s\n\n",
          line.c_str());
      throw gridpack::Exception(buf);
    }
    // Get number of steps and allocate data arrays
    nsteps = atoi(token.c_str());
    double *time = new double[nsteps];
    for (i=0; i<ndata; i++) {
      double *dptr = new double[nsteps];
      series.push_back(dptr);
    }
    // Get bus IDs
    for (i=1; i<=ndata; i++) {
      idx = tokens[i].find_last_of('-');
      if (idx != std::string::npos) {
        len = tokens[i].length()-idx;
        token = tokens[i].substr(idx+1,len);
        keys.push_back(atoi(token.c_str()));
      } else {
        char buf[512];
        sprintf(buf,"Indecipherable bus ID: %s\n\n",
            tokens[i].c_str());
        throw gridpack::Exception(buf);
      }
    }
    // read remaining data
    int ncnt = 0;
    while(std::getline(input,line) && ncnt < nsteps) {
      boost::split(tokens, line, boost::algorithm::is_any_of(","),
          boost::token_compress_on);
      time[ncnt] = atof(tokens[0].c_str());
      for (i=1; i<=ndata; i++) {
        (series[i-1])[ncnt] = atof(tokens[i].c_str());
      }
      ncnt++;
    }
    if (nsteps > 1) {
      delta_t = (time[nsteps-1]-time[0])/(static_cast<double>(nsteps-1));
    }
    delete [] time;
  } else {
    delta_t = 0.0;
    nsteps = 0;
  }
  p_comm.sum(&delta_t, 1);
  p_comm.sum(&nsteps, 1);
}

/**
 * Set up time series data for all buses
 * @param network pointer to KalmanNetwork object
 * @param cursor pointer to data in input deck
 */
void gridpack::kalman_filter::KalmanApp::setTimeData(
    boost::shared_ptr<KalmanNetwork> &network,
    gridpack::utility::Configuration::CursorPtr cursor)
{
  std::string filename;
  if (!cursor->get("KalmanAngData",&filename)) {
     printf("No Kalman angle time series data specified\n");
     return;
  }
  std::vector<double*> data;
  std::vector<int> keys;
  int nsteps_ang;
  double delta_t_ang;
  gridpack::kalman_filter::KalmanApp::getTimeSeries(filename,
    data, delta_t_ang, nsteps_ang,
    keys);
  int me = network->communicator().rank();
  // Create hash distribution object
  gridpack::hash_distr::HashDistribution<KalmanNetwork,double,double>
    hash(network);
  // distribute angle time series data
  hash.distributeBusValues(keys, data, nsteps_ang);
  int i;
  gridpack::kalman_filter::KalmanBus *bus;
  int ndata = keys.size();
  for (i=0; i<ndata; i++) {
    bus = dynamic_cast<gridpack::kalman_filter::KalmanBus*>
      (network->getBus(keys[i]).get());
    bus->setVAngSeries(data[i]);
  }
  data.clear();
  keys.clear();
  if (!cursor->get("KalmanMagData",&filename)) {
     printf("No Kalman magnitude time series data specified\n");
     return;
  }
  double delta_t_mag;
  int nsteps_mag;
  gridpack::kalman_filter::KalmanApp::getTimeSeries(filename,
    data, delta_t_mag, nsteps_mag,
    keys);
  if (delta_t_ang != delta_t_mag) {
    // Some kind of error
  }
  if (nsteps_ang != nsteps_mag) {
    // Some kind of error
  }
  p_nsteps = nsteps_ang;
  p_delta_t = delta_t_ang;
  hash.distributeBusValues(keys, data, nsteps_mag);
  ndata = keys.size();
  for (i=0; i<ndata; i++) {
    bus = dynamic_cast<gridpack::kalman_filter::KalmanBus*>
      (network->getBus(keys[i]).get());
    bus->setVMagSeries(data[i]);
  }
  int nbus = network->numBuses();
  for (i=0; i<nbus; i++) {
    bus = dynamic_cast<gridpack::kalman_filter::KalmanBus*>
      (network->getBus(i).get());
    bus->setTimeSteps(delta_t_ang, nsteps_ang);
  }
//  for (i=0; i<nbus; i++) {
//    bus = dynamic_cast<gridpack::kalman_filter::KalmanBus*>
//      (network->getBus(i).get());
//    if (network->getActiveBus(i)) {
//      bus->printT();
//    }
//  }
}

/**
 * Execute application
 */
void gridpack::kalman_filter::KalmanApp::execute(int argc, char** argv)
{

  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_Total = timer->createCategory("App:Total");
  int t_PF = timer->createCategory("PF: Total");
  int t_KF = timer->createCategory("KF: Time Loop");
  int t_In = timer->createCategory("KF: Input and Initialization");
  timer->start(t_Total);
  timer->start(t_PF);
  
  gridpack::parallel::Communicator p_comm;
  // Initialize Kalman filter calculation by first running a powerflow
  // simulation
  boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network(new gridpack::powerflow::PFNetwork(p_comm));

  // Read configuration file
  gridpack::utility::Configuration *config
    = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,p_comm);
  } else {
    config->open("input.xml",p_comm);
  }

  // run powerflow calculation and save data to data collection objects
  gridpack::powerflow::PFAppModule pf_app;
  pf_app.readNetwork(pf_network,config);
  pf_app.initialize();
  pf_app.solve();
  pf_app.write();
  pf_app.saveData();
  
  // Copy data collection objects to new Kalman filter network
  // Kalman filter network
  boost::shared_ptr<gridpack::kalman_filter::KalmanNetwork>
    network(new gridpack::kalman_filter::KalmanNetwork(p_comm));
  pf_network->clone<gridpack::kalman_filter::KalmanBus,
    gridpack::kalman_filter::KalmanBranch>(network);

  timer->stop(t_PF);
  timer->start(t_In);

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<KalmanNetwork> busIO(512, network);
  gridpack::serial_io::SerialBusIO<KalmanNetwork> deltaIO(32, network);
  gridpack::serial_io::SerialBusIO<KalmanNetwork> omegaIO(32, network);
  gridpack::serial_io::SerialBranchIO<KalmanNetwork> branchIO(512, network);
  char ioBuf[128];
  
  // Read in Dynamic simulation and Kalman filter related info
  gridpack::utility::Configuration::CursorPtr cursor, secursor;
  cursor   = config->getCursor("Configuration.Dynamic_simulation");
  secursor = config->getCursor("Configuration.Kalman_filter");

  // Read in information of simulation time and timestep
  double sim_time = cursor->get("simulationTime",0.0);
  double time_step = cursor->get("timeStep",0.0);
  int KnownFault = cursor->get("KnownFault",0); 
  int TimeOffset = cursor->get("TimeOffset",0);
  int CheckEqn = cursor->get("CheckEqn",0);
  if (CheckEqn) {
    sprintf(ioBuf,"\nOnly DAE equations, no EnKF analysis!\n"); busIO.header(ioBuf);
  } else {
    sprintf(ioBuf,"\nStart EnKF analysis......\n"); busIO.header(ioBuf);
  }
  sprintf(ioBuf,"\nSimulation Time: %8.4f\n",sim_time); busIO.header(ioBuf);
  sprintf(ioBuf,"\nTime Step: %8.4f\n",time_step); busIO.header(ioBuf);
  sprintf(ioBuf,"\nTime Offset: %d\n",TimeOffset); busIO.header(ioBuf);

  // load input *.dyr file
  gridpack::parser::PTI23_parser<KalmanNetwork> parser(network);
  std::string filename = secursor->get("generatorParameters","");
  if (filename.size() > 0) {
    parser.externalParse(filename.c_str());
  } else {
    printf("No generator parameters found for file (%s)\n",
        filename.c_str());
    return;
  }

  int idx;
  // Read in information about fault events and store them in internal data
  // structure
  cursor = config->getCursor("Configuration.Dynamic_simulation.faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  if (cursor) cursor->children(events);
  std::vector<gridpack::kalman_filter::KalmanBranch::Event>
     faults = setFaultEvents(events);
  // Index that if the fault event is known
  if (KnownFault) {
  sprintf(ioBuf,"\nFault event is known\n"); busIO.header(ioBuf);
    for (idx=0; idx<events.size(); idx++) {
      sprintf(ioBuf,"\nFault Events No.: %d\n",idx+1); busIO.header(ioBuf);
      sprintf(ioBuf,"Begin Fault: %8.4f\n",faults[idx].start); busIO.header(ioBuf);
      sprintf(ioBuf,"End Fault: %8.4f\n",faults[idx].end); busIO.header(ioBuf);
      sprintf(ioBuf,"From Branch: %d  To Branch: %d\n",faults[idx].from_idx,faults[idx].to_idx); busIO.header(ioBuf);
      sprintf(ioBuf,"Time Step: %8.4f\n",faults[idx].step); busIO.header(ioBuf);
    }
  } else {
  sprintf(ioBuf,"\nFault event is unknown\n"); busIO.header(ioBuf);
  }
  
  // Convergence and iteration parameters
  double tolerance = secursor->get("tolerance",1.0e-6);
  int max_iteration = secursor->get("maxIteration",20);
  sprintf(ioBuf,"\nTolerance: %16.8f\n",tolerance);
  busIO.header(ioBuf);
  sprintf(ioBuf,"Maximum Iterations: %d\n",max_iteration);
  busIO.header(ioBuf);

  // Ensemble parameters
  int nsize = secursor->get("ensembleSize",20);
  double sigma = secursor->get("gaussianWidth",0.1);
  double noise = secursor->get("noiseScale",0.1);
  int iseed = secursor->get("randomSeed",11238);
  int maxstep = secursor->get("maxSteps",0);
  double Rm1 = 1.0/(noise*noise);
  double Rm1n, N_inv;
  if (CheckEqn) {
    nsize = 1; sigma = 0.0; noise = 0.0;
  }
  if (nsize > 1) {
    Rm1n= Rm1/static_cast<double>(nsize-1);
    N_inv = 1.0/static_cast<double>(nsize-1);
  } else {
    Rm1n = 0.0;
    N_inv = 0.0;
  }
  sprintf(ioBuf,"Ensemble Size: %d\n",nsize); busIO.header(ioBuf);
  sprintf(ioBuf,"Gaussian Width: %16.8f\n",sigma); busIO.header(ioBuf);
  sprintf(ioBuf,"Noise Scale: %16.8f\n",noise); busIO.header(ioBuf);
  sprintf(ioBuf,"Random Number Seed: %d\n",iseed); busIO.header(ioBuf);
  
  // Initialize random number generator
  gridpack::random::Random random;
  random.seed(iseed);
 
  // create factory
  gridpack::kalman_filter::KalmanFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  factory.setExchange();

  printf("p[%d] Got to 1 nsize: %d\n",p_comm.rank(),nsize);
  // Set ensemble parameters
  factory.setEnsembleSize(nsize);
  factory.setGaussianWidth(sigma,noise);
  printf("p[%d] Got to 2\n",p_comm.rank());


  // Get time series data
  setTimeData(network, secursor);
  sprintf(ioBuf,"Number of Time Steps: %d\n",p_nsteps);
  busIO.header(ioBuf);
  if (maxstep > 0 && maxstep < p_nsteps) {
    p_nsteps = maxstep;
    sprintf(ioBuf,"Maximum Steps: %d\n",p_nsteps);
    busIO.header(ioBuf);
  }

  // Create bus data exchange
  network->initBusUpdate();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  printf("p[%d] Got to 3\n",p_comm.rank());
  // create ensembles
  factory.createEnsemble();
  printf("p[%d] Got to 4\n",p_comm.rank());

  // Create the Kalman ensemble X matrix
  factory.setMode(EnsembleX);
  gridpack::mapper::GenSlabMap<KalmanNetwork> xSlab(network);
  printf("p[%d] Got to 5\n",p_comm.rank());
  boost::shared_ptr<gridpack::math::Matrix> X = xSlab.mapToMatrix();
  printf("p[%d] Got to 6\n",p_comm.rank());
  //X -> print();
  //X->save("X.m");

  // Create initial Y matrix
  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<KalmanNetwork> ybusMap(network);
  factory.setMode(RefreshY);
  boost::shared_ptr<gridpack::math::Matrix> Y_init = ybusMap.mapToMatrix();
  //Y_init -> print();
  //Y_init->save("Y_init.m");
  // Use linear solution to get inverse of Y matrix
  gridpack::math::LinearMatrixSolver solver1(*Y_init);
  secursor = config->getCursor("Configuration.Kalman_filter");
  solver1.configure(secursor);
  // Create Y_c matrix
  factory.setMode(YC);
  gridpack::mapper::FullMatrixMap<KalmanNetwork> ycMap(network);

  // Create dense version of Y_c
  boost::shared_ptr<gridpack::math::Matrix> Y_c = ycMap.mapToMatrix(true);
  //Y_c -> print();
  //Y_c->save("Y_c.m");
  // Evaluate RecV_0(pre-fault) matrix by solving Y_init*RecV_0 = Y_c.

  boost::shared_ptr<gridpack::math::Matrix> RecV_0(solver1.solve(*Y_c));
//  boost::shared_ptr<gridpack::math::Matrix> RecV(solver1.solve(*Y_c));
  gridpack::math::Matrix *RecV;
  //RecV_0 -> print();
  //RecV_0->save("RecV_0.m");
  // Get RecV_1(fault)
  int sw2_2 = faults[0].from_idx - 1;
  int sw3_2 = faults[0].to_idx - 1;
  boost::shared_ptr<gridpack::math::Matrix> Yd_1(Y_init->clone());
  //gridpack::ComplexType x(0.0, -1e7);
  factory.setEvent(faults[0]);
  factory.setMode(onFY);
  ybusMap.overwriteMatrix(Yd_1); 
  gridpack::math::LinearMatrixSolver solver2(*Yd_1);
  boost::shared_ptr<gridpack::math::Matrix> RecV_1(solver2.solve(*Y_c));
  //RecV_1->save("RecV_1.m");
/*  //Get RecV_2(post-fault)
  boost::shared_ptr<gridpack::math::Matrix> Yd_2(Y_init->clone());
  factory.setMode(posFY);
  ybusMap.incrementMatrix(Yd_2);
//  Yd_2->save("Yd_2.m");
  gridpack::math::LinearMatrixSolver solver3(*Yd_2);
  boost::shared_ptr<gridpack::math::Matrix> RecV_2(solver3.solve(*Y_c));
//  RecV_2->print();
//  RecV_2->save("RecV_2.m"); */

  // Set up io routines
  deltaIO.open("delta.dat");
  omegaIO.open("omega.dat");

  // Create eSlab mapper
  factory.setMode(E_Ensemble1);
  gridpack::mapper::GenSlabMap<KalmanNetwork> eSlab(network);

  // Create hxSlab mapper
  factory.setMode(HX);
  gridpack::mapper::GenSlabMap<KalmanNetwork> hxSlab(network);

  // Create vSlab mapper
  factory.setMode(V1);
  gridpack::mapper::GenSlabMap<KalmanNetwork> vSlab(network);

  // Create v3Slab mapper
  factory.setMode(V3);
  gridpack::mapper::GenSlabMap<KalmanNetwork> v3Slab(network);

  // Create measurement matrix for timestep = 0. You can use mapper for the
  // ensemble HX matrix
  factory.setMode(Measurements);
  factory.setCurrentTimeStep(TimeOffset+1);
  boost::shared_ptr<gridpack::math::Matrix> D = hxSlab.mapToMatrix();

  sprintf(ioBuf,"%12.6f",static_cast<double>(0.0));
  deltaIO.header(ioBuf);
  deltaIO.write("delta");
  deltaIO.header("\n");
  omegaIO.header(ioBuf);
  omegaIO.write("omega");
  omegaIO.header("\n");

  timer->stop(t_In);
  timer->start(t_KF);

  //-----------------------------------------------------------------------
  // Ensemble Kalman filters 
  //-----------------------------------------------------------------------
  // Time-step related variables
  int simu_k;
  int t_step[20];
  int flagF;
  int steps3, steps2, steps1;
  int I_Steps;

  const double sysFreq = 60.0;
  double pi = 4.0*atan(1.0);
  const double basrad = 2.0 * pi * sysFreq;
  gridpack::ComplexType jay(0.0, 1.0);

  // Set up switch info
  int nswtch = 4;
  static double sw1[4];
  static double sw7[4];
  sw1[0] = 0.0;
  sw1[1] = faults[0].start;
  sw1[2] = faults[0].end;
  sw1[3] = sim_time;
  sw7[0] = time_step;
  sw7[1] = faults[0].step;
  sw7[2] = time_step;
  sw7[3] = time_step;
  simu_k = 0;
  for (int i = 0; i < nswtch-1; i++) {
    t_step[i] = (int) (round((sw1[i+1] -sw1[i]) / sw7[i]));
    simu_k += t_step[i];
  }
  simu_k++;
  if (simu_k > p_nsteps) simu_k = p_nsteps;

  steps3 = t_step[0] + t_step[1] + t_step[2];
  steps2 = t_step[0] + t_step[1];
  steps1 = t_step[0];
  
  sprintf(ioBuf,"Start Time Step of Fault: %d\n",steps1); busIO.header(ioBuf);
  sprintf(ioBuf,"End Time Step of Fault: %d\n",steps2); busIO.header(ioBuf);
   

  for (I_Steps = 2; I_Steps < simu_k; I_Steps++) { // Simulation Steps

    int t_onlyDAE= timer->createCategory("KF: In-Loop Only DAE");
    int t_EnKF = timer->createCategory("KF: In-Loop EnKF");
    int t_Output = timer->createCategory("KF: In-Loop Output");
     
    timer->start(t_onlyDAE);

    int t_selectRecV = timer->createCategory("KF: In-Loop Select RecV");
    timer->start(t_selectRecV);    

    if ((I_Steps <= steps1) || (I_Steps > steps2 + 1)) {
      flagF = 0; // pre and post fault
    } else if ((I_Steps > steps1 + 1) && (I_Steps <= steps2)) {
      flagF = 1; // fault
    } else if (I_Steps == steps1 + 1) {
      flagF = 2; // first switch
    } else if (I_Steps == steps2 + 1) {
      flagF = 3; // second switch
    } else {
      flagF = 0;
    }

    if (KnownFault) {
      if (flagF == 0) {
        RecV = RecV_0.get(); 
      } else if (flagF == 1) {
        RecV = RecV_1.get(); 
      } else if (flagF == 2) {
        RecV = RecV_0.get(); 
      } else if (flagF == 3) {
        RecV = RecV_1.get();
      }
    } else {
      RecV = RecV_0.get();
    }

    timer->stop(t_selectRecV);
    
    // Create E_ensemble 1 matrix
    factory.setMode(E_Ensemble1);
    boost::shared_ptr<gridpack::math::Matrix> E_ensmb1 = eSlab.mapToMatrix();
    //E_ensmb1 -> print();    

    // Create V1
    boost::shared_ptr<gridpack::math::Matrix> v1(multiply(*RecV, *E_ensmb1));
    //v1 -> print();
    
    // Push elements of V1 back onto buses
    factory.setMode(V1);
    vSlab.mapToNetwork(v1);

    // Create elements of X2
    factory.evaluateX2();

    if (KnownFault) {
      if (flagF == 0) {
        RecV = RecV_0.get();
      } else if (flagF == 1) {
        RecV = RecV_1.get();
      } else if (flagF == 2) {
        RecV = RecV_1.get();
      } else if (flagF == 3) {
        RecV = RecV_0.get();
      }
    } else {
      RecV = RecV_0.get();
    }

    // Create E_ensemble 2 matrix
    factory.setMode(E_Ensemble2);
    boost::shared_ptr<gridpack::math::Matrix> E_ensmb2 = eSlab.mapToMatrix();

    // Create V2
    boost::shared_ptr<gridpack::math::Matrix> v2(multiply(*RecV, *E_ensmb2));

    // Push elements of V2 back onto buses
    factory.setMode(V2);
    vSlab.mapToNetwork(v2);

    // Create elements of X3
    factory.evaluateX3();
    
    timer->stop(t_onlyDAE);
    timer->start(t_EnKF);
    
if (!(CheckEqn)) {
    int t_A = timer->createCategory("KF: In-Loop EnKF A");
    timer->start(t_A);    
    // Create perturbation matrix for X3
    factory.setMode(Perturbation);
    boost::shared_ptr<gridpack::math::Matrix> A = xSlab.mapToMatrix();
    timer->stop(t_A);

    int t_ensmb3 = timer->createCategory("KF: In-Loop EnKF E_ensmb3");
    timer->start(t_ensmb3);
    // Create E_ensemble 3 matrix
    factory.setMode(E_Ensemble3);
    boost::shared_ptr<gridpack::math::Matrix> E_ensmb3 = eSlab.mapToMatrix();
    timer->stop(t_ensmb3);

    int t_V3 = timer->createCategory("KF: In-Loop EnKF V3");
    timer->start(t_V3);
    int t_V3m = timer->createCategory("KF: In-Loop EnKF RecV*Ensmb3");
    timer->start(t_V3m);
    // Create V3
    boost::shared_ptr<gridpack::math::Matrix> v3(multiply(*RecV, *E_ensmb3));
    timer->stop(t_V3m);

    // Push elements of V3 back onto buses
    factory.setMode(V3);
    v3Slab.mapToNetwork(v3);
    timer->stop(t_V3);

    // Create HX matrix
    int t_HX = timer->createCategory("KF: In-Loop EnKF HX");
    timer->start(t_HX);
    factory.setMode(HX);
    boost::shared_ptr<gridpack::math::Matrix> HX = hxSlab.mapToMatrix();
    timer->stop(t_HX);

    // Create Y = D-HX
    int t_Y = timer->createCategory("KF: In-Loop EnKF Y");
    timer->start(t_Y);
    boost::shared_ptr<gridpack::math::Matrix> Y(D->clone());
    HX->scale(-1.0);
    Y->add(*HX);
    HX->scale(-1.0);
    timer->stop(t_Y);
    
    int t_Q = timer->createCategory("KF: In-Loop EnKF Q");
    timer->start(t_Q);
    // Create HA matrix
    int t_setHA = timer->createCategory("KF: In-Loop EnKF setHA");
    timer->start(t_setHA);
    factory.setMode(HA);
    timer->stop(t_setHA);
    int t_HA = timer->createCategory("KF: In-Loop EnKF HA");
    timer->start(t_HA);
    boost::shared_ptr<gridpack::math::Matrix> HA = hxSlab.mapToMatrix();
    timer->stop(t_HA);
    int t_HAt = timer->createCategory("KF: In-Loop EnKF HAt");
    timer->start(t_HAt);
    boost::shared_ptr<gridpack::math::Matrix> HA_t(transpose(*HA));
    timer->stop(t_HAt);
    int t_HAtHA = timer->createCategory("KF: In-Loop EnKF HAtHA");
    timer->start(t_HAtHA);
    boost::shared_ptr<gridpack::math::Matrix> Q(multiply(*HA_t,*HA));
    timer->stop(t_HAtHA);
    int t_ScaleQ = timer->createCategory("KF: In-Loop EnKF ScaleQ");
    timer->start(t_ScaleQ);
    Q->scale(Rm1n);
    boost::shared_ptr<gridpack::math::Matrix> H1(Q->clone());
    gridpack::ComplexType z_one(1.0,0.0);
    Q->addDiagonal(z_one);
    timer->stop(t_ScaleQ);
    timer->stop(t_Q);

    int t_Z1 = timer->createCategory("KF: In-Loop EnKF Z1");
    timer->start(t_Z1);
    // Create Z1 matrix
    boost::shared_ptr<gridpack::math::Matrix> Z1(multiply(*HA_t,*Y));
    Z1->scale(Rm1);
    timer->stop(t_Z1);

    int t_W = timer->createCategory("KF: In-Loop EnKF Solve W"); 
    timer->start(t_W);
    // Create W by solving Q*W = Z1
    boost::shared_ptr<gridpack::math::Matrix>
      Q_sparse(gridpack::math::storageType(*Q,
            gridpack::math::Sparse));
    gridpack::math::LinearMatrixSolver solver2(*Q_sparse);
    secursor = config->getCursor("Configuration.Kalman_filter");
    solver2.configure(secursor);
    boost::shared_ptr<gridpack::math::Matrix> W(solver2.solve(*Z1));
    timer->stop(t_W);

    int t_Z2 = timer->createCategory("KF: In-Loop EnKF Z2");
    timer->start(t_Z2);
    // Evaluate Z2 = Z1 - H1*W
    boost::shared_ptr<gridpack::math::Matrix> Z2(multiply(*H1,*W));
    Z2->scale(-1);
    Z2->add(*Z1);
    timer->stop(t_Z2);
    
    int t_Update = timer->createCategory("KF: In-Loop EnKF X Update");
    timer->start(t_Update);
    int t_X_inc = timer->createCategory("KF: In-Loop EnKF X_inc");
    timer->start(t_X_inc);
    // Evaluate X_inc
    boost::shared_ptr<gridpack::math::Matrix> X_inc(multiply(*A,*Z2));
    X_inc->scale(N_inv);
    timer->stop(t_X_inc);

    // Push results back onto buses and update values of rotor angle and speed
    factory.setMode(X_INC);
    xSlab.mapToNetwork(X_inc);
    timer->stop(t_Update);
} else {
    factory.setMode(X_Update);
    xSlab.mapToNetwork(X);
}
    factory.setCurrentTimeStep(TimeOffset+I_Steps);
    timer->stop(t_EnKF);
    timer->start(t_Output);

    char buf[128];
    sprintf(buf,"\n\n  Results for Timestep %d\n",I_Steps);
    busIO.header(buf);
    sprintf(buf,"\n    Bus ID     Gen. ID         Delta             Omega\n\n");
    busIO.header(buf);
    busIO.write("xnew");
    
    // Create measurement matrix for next timestep
    factory.setMode(Measurements);
    hxSlab.mapToMatrix(D);

    sprintf(ioBuf,"%12.6f",static_cast<double>(I_Steps-1)*p_delta_t);
    deltaIO.header(ioBuf);
    deltaIO.write("delta");
    deltaIO.header("\n");
    omegaIO.header(ioBuf);
    omegaIO.write("omega");
    omegaIO.header("\n");
    timer->stop(t_Output);
  }
  deltaIO.close();
  omegaIO.close();
  sprintf(ioBuf,"\nEnd EnKF analysis......\n"); busIO.header(ioBuf);
  
  timer->stop(t_KF);
  timer->stop(t_Total);
  timer->dump();
}


/**
 * Utility function to convert faults that are in event list into
 * internal data structure that can be used by code
 * @param cursors list of cursors pointing to individual events in input
 * deck
 * @return list of event data structures
 */
std::vector<gridpack::kalman_filter::KalmanBranch::Event>
   gridpack::kalman_filter::KalmanApp::setFaultEvents(
   std::vector<gridpack::utility::Configuration::CursorPtr > events)
{
  int size = events.size();
  int idx;
  std::vector<gridpack::kalman_filter::KalmanBranch::Event> faults;
  // Parse fault events
  for (idx=0; idx<size; idx++) {
    gridpack::kalman_filter::KalmanBranch::Event event;
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
    } else {
      event.from_idx = 0;
      event.to_idx = 0;
    }
    event.step = events[idx]->get("timeStep",0.0);
    if (event.step != 0.0 && event.end != 0.0 && event.from_idx != event.to_idx) {
      faults.push_back(event);
    }
  }
  return faults;
}
