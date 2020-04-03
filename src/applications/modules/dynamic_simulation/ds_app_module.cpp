/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shuangshuang Jin
 * @date   2015-03-06 14:51:59 d3g096
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "ds_app_module.hpp"

// Calling program for dynamic simulation application

/**
 * Basic constructor
 */
gridpack::dynamic_simulation_r::DSAppModule::DSAppModule(void)
{
  p_generatorWatch = false;
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation_r::DSAppModule::~DSAppModule(void)
{
}

/**
 * Read in and partition the powerflow network. The input file is read
 * directly from the Dynamic_simulation block in the configuration file so no
 * external file names or parameters need to be passed to this routine
 * @param network pointer to a DSNetwork object. This should not have any
 * buses or branches defined on it.
 * @param config pointer to open
 * configuration file
 */
void gridpack::dynamic_simulation_r::DSAppModule::readNetwork(
    boost::shared_ptr<DSNetwork> &network,
    gridpack::utility::Configuration *config)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);
  p_comm = network->communicator();
  p_network = network;
  p_config = config;

  int t_setup = timer->createCategory("Dynamic Simulation: Read in Data");
  timer->start(t_setup);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  std::string filename;
  if (!cursor->get("networkConfiguration",&filename)) {
      printf("No network configuration specified\n");
  }
  p_sim_time = cursor->get("simulationTime",0.0);
  if (p_sim_time == 0.0) {
    // TODO: some kind of error
  }
  p_time_step = cursor->get("timeStep",0.0);
  if (p_time_step == 0.0) {
    // TODO: some kind of error
  }

  // load input file
  gridpack::parser::PTI23_parser<DSNetwork> parser(network);
  parser.parse(filename.c_str());
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  filename = cursor->get("generatorParameters","");
  if (filename.size() > 0) parser.parse(filename.c_str());
  timer->stop(t_setup);

  int t_part = timer->createCategory("Dynamic Simulation: Partition Network");
  timer->start(t_part);
  // partition network
  network->partition();
  timer->stop(t_part);

  // Create serial IO object to export data from buses or branches
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<DSNetwork>(2048, network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<DSNetwork>(128, network));
  timer->stop(t_total);
}

/**
 * Read faults from external file and form a list of faults
 * @param cursor pointer to open file contain fault or faults
 * @return a list of fault events
 */
std::vector<gridpack::dynamic_simulation_r::DSBranch::Event>
  gridpack::dynamic_simulation_r::DSAppModule::
  getFaults(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  std::vector<gridpack::dynamic_simulation_r::DSBranch::Event> ret;
  if (list) {
    list->children(events);
    int size = events.size();
    int idx;
    // Parse fault events
    for (idx=0; idx<size; idx++) {
      gridpack::dynamic_simulation_r::DSBranch::Event event;
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
        ret.push_back(event);
      }
    }
  }
  return ret;
}

/**
 * Read generator parameters. These will come from a separate file (most
 * likely). The name of this file comes from the input configuration file.
 */
void gridpack::dynamic_simulation_r::DSAppModule::readGenerators()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  gridpack::parser::PTI23_parser<DSNetwork> parser(p_network);
  std::string filename = cursor->get("generatorParameters","");
  if (filename.size() > 0) parser.externalParse(filename.c_str());
}

void gridpack::dynamic_simulation_r::DSAppModule::setNetwork(
    boost::shared_ptr<DSNetwork> &network,
    gridpack::utility::Configuration *config)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);
  p_comm = network->communicator();
  p_network = network;
  p_config = config;

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
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
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<DSNetwork>(2048, network));
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<DSNetwork>(128, network));
  timer->stop(t_total);
}

/**
 * Set up exchange buffers and other internal parameters and
 * initialize
 * network components using data from data collection
 */
void gridpack::dynamic_simulation_r::DSAppModule::initialize()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);
  int t_config = timer->createCategory("Dynamic Simulation: Configure Network");
  timer->start(t_config);
  // create factory
  p_factory.reset(new gridpack::dynamic_simulation_r::DSFactoryModule(p_network));
  p_factory->load();

  // set network components using factory
  p_factory->setComponents();
  timer->stop(t_config);

  // set YBus components so that you can create Y matrix  
  p_factory->setYBus();

  if (!p_factory->checkGen()) {
    p_busIO->header("Missing generators on at least one processor\n");
    timer->stop(t_total);
    return;
  }
  timer->stop(t_total);
}

void gridpack::dynamic_simulation_r::DSAppModule::solve(
    gridpack::dynamic_simulation_r::DSBranch::Event fault)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
#if 0
  timer->start(t_total);
  int t_matset = timer->createCategory("Dynamic Simulation: Matrix Setup");
  timer->start(t_matset);
  p_factory->setMode(YBUS);
  gridpack::mapper::FullMatrixMap<DSNetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  ///p_branchIO.header("\n=== orginal ybus: ============\n");
  ///orgYbus->print();

  // Form constant impedance load admittance yl for all buses and add it to 
  // system Y matrix: ybus = ybus + yl
  p_factory->setMode(YL);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  ///p_branchIO.header("\n=== ybus after added yl: ============\n");
  ///ybus->print();
 
  // Construct permutation matrix perm by checking for multiple generators at a bus
  p_factory->setMode(PERM);
  gridpack::mapper::FullMatrixMap<DSNetwork> permMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> perm = permMap.mapToMatrix();
  timer->stop(t_matset);
  ///p_busIO->header("\n=== perm: ============\n");
  ///perm->print(); 

  // Form a transposed matrix of perm
  int t_trans = timer->createCategory("Dynamic Simulation: Matrix Transpose");
  timer->start(t_trans);
   boost::shared_ptr<gridpack::math::Matrix> permTrans(transpose(*perm));
  timer->stop(t_trans);
  ///p_busIO->header("\n=== permTrans: ============\n");
  ///permTrans->print();

  // Construct matrix Y_a using extracted xd and ra from gen data, 
  // and construct its diagonal matrix diagY_a
  timer->start(t_matset);
  p_factory->setMode(YA);
  gridpack::mapper::FullMatrixMap<DSNetwork> yaMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> diagY_a = yaMap.mapToMatrix();
  ///p_busIO->header("\n=== diagY_a: ============\n");
  ///diagY_a->print(); 
  // Convert diagY_a from sparse matrix to dense matrix Y_a so that SuperLU_DIST can solve
  gridpack::math::MatrixStorageType denseType = gridpack::math::Dense;
  boost::shared_ptr<gridpack::math::Matrix> Y_a(gridpack::math::storageType(*diagY_a, denseType));
  timer->stop(t_matset);

  // Construct matrix Ymod: Ymod = diagY_a * permTrans
  int t_matmul = timer->createCategory("Dynamic Simulation: Matrix Multiply");
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> Ymod(multiply(*diagY_a, *permTrans));
  timer->stop(t_matmul);
  ///p_busIO->header("\n=== Ymod: ============\n");
  ///Ymod->print(); 
 
  // Form matrix Y_b: Y_b(1:ngen, jg) = -Ymod, where jg represents the 
  // corresponding index sets of buses that the generators are connected to. 
  // Then construct Y_b's transposed matrix Y_c: Y_c = Y_b'
  // This two steps can be done by using a P matrix to get Y_c directly.
  timer->start(t_matset);
  p_factory->setMode(PMatrix);
  gridpack::mapper::FullMatrixMap<DSNetwork> pMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> P = pMap.mapToMatrix();
  ///p_busIO->header("\n=== P: ============\n");
  ///P->print();
  p_factory->setMode(YC);
  gridpack::mapper::FullMatrixMap<DSNetwork> cMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Y_c = cMap.mapToMatrix();
  Y_c->scale(-1.0);
  ///p_busIO->header("\n=== Y_c: ============\n");
  ///Y_c->print();
  p_factory->setMode(YB);
  gridpack::mapper::FullMatrixMap<DSNetwork> bMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Y_b = bMap.mapToMatrix();
  timer->stop(t_matset);
  ///p_busIO->header("\n=== Y_b: ============\n");
  ///Y_b->print();
  
  Y_c->scale(-1.0); // scale Y_c by -1.0 for linear solving
  // Convert Y_c from sparse matrix to dense matrix Y_cDense so that SuperLU_DIST can solve
  //gridpack::math::Matrix::StorageType denseType = gridpack::math::Matrix::Dense;
  timer->start(t_matset);
  boost::shared_ptr<gridpack::math::Matrix> Y_cDense(gridpack::math::storageType(*Y_c, denseType));
   
  // Form matrix permYmod
  p_factory->setMode(permYMOD);
  gridpack::mapper::FullMatrixMap<DSNetwork> pymMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix>  permYmod= pymMap.mapToMatrix();
  ///p_busIO->header("\n=== permYmod: ============\n");
  ///permYmod->print();

  // Update ybus: ybus = ybus+permYmod (different dimension) => prefy11ybus
  p_factory->setMode(updateYbus);

  boost::shared_ptr<gridpack::math::Vector> vPermYmod(diagonal(*permYmod));
  ///p_busIO->header("\n=== vPermYmod: ============\n");
  ///vPermYmod->print();
  gridpack::mapper::BusVectorMap<DSNetwork> permYmodMap(p_network);
  permYmodMap.mapToBus(vPermYmod);

  boost::shared_ptr<gridpack::math::Matrix> prefy11ybus = ybusMap.mapToMatrix();
  timer->stop(t_matset);
  ///p_branchIO.header("\n=== prefy11ybus: ============\n");
  ///prefy11ybus->print();

  // Solve linear equations of ybus * X = Y_c
  //gridpack::math::LinearSolver solver1(*prefy11ybus);
  //solver1.configure(cursor);
  int t_solve = timer->createCategory("Dynamic Simulation: Solve Linear Equation");
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver1(*prefy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> prefy11X(solver1.solve(*Y_cDense));
  timer->stop(t_solve);
  ///p_branchIO.header("\n=== prefy11X: ============\n");
  ///prefy11X->print(); 
  
  //-----------------------------------------------------------------------
  // Compute prefy11
  //-----------------------------------------------------------------------
  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> prefy11(multiply(*Y_b, *prefy11X)); 
  timer->stop(t_matmul);
  // Update prefy11: prefy11 = Y_a + prefy11
  prefy11->add(*Y_a);
  ///p_branchIO.header("\n=== Reduced Ybus: prefy11: ============\n");
  ///prefy11->print();

  //-----------------------------------------------------------------------
  // Compute fy11
  // Update ybus values at fault stage
  //-----------------------------------------------------------------------
  // Read the switch info from faults Event from input.xml
  int sw2_2 = fault.from_idx - 1;
  int sw3_2 = fault.to_idx - 1;

  boost::shared_ptr<gridpack::math::Matrix> fy11ybus(prefy11ybus->clone());
  ///p_branchIO.header("\n=== fy11ybus(original): ============\n");
  ///fy11ybus->print();
  gridpack::ComplexType x(0.0, -1e7);
#if 0
  fy11ybus->setElement(sw2_2, sw2_2, -x);
  fy11ybus->ready();
#else
  timer->start(t_matset);
  p_factory->setEvent(fault);
  p_factory->setMode(onFY);
  ybusMap.overwriteMatrix(fy11ybus);
  timer->stop(t_matset);
#endif
  ///p_branchIO.header("\n=== fy11ybus: ============\n");
  ///fy11ybus->print();

  // Solve linear equations of fy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver2(*fy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> fy11X(solver2.solve(*Y_cDense)); 
  timer->stop(t_solve);
  ///p_branchIO.header("\n=== fy11X: ============\n");
  ///fy11X->print();
  
  // Form reduced admittance matrix fy11: fy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> fy11(multiply(*Y_b, *fy11X)); 
  timer->stop(t_matmul);
  // Update fy11: fy11 = Y_a + fy11
  fy11->add(*Y_a);
  ///p_branchIO.header("\n=== Reduced Ybus: fy11: ============\n");
  ///fy11->print();

  //-----------------------------------------------------------------------
  // Compute posfy11
  // Update ybus values at clear fault stage
  //-----------------------------------------------------------------------
  // Get the updating factor for posfy11 stage ybus
#if 0
  gridpack::ComplexType myValue = p_factory->setFactor(sw2_2, sw3_2);
#endif
  timer->start(t_matset);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus(prefy11ybus->clone());
  timer->stop(t_matset);
  ///p_branchIO.header("\n=== posfy11ybus (original): ============\n");
  ///posfy11ybus->print();
#if 0
  gridpack::ComplexType big(0.0, 1e7);
  gridpack::ComplexType x11 = big - myValue;
  posfy11ybus->addElement(sw2_2, sw2_2, x+x11);

  gridpack::ComplexType x12 = myValue;
  posfy11ybus->addElement(sw2_2, sw3_2, x12);

  gridpack::ComplexType x21 = myValue;
  posfy11ybus->addElement(sw3_2, sw2_2, x21);

  gridpack::ComplexType x22 = -myValue; 
  posfy11ybus->addElement(sw3_2, sw3_2, x22);
  
  posfy11ybus->ready(); 
#else
  timer->start(t_matset);
  p_factory->setMode(posFY);
  ybusMap.incrementMatrix(posfy11ybus);
  timer->stop(t_matset);
#endif
  ///p_branchIO.header("\n=== posfy11ybus: ============\n");
  ///posfy11ybus->print();
    
  // Solve linear equations of posfy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver3(*posfy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> posfy11X(solver3.solve(*Y_cDense)); 
  timer->stop(t_solve);
  ///p_branchIO.header("\n=== posfy11X: ============\n");
  ///posfy11X->print();
  
  // Form reduced admittance matrix posfy11: posfy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> posfy11(multiply(*Y_b, *posfy11X)); 
  timer->stop(t_matmul);
  // Update posfy11: posfy11 = Y_a + posfy11
  posfy11->add(*Y_a);
  ///p_branchIO.header("\n=== Reduced Ybus: posfy11: ============\n");
  ///posfy11->print();

  //-----------------------------------------------------------------------
  // Integration implementation (Modified Euler Method)
  //-----------------------------------------------------------------------
 
  // Map to create vector pelect  
  int t_vecset = timer->createCategory("Dynamic Simulation: Setup Vector");
  timer->start(t_vecset);
  p_factory->setMode(init_pelect);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap1(p_network);
  boost::shared_ptr<gridpack::math::Vector> pelect = XMap1.mapToVector();
  ///p_busIO->header("\n=== pelect: ===\n");
  ///pelect->print();

  // Map to create vector eprime_s0 
  p_factory->setMode(init_eprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap2(p_network);
  boost::shared_ptr<gridpack::math::Vector> eprime_s0 = XMap2.mapToVector();
  ///p_busIO->header("\n=== eprime: ===\n");
  ///eprime_s0->print();

  // Map to create vector mac_ang_s0
  p_factory->setMode(init_mac_ang);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap3(p_network);
  boost::shared_ptr<gridpack::math::Vector> mac_ang_s0 = XMap3.mapToVector();
  ///p_busIO->header("\n=== mac_ang_s0: ===\n");
  ///mac_ang_s0->print();

  // Map to create vector mac_spd_s0
  p_factory->setMode(init_mac_spd);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap4(p_network);
  boost::shared_ptr<gridpack::math::Vector> mac_spd_s0 = XMap4.mapToVector();
  ///p_busIO->header("\n=== mac_spd_s0: ===\n");
  ///mac_spd_s0->print();

  // Map to create vector eqprime
  p_factory->setMode(init_eqprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap5(p_network);
  boost::shared_ptr<gridpack::math::Vector> eqprime = XMap5.mapToVector();
  ///p_busIO->header("\n=== eqprime: ===\n");
  ///eqprime->print();

  // Map to create vector pmech  
  p_factory->setMode(init_pmech);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap6(p_network);
  boost::shared_ptr<gridpack::math::Vector> pmech = XMap6.mapToVector();
  ///p_busIO->header("\n=== pmech: ===\n");
  ///pmech->print();

  // Map to create vector mva
  p_factory->setMode(init_mva);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap7(p_network);
  boost::shared_ptr<gridpack::math::Vector> mva = XMap7.mapToVector();
  ///p_busIO->header("\n=== mva: ===\n");
  ///mva->print();

  // Map to create vector d0
  p_factory->setMode(init_d0);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap8(p_network);
  boost::shared_ptr<gridpack::math::Vector> d0 = XMap8.mapToVector();
  ///p_busIO->header("\n=== d0: ===\n");
  ///d0->print();

  // Map to create vector h
  p_factory->setMode(init_h);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap9(p_network);
  boost::shared_ptr<gridpack::math::Vector> h = XMap9.mapToVector();
  ///p_busIO->header("\n=== h: ===\n");
  ///h->print();

  ///p_busIO->header("\n============Start of Simulation:=====================\n");

  // Declare vector mac_ang_s1, mac_spd_s1
  boost::shared_ptr<gridpack::math::Vector> mac_ang_s1(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> mac_spd_s1(mac_spd_s0->clone()); 
 
  // Declare vector eprime_s1
  boost::shared_ptr<gridpack::math::Vector> eprime_s1(mac_ang_s0->clone());
 
  // Declare vector dmac_ang_s0, dmac_spd_s0, dmac_ang_s1, dmac_spd_s1
  boost::shared_ptr<gridpack::math::Vector> dmac_ang_s0(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_ang_s1(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_spd_s0(mac_spd_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_spd_s1(mac_spd_s0->clone()); 

  // Declare vector curr
  boost::shared_ptr<gridpack::math::Vector> curr(mac_ang_s0->clone());
  timer->stop(t_vecset);

  timer->start(t_trans);
  boost::shared_ptr<gridpack::math::Matrix> trans_prefy11(transpose(*prefy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_fy11(transpose(*fy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_posfy11(transpose(*posfy11));
  timer->stop(t_trans);

  timer->start(t_vecset);
  boost::shared_ptr<gridpack::math::Vector> vecTemp(mac_ang_s0->clone());
  timer->stop(t_vecset);

  // Simulation related variables
  int simu_k;
  int t_step[20];
  double t_width[20];
  int flagF;
  int S_Steps;
  int last_S_Steps;
  int steps3, steps2, steps1;
  double h_sol1, h_sol2;
  int flagF1, flagF2;
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

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagF1 = 0; 
  flagF2 = 0; 
  p_S_Steps = 1; 
  last_S_Steps = -1;

  if (p_generatorWatch) {
    p_generatorIO->header("t");
    p_generatorIO->write("watch_header");
    p_generatorIO->header("\n");
  }
  for (I_Steps = 0; I_Steps < simu_k+1; I_Steps++) {
    if (I_Steps < steps1) {
      p_S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 0;
    } else if (I_Steps == steps1) { 
      p_S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 1;
    } else if (I_Steps == steps1+1) {
      p_S_Steps = I_Steps;
      flagF1 = 1;
      flagF2 = 1;
    } else if ((I_Steps>steps1+1) && (I_Steps<steps2+1)) {
      p_S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 1;
    } else if (I_Steps==steps2+1) {
      p_S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 2;
    } else if (I_Steps==steps2+2) {
      p_S_Steps = I_Steps - 1;
      flagF1 = 2;
      flagF2 = 2;
    } else if (I_Steps>steps2+2) {
      p_S_Steps = I_Steps - 2;
      flagF1 = 2;
      flagF2 = 2;
    }

    if (I_Steps !=0 && last_S_Steps != p_S_Steps) {
      mac_ang_s0->equate(*mac_ang_s1); 
      mac_spd_s0->equate(*mac_spd_s1); 
      eprime_s0->equate(*eprime_s1); 
    }

    vecTemp->equate(*mac_ang_s0); 
    vecTemp->scale(jay); 
    vecTemp->exp(); 
    vecTemp->elementMultiply(*eqprime); 
    eprime_s0->equate(*vecTemp); 
     
    // ---------- CALL i_simu_innerloop(k,p_S_Steps,flagF1): ----------
    int t_trnsmul = timer->createCategory("Dynamic Simulation: Transpose Multiply");
    timer->start(t_trnsmul);
    if (flagF1 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s0)); //MatMultTranspose(prefy11, eprime_s0, curr);
      //transposeMultiply(*prefy11,*eprime_s0,*curr);
    } else if (flagF1 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s0)); //MatMultTranspose(fy11, eprime_s0, curr);
      //transposeMultiply(*fy11,*eprime_s0,*curr);
    } else if (flagF1 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s0)); //MatMultTranspose(posfy11, eprime_s0, curr);
      //transposeMultiply(*posfy11,*eprime_s0,*curr);
    } 
    timer->stop(t_trnsmul);
    
    // ---------- CALL mac_em2(k,p_S_Steps): ----------
    // ---------- pelect: ----------
    curr->conjugate(); 
    vecTemp->equate(*eprime_s0);
    vecTemp->elementMultiply(*curr);
    pelect->equate(*vecTemp); 
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s0); 
    vecTemp->add(-1.0); 
    vecTemp->scale(basrad); 
    dmac_ang_s0->equate(*vecTemp); 
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); 
    dmac_spd_s0->equate(*pmech); 
    dmac_spd_s0->add(*vecTemp, -1.0); 
    vecTemp->equate(*mac_spd_s0); 
    vecTemp->add(-1.0); 
    vecTemp->elementMultiply(*d0); 
    dmac_spd_s0->add(*vecTemp, -1.0); 
    vecTemp->equate(*h); 
    vecTemp->scale(2.0); 
    dmac_spd_s0->elementDivide(*vecTemp); 

    mac_ang_s1->equate(*mac_ang_s0); 
    vecTemp->equate(*dmac_ang_s0);
    mac_ang_s1->add(*dmac_ang_s0, h_sol1); 
    mac_spd_s1->equate(*mac_spd_s0); 
    vecTemp->equate(*dmac_spd_s0);
    mac_spd_s1->add(*vecTemp, h_sol1); 

    vecTemp->equate(*mac_ang_s1); 
    vecTemp->scale(jay); 
    vecTemp->exp(); 
    vecTemp->elementMultiply(*eqprime);
    eprime_s1->equate(*vecTemp); 

    // ---------- CALL i_simu_innerloop2(k,p_S_Steps+1,flagF2): ----------
    timer->start(t_trnsmul);
    if (flagF2 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s1)); //MatMultTranspose(prefy11, eprime_s1, curr);
      //transposeMultiply(*prefy11,*eprime_s1,*curr);
    } else if (flagF2 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s1)); //MatMultTranspose(fy11, eprime_s1, curr);
      //transposeMultiply(*fy11,*eprime_s1,*curr);
    } else if (flagF2 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s1)); //MatMultTranspose(posfy11, eprime_s1, curr);
      //transposeMultiply(*posfy11,*eprime_s1,*curr);
    }
    timer->stop(t_trnsmul);

    // ---------- CALL mac_em2(k,p_S_Steps+1): ---------- 
    // ---------- pelect: ----------
    curr->conjugate(); 
    vecTemp->equate(*eprime_s1);
    vecTemp->elementMultiply(*curr); 
    pelect->equate(*vecTemp); 
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s1); 
    vecTemp->add(-1.0); 
    vecTemp->scale(basrad); 
    dmac_ang_s1->equate(*vecTemp); 
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); 
    dmac_spd_s1->equate(*pmech); 
    dmac_spd_s1->add(*vecTemp, -1.0); 
    vecTemp->equate(*mac_spd_s1); 
    vecTemp->add(-1.0); 
    vecTemp->elementMultiply(*d0); 
    dmac_spd_s1->add(*vecTemp, -1.0); 
    vecTemp->equate(*h); 
    vecTemp->scale(2.0); 
    dmac_spd_s1->elementDivide(*vecTemp); 

    vecTemp->equate(*dmac_ang_s0); 
    vecTemp->add(*dmac_ang_s1); 
    vecTemp->scale(0.5); 
    mac_ang_s1->equate(*mac_ang_s0); 
    mac_ang_s1->add(*vecTemp, h_sol2); 
    vecTemp->equate(*dmac_spd_s0); 
    vecTemp->add(*dmac_spd_s1); 
    vecTemp->scale(0.5); 
    mac_spd_s1->equate(*mac_spd_s0); 
    mac_spd_s1->add(*vecTemp, h_sol2); 
    if (p_generatorWatch && I_Steps%p_watchFrequency == 0) {
      p_factory->setMode(init_mac_ang);
      XMap3.mapToBus(mac_ang_s1);
      p_factory->setMode(init_mac_spd);
      XMap3.mapToBus(mac_spd_s1);
      char tbuf[32];
      sprintf(tbuf,"%8.4f",static_cast<double>(I_Steps)*p_time_step);
      p_generatorIO->header(tbuf);
      p_generatorIO->write("watch");
      p_generatorIO->header("\n");
    }

    // Print to screen
    if (last_S_Steps != p_S_Steps) {
      //sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", p_S_Steps);
      //p_busIO->header(ioBuf);
      //mac_ang_s0->print();  
      //mac_spd_s0->print();  
      //pmech->print();
      //pelect->print();
      //sprintf(ioBuf, "========================End of S_Steps = %d=========================\n\n", p_S_Steps);
      //p_busIO->header(ioBuf);
    }
    if (I_Steps == simu_k) {
      p_factory->setMode(init_mac_ang);
      XMap3.mapToBus(mac_ang_s1);
      p_factory->setMode(init_mac_spd);
      XMap3.mapToBus(mac_spd_s1);
      p_factory->setMode(init_pmech);
      XMap6.mapToBus(pmech);
      p_factory->setMode(init_pelect);
      XMap1.mapToBus(pelect);
      //mac_ang_s1->print();  
      //mac_spd_s1->print();  
      //pmech->print();
      //pelect->print();
    } // End of Print to screen 
    
    last_S_Steps = p_S_Steps;
  }
  closeGeneratorWatchFile();
  timer->stop(t_total);
#else
  timer->start(t_total);
  int t_matset = timer->createCategory("Matrix Setup");
  timer->start(t_matset);
  p_factory->setMode(YBUS);
  gridpack::mapper::FullMatrixMap<DSNetwork> ybusMap(p_network);
  timer->stop(t_matset);

  int t_trans = timer->createCategory("Matrix Transpose");
  // Construct matrix diagY_a using xd and ra extracted from gen data, 
  timer->start(t_matset);
  p_factory->setMode(YA);
  gridpack::mapper::FullMatrixMap<DSNetwork> yaMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> diagY_a = yaMap.mapToMatrix();
  // Convert diagY_a from sparse matrix to dense matrix Y_a so that SuperLU_DIST can solve
  gridpack::math::MatrixStorageType denseType = gridpack::math::Dense;
  boost::shared_ptr<gridpack::math::Matrix> Y_a(gridpack::math::storageType(*diagY_a, denseType));
  timer->stop(t_matset);

  int t_matmul = timer->createCategory("Matrix Multiply");
  // Construct Y_c as a dense matrix
  timer->start(t_matset);
  p_factory->setMode(YC);
  gridpack::mapper::FullMatrixMap<DSNetwork> cMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Y_cDense = cMap.mapToMatrix(true);

  // Construct Y_b
  p_factory->setMode(YB);
  gridpack::mapper::FullMatrixMap<DSNetwork> bMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> Y_b = bMap.mapToMatrix();
  timer->stop(t_matset);
  
  timer->start(t_matset);
  p_factory->setMode(updateYbus);
  boost::shared_ptr<gridpack::math::Matrix> prefy11ybus = ybusMap.mapToMatrix();
  timer->stop(t_matset);

  // Solve linear equations ybus * X = Y_c
  int t_solve = timer->createCategory("Solve Linear Equation");
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver1(*prefy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> prefy11X(solver1.solve(*Y_cDense));
  timer->stop(t_solve);
  
  //-----------------------------------------------------------------------
  // Compute prefy11
  //-----------------------------------------------------------------------
  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> prefy11(multiply(*Y_b, *prefy11X)); 
  timer->stop(t_matmul);
  // Update prefy11: prefy11 = Y_a + prefy11
  prefy11->add(*Y_a);

  //-----------------------------------------------------------------------
  // Compute fy11
  // Update ybus values at beginning of fault
  //-----------------------------------------------------------------------
  // Read the switch info from fault. Event from input.xml
  int sw2_2 = fault.from_idx - 1;
  int sw3_2 = fault.to_idx - 1;

  boost::shared_ptr<gridpack::math::Matrix> fy11ybus(prefy11ybus->clone());
  gridpack::ComplexType x(0.0, -1e7);
  timer->start(t_matset);
  p_factory->setEvent(fault);
  p_factory->setMode(onFY);
  ybusMap.overwriteMatrix(fy11ybus);
  timer->stop(t_matset);

  // Solve linear equations of fy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver2(*fy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> fy11X(solver2.solve(*Y_cDense)); 
  timer->stop(t_solve);
  
  // Form reduced admittance matrix fy11: fy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> fy11(multiply(*Y_b, *fy11X)); 
  timer->stop(t_matmul);
  // Update fy11: fy11 = Y_a + fy11
  fy11->add(*Y_a);

  //-----------------------------------------------------------------------
  // Compute posfy11
  // Update ybus values at clear fault stage
  //-----------------------------------------------------------------------
  // Get the updating factor for posfy11 stage ybus
  timer->start(t_matset);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus(prefy11ybus->clone());
  timer->stop(t_matset);
  timer->start(t_matset);
  p_factory->setMode(posFY);
  ybusMap.incrementMatrix(posfy11ybus);
  timer->stop(t_matset);
    
  // Solve linear equations of posfy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver3(*posfy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> posfy11X(solver3.solve(*Y_cDense)); 
  timer->stop(t_solve);
  
  // Form reduced admittance matrix posfy11: posfy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> posfy11(multiply(*Y_b, *posfy11X)); 
  timer->stop(t_matmul);
  // Update posfy11: posfy11 = Y_a + posfy11
  posfy11->add(*Y_a);

  //-----------------------------------------------------------------------
  // Integration implementation (Modified Euler Method)
  //-----------------------------------------------------------------------
 
  p_factory->setDSParams();
  p_factory->setMode(Eprime0);
  // create vectors used in solution method
  gridpack::mapper::BusVectorMap<DSNetwork> Emap(p_network);
  boost::shared_ptr<gridpack::math::Vector> eprime_s0 = Emap.mapToVector();
  boost::shared_ptr<gridpack::math::Vector> eprime_s1(eprime_s0->clone());
  boost::shared_ptr<gridpack::math::Vector> curr(eprime_s0->clone());

  timer->start(t_trans);
  boost::shared_ptr<gridpack::math::Matrix> trans_prefy11(transpose(*prefy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_fy11(transpose(*fy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_posfy11(transpose(*posfy11));
  timer->stop(t_trans);

  // Simulation related variables
  int simu_k;
  int t_step[20];
  double t_width[20];
  int flagF;
  int S_Steps;
  int last_S_Steps;
  int steps3, steps2, steps1;
  double h_sol1, h_sol2;
  int flagF1, flagF2;
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

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagF1 = 0; 
  flagF2 = 0; 
  S_Steps = 1; 
  last_S_Steps = -1;

  if (p_generatorWatch) {
    p_generatorIO->header("t");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    p_generatorIO->header("\n");
  }
  for (I_Steps = 0; I_Steps < simu_k+1; I_Steps++) {
    if (I_Steps < steps1) {
      S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 0;
    } else if (I_Steps == steps1) { 
      S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 1;
    } else if (I_Steps == steps1+1) {
      S_Steps = I_Steps;
      flagF1 = 1;
      flagF2 = 1;
    } else if ((I_Steps>steps1+1) && (I_Steps<steps2+1)) {
      S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 1;
    } else if (I_Steps==steps2+1) {
      S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 2;
    } else if (I_Steps==steps2+2) {
      S_Steps = I_Steps - 1;
      flagF1 = 2;
      flagF2 = 2;
    } else if (I_Steps>steps2+2) {
      S_Steps = I_Steps - 2;
      flagF1 = 2;
      flagF2 = 2;
    }

    if (I_Steps !=0 && last_S_Steps != S_Steps) {
      p_factory->initDSStep(false);
    } else {
      p_factory->initDSStep(true);
    }
    p_factory->setMode(Eprime0);
    Emap.mapToVector(eprime_s0);
     
    // ---------- CALL i_simu_innerloop(k,S_Steps,flagF1): ----------
    int t_trnsmul = timer->createCategory("Transpose Multiply");
    timer->start(t_trnsmul);
    if (flagF1 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s0)); //MatMultTranspose(prefy11, eprime_s0, curr);
      //transposeMultiply(*prefy11,*eprime_s0,*curr);
    } else if (flagF1 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s0)); //MatMultTranspose(fy11, eprime_s0, curr);
      //transposeMultiply(*fy11,*eprime_s0,*curr);
    } else if (flagF1 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s0)); //MatMultTranspose(posfy11, eprime_s0, curr);
      //transposeMultiply(*posfy11,*eprime_s0,*curr);
    } 
    timer->stop(t_trnsmul);
    
    p_factory->setMode(Current);
    Emap.mapToBus(curr);
    p_factory->predDSStep(h_sol1);
    p_factory->setMode(Eprime1);
    Emap.mapToVector(eprime_s1);

    // ---------- CALL i_simu_innerloop2(k,S_Steps+1,flagF2): ----------
    timer->start(t_trnsmul);
    if (flagF2 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s1)); //MatMultTranspose(prefy11, eprime_s1, curr);
      //transposeMultiply(*prefy11,*eprime_s1,*curr);
    } else if (flagF2 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s1)); //MatMultTranspose(fy11, eprime_s1, curr);
      //transposeMultiply(*fy11,*eprime_s1,*curr);
    } else if (flagF2 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s1)); //MatMultTranspose(posfy11, eprime_s1, curr);
      //transposeMultiply(*posfy11,*eprime_s1,*curr);
    }
    timer->stop(t_trnsmul);

    p_factory->setMode(Current);
    Emap.mapToBus(curr);
    p_factory->corrDSStep(h_sol2);

    if (p_generatorWatch && I_Steps%p_watchFrequency == 0) {
      char tbuf[32];
      sprintf(tbuf,"%8.4f",static_cast<double>(I_Steps)*p_time_step);
      p_generatorIO->header(tbuf);
      p_generatorIO->write("watch");
      p_generatorIO->header("\n");
    }
    // Print to screen
    if (I_Steps == simu_k) {
      char ioBuf[128];
      sprintf(ioBuf, "\n========================S_Steps = %d=========================\n",
          S_Steps+1);
      p_busIO->header(ioBuf);
      sprintf(ioBuf, "\n         Bus ID     Generator ID"
          "    mac_ang         mac_spd         mech            elect\n\n");
      p_busIO->header(ioBuf);
      p_busIO->write();
      sprintf(ioBuf, "\n========================End of S_Steps = %d=========================\n\n",
          S_Steps+1);
      p_busIO->header(ioBuf);
    } // End of Print to screen 
    
    last_S_Steps = S_Steps;
  }
  closeGeneratorWatchFile();
  timer->stop(t_total);
#endif
}
/**
 * Write out final results of dynamic simulation calculation to
 * standard output
 */
void gridpack::dynamic_simulation_r::DSAppModule::write()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);
  char ioBuf[128];
  sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", p_S_Steps+1);
  p_busIO->header(ioBuf);
  sprintf(ioBuf, "\n         Bus ID     Generator ID"
      "    mac_ang         mac_spd         mech            elect\n\n");
  p_busIO->header(ioBuf);
  p_busIO->write();
  sprintf(ioBuf, "\n========================End of S_Steps = %d=========================\n\n", p_S_Steps+1);
  p_busIO->header(ioBuf);
  timer->stop(t_total);
}

/**
 * Read in generators that should be monitored during simulation
 */
void gridpack::dynamic_simulation_r::DSAppModule::setGeneratorWatch()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
  if (!cursor->get("generatorWatchFrequency",&p_watchFrequency)) {
    p_watchFrequency = 1;
  }
  char buf[128];
  cursor = p_config->getCursor("Configuration.Dynamic_simulation.generatorWatch");
  gridpack::utility::Configuration::ChildCursors generators;
  if (cursor) cursor->children(generators);
  int i, j, idx, id, len;
  int ncnt = generators.size();
  std::string generator, tag, clean_tag;
  gridpack::dynamic_simulation_r::DSBus *bus;
  if (ncnt > 0) p_busIO->header("Monitoring generators:\n");
  for (i=0; i<ncnt; i++) {
    // Parse contents of "generator" field to get bus ID and generator tag
    generators[i]->get("busID",&id);
    generators[i]->get("generatorID",&tag);
    gridpack::utility::StringUtils util;
    clean_tag = util.clean2Char(tag);
    // Find local bus indices for generator
    std::vector<int> local_ids = p_network->getLocalBusIndices(id);
    for (j=0; j<local_ids.size(); j++) {
      bus = dynamic_cast<gridpack::dynamic_simulation_r::DSBus*>
        (p_network->getBus(local_ids[j]).get());
      bus->setWatch(clean_tag,true);
    }
    sprintf(buf,"  Bus: %8d Generator ID: %2s\n",id,clean_tag.c_str());
    p_busIO->header(buf);
  }
  if (ncnt > 0) {
    p_generatorWatch = true;
    sprintf(buf,"Generator Watch Frequency: %d\n",p_watchFrequency);
    p_busIO->header(buf);
    openGeneratorWatchFile();
  }
}

/**
 * Close file contain generator watch results
 */
void gridpack::dynamic_simulation_r::DSAppModule::openGeneratorWatchFile()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");
#ifndef USEX_GOSS
  std::string filename;
  if (cursor->get("generatorWatchFileName",&filename)) {
    p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSNetwork>(128,
          p_network));
    p_generatorIO->open(filename.c_str());
  } else {
    p_busIO->header("No Generator Watch File Name Found\n");
    p_generatorWatch = false;
  }
#else
  std::string topic, URI, username, passwd;
  bool ok = true;
  ok = ok && cursor->get("channelTopic",&topic);
  ok = ok && cursor->get("channelURI",&URI);
  ok = ok && cursor->get("username",&username);
  ok = ok && cursor->get("password",&passwd);
  if (ok) {
    p_generatorIO.reset(new gridpack::serial_io::SerialBusIO<DSNetwork>(128,
          p_network));
    p_generatorIO->openChannel(topic.c_str());
  } else {
    p_busIO->header("Unable to open channel\n");
    p_generatorWatch = false;
  }
#endif
}

/**
 * Close file contain generator watch results
 */
void gridpack::dynamic_simulation_r::DSAppModule::closeGeneratorWatchFile()
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
 * Utility function to convert faults that are in event list into
 * internal data structure that can be used by code
 */
void gridpack::dynamic_simulation_r::DSAppModule::setFaultEvents()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation.faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  if (cursor) cursor->children(events);
  int size = events.size();
  int idx;
  // Parse fault events
  for (idx=0; idx<size; idx++) {
    gridpack::dynamic_simulation_r::DSBranch::Event event;
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
      p_faults.push_back(event);
    }
  }
}
