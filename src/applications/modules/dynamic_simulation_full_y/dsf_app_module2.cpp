#include "dsf_app_module.hpp"
#include <iostream>
#include <string>
#include <vector>
#include "gridpack/utilities/string_utils.hpp"

/*
  Get the time-step
*/
double gridpack::dynamic_simulation::DSFullApp::getTimeStep()
{
  return p_time_step;
}

/*
  Set the time-step
*/
void gridpack::dynamic_simulation::DSFullApp::setTimeStep(double time_step)
{
  p_time_step = time_step;
}

/*
  Set simulation end time
*/
void gridpack::dynamic_simulation::DSFullApp::setFinalTime(double final_time)
{
  p_sim_time = final_time;
}

/*
  Get simulation end time
*/
double gridpack::dynamic_simulation::DSFullApp::getFinalTime()
{
  return p_sim_time;
}

/*
  Get current time
*/
double gridpack::dynamic_simulation::DSFullApp::getCurrentTime()
{
  return p_current_time;
}

/**
 * Reset data structures
 */
void gridpack::dynamic_simulation::DSFullApp::reset()
{
  orgYbus.reset();
  ybusyl.reset();
  ybuspg.reset();
  ybus_jxd.reset();
  ybus.reset();

  ybusMap_sptr.reset();
  nbusMap_sptr.reset();

  INorton_full.reset();
  INorton_full_chk.reset();
  volt_full.reset();

  solver_sptr.reset();

  p_factory->load();
}


/*
  solve network equations
  predcorrflag = 0 => Predictor stage
  predcorrflag = 1 => Corrector stage
*/
bool gridpack::dynamic_simulation::DSFullApp::solveNetwork(int predcorrflag)
{
  bool converged = false;
  
  //  volt_full->zero();

  int its=0;
  //bool p_iterative_network_debug = false;
  if (p_biterative_solve_network) {
    /* Iterative network solution */
    its = 0;
    while (!converged &&  its <= MAX_ITR_NO ) {
      /* Solve network equations for volt_full */
      solver_sptr->solve(*INorton_full, *volt_full);
      
      /* Copy over INorton_full vector (INorton_{i-1}) */
      INorton_full_chk->equate(*INorton_full);
      
      /* Push voltages on to buses */
      nbusMap_sptr->mapToBus(volt_full);
      p_factory->setVolt(false);

      getCurrent(predcorrflag);
      
      /* Compute current change (INorton_i - INorton_{i-1}) */
      INorton_full_chk->add(*INorton_full, -1.0);
      
      /* Norm of current change */
      max_INorton_full=abs(INorton_full_chk->normInfinity());
      
      /* Convergence check */
      if (max_INorton_full < ITER_TOL) {
	/* Converged */
	converged = true;
      } else {
	its += 1;
      }
    }// end of while
  }else {// p_biterative_solve_network = false
    /* Non-iterative solution */
   solver_sptr->solve(*INorton_full, *volt_full);
   /* Push voltage to buses */
   nbusMap_sptr->mapToBus(volt_full);
   p_factory->setVolt(false);

   converged = true;
 } // end of if (p_biterative_solve_network)

  return converged;
}

/*
  Update Norton current injected in the network
  predcorrflag = 0 => Predictor stage
  predcorrflag = 1 => Corrector stage
*/
void gridpack::dynamic_simulation::DSFullApp::getCurrent(int predcorrflag)
{
  if(predcorrflag == 0) { /* Predictor */
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->predictor_currentInjection(false);
    } else {
      p_factory->predictor_currentInjection(true);
    }    
  } else { /* Corrector */
    if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
      p_factory->corrector_currentInjection(false);
    } else {
      p_factory->corrector_currentInjection(true);
    }
  }
  /* Norton current injection vector */
  p_factory->setMode(make_INorton_full);
  nbusMap_sptr->mapToVector(INorton_full);  
}
/**
 * initialization before the time step integration starts 
 */
void gridpack::dynamic_simulation::DSFullApp::setup()
{
  p_current_time = 0.0;
  
  // Get cursor for setting solver options
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");

  // Get events from input file
  p_events = getEvents();

  p_factory->setMode(YBUS);
  
  ybusMap_sptr.reset(new gridpack::mapper::FullMatrixMap<DSFullNetwork> (p_network));
  orgYbus = ybusMap_sptr->mapToMatrix();
  
  // Form constant impedance load admittance yl for all buses and add it to
  // system Y matrix: ybus = ybus + yl
  p_factory->setMode(YL);
  ybusyl = ybusMap_sptr->mapToMatrix();

  p_factory->setMode(PG);
  ybuspg = ybusMap_sptr->mapToMatrix();

  // Add j*Xd' to system Y matrix:
  // Extract appropriate xdprime and xdpprime from machine data
  p_factory->setMode(jxd);
  ybus_jxd = ybusMap_sptr->mapToMatrix();
  
  // Add dynamic load impedance to system Y matrix:
  p_factory->setMode(YDYNLOAD);
  ybus = ybusMap_sptr->mapToMatrix();

  // Initialize vectors for integration 
  p_factory->initDSVect(p_time_step);
  
  p_factory->setGeneratorObPowerBaseFlag(p_generator_observationpower_systembase);

  /* Create mapper and vectors */
  ngenMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork> (p_network));
  
  p_insecureAt = -1;

  p_factory->setMode(make_INorton_full);
  nbusMap_sptr.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork>(p_network));

  /* Norton current */
  INorton_full = nbusMap_sptr->mapToVector();

  /* Copy of Norton current vector used during network iterations */
  INorton_full_chk = nbusMap_sptr->mapToVector();
  max_INorton_full = 0.0;

  /* Voltage vector */
  volt_full.reset(INorton_full->clone());

  /* Linear solver */
  solver_sptr.reset(new gridpack::math::LinearSolver (*ybus));
  solver_sptr->configure(cursor);

  
  if (!p_suppress_watch_files) {
    /* Create CSV file header */
    if (p_generatorWatch) p_generatorIO->header("t");
    if (p_generatorWatch) p_generatorIO->write("watch_header");
    if (p_generatorWatch) p_generatorIO->header("\n");

    if (p_loadWatch) p_loadIO->header("t");
    if (p_loadWatch) p_loadIO->write("load_watch_header");
    if (p_loadWatch) p_loadIO->header("\n");
  }

  p_frequencyOK = true;
}

/**
 * Execute only one simulation time step 
 */
void gridpack::dynamic_simulation::DSFullApp::runonestep()
{
  bool converged;
  
  S_Steps = Simu_Current_Step;

  /* Predictor current injection */
  getCurrent(0);

  /* Solve Network equations */
  converged = solveNetwork(0);
  
  if ( Simu_Current_Step==0 ) {
    //printf("enter the initial update oldbusvoltage, Timestep: %d \n", Simu_Current_Step);
    p_factory->updateoldbusvoltage(); //renke add, first timestep, copy volt_full to volt_full_old
  }

  /* Update frequency */
  p_factory->updateBusFreq(p_time_step);
	
  std::vector <double> vwideareafreqs;
  vwideareafreqs = p_factory->grabWideAreaFreq();
  int tmp = vwideareafreqs.size();
  double widearea_deltafreq = vwideareafreqs[tmp-1];

  bool flagBus = p_factory->updateBusRelay(false, p_time_step);
  bool flagBranch = p_factory->updateBranchRelay(false, p_time_step);
	
  // update dynamic load internal relay functions here
  p_factory->dynamicload_post_process(p_time_step, false);
    
  // if bus relay trips, modify the corresponding Ymatrix
  if (flagBus) {
    p_factory->setMode(bus_relay);
    ybusMap_sptr->overwriteMatrix(ybus);
  }
	
  // if branch relay trips, modify the corresponding Ymatrix
  if (flagBranch) {
    p_factory->setMode(branch_relay);
    ybusMap_sptr->incrementMatrix(ybus);
  }
	
  // Update old voltage (??)
  p_factory->updateoldbusvoltage();

  /* Predictor */
  if (Simu_Current_Step !=0 && last_S_Steps != S_Steps) {
    p_factory->predictor(p_time_step, false);
  } else { 
    p_factory->predictor(p_time_step, true);
  }

  /* Network current injection */
  getCurrent(1);

  /* Solve network */
  converged = solveNetwork(1);

  /* Update frequency */
  p_factory->updateBusFreq(p_time_step);

  /* Correct update */
  if (last_S_Steps != S_Steps) {
    p_factory->corrector(p_time_step, false);
  } else {
    p_factory->corrector(p_time_step, true);
  }

  if (p_generatorWatch && Simu_Current_Step%p_generatorWatchFrequency == 0) {
    char tbuf[32];
    if (!p_suppress_watch_files) {
        sprintf(tbuf,"%8.4f",p_current_time);
        if (p_generatorWatch) p_generatorIO->header(tbuf);
        if (p_generatorWatch) p_generatorIO->write("watch");
        if (p_generatorWatch) p_generatorIO->header("\n");
    }
  }

  if (p_loadWatch && Simu_Current_Step%p_loadWatchFrequency == 0) {
    char tbuf[32];
    if (!p_suppress_watch_files) {
      sprintf(tbuf,"%8.4f",p_current_time);
      if (p_loadWatch) p_loadIO->header(tbuf);
      if (p_loadWatch) p_loadIO->write("load_watch");
      if (p_loadWatch) p_loadIO->header("\n");
    }
  }
  
  if ((!p_factory->securityCheck()) && p_insecureAt == -1)  
    p_insecureAt = Simu_Current_Step;

  last_S_Steps = S_Steps;
  
  if (p_monitorGenerators) {
    p_frequencyOK = p_frequencyOK && checkFrequency(0.5,p_current_time);
    if (!p_frequencyOK) Simu_Current_Step = simu_total_steps;
  }

  /* Update steps and current time */
  Simu_Current_Step++;
  p_current_time = p_current_time + p_time_step;
}

/**
   setLineStatus - Sets the line status and updates the associated
   branch and bus objects. 
   
   @param: from_idx - from bus number
   @param: to_idx - to bus number
   @param: ckt_id - circuit id
   @param: status - new line status
   
   Note: This method is called by handleEvents method to
   update the branch status and update the bus/branch
   objects. It sets up values in the bus and branch objects
   so that incrementMatrix method called on the network Ybus
   uses these values to remove the branch contributions from
   the Y-bus matrix
**/
void gridpack::dynamic_simulation::DSFullApp::setLineStatus(int from_idx, int to_idx, std::string ckt_id, int status)
{
  /* Get Branch */
  std::vector<int> vec_branchintidx;
  vec_branchintidx = p_network->getLocalBranchIndices(from_idx,to_idx);
  int ibr, nbr;
  gridpack::dynamic_simulation::DSFullBranch *pbranch;	
  nbr = vec_branchintidx.size();
  pbranch = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>
    (p_network->getBranch(vec_branchintidx[0]).get());
  if(pbranch) {
    pbranch->setLineStatus(ckt_id,status);
  }
}

/**
   setGenStatus - Sets the generator status and updates the associated
   bus objects. 

   @param: bus_idx - bus number
   @param: gen_id - generator id
   @param: status - new generator status
   
   Note: This method is called by handleEvents method to
   update the generator status and update the bus
   object. It sets up values in the bus objects
   so that incrementMatrix method called on the network Ybus
   uses these values to remove the generator contributions from
   the Y-bus matrix
**/
void gridpack::dynamic_simulation::DSFullApp::setGenStatus(int bus_idx, std::string gen_id, int status)
{
  std::vector<int> bus_internal_idx;
  gridpack::dynamic_simulation::DSFullBus *bus;
  bus_internal_idx = p_network->getLocalBusIndices(bus_idx);
  if(bus_internal_idx.size()) {
    bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_network->getBus(bus_internal_idx[0]).get());
    bus->setGenStatus(gen_id,status);
  }
}



/**
 ** Run till time tend
**/
void gridpack::dynamic_simulation::DSFullApp::run(double tend)
{
  while(fabs(tend - p_current_time) > 1e-6) {

    // Process events
    handleEvents();
    
    // advance one step
    runonestep();

    if(!p_comm.rank())
      printf("Time = %5.4f\n",p_current_time);
  }
}

/**
 ** Run till end time
**/
void gridpack::dynamic_simulation::DSFullApp::run()
{
  run(p_sim_time);
}
