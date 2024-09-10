/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emt.hpp
 *
 * @brief  Header file for the EMT application
 *
 *
 */
// -------------------------------------------------------------

#ifndef _emt_h_
#define _emt_h_
     
#include <gridpack/applications/modules/emt/emtutils.hpp>
#include <gridpack/applications/modules/emt/emtnetwork.hpp>
#include <gridpack/applications/modules/emt/emtfactory.hpp>
#include <gridpack/math/dae_solver.hpp>

class Emt
{
public:
  
  // Typedef for some objects
  
  typedef gridpack::math::RealDAESolver RealDAESolver;
  typedef RealDAESolver::VectorType VectorType;
  typedef RealDAESolver::Event Event;
  typedef RealDAESolver::EventPtr EventPtr;
  typedef RealDAESolver::EventManager EventManager;
  typedef RealDAESolver::EventManagerPtr EventManagerPtr;
  
  /**
     Basic Constructor
  */
  Emt(void); 

  /**
   * Basic constructor with commmunicator argument
   * @param comm communicator that application object is restricted to
   */
  Emt(gridpack::parallel::Communicator comm);

  /**
   * Basic destructor
   */
  ~Emt(void);

  /**
   *   Get the communicator rank
   */
  int rank() {return p_comm.rank(); };

  /**
   *   Get the communicator size
   */

  int size() {return p_comm.size(); };

  /**
   *   Set the configuration file
   */

  void setconfigurationfile(const char*);

  void setup(void);

  void initialize(void);

  void solve(void);

  void readnetworkdatafromconfig();

  void solvepowerflow(void);

  void transferPFtoEMT(boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network,boost::shared_ptr<EmtNetwork> emt_network);

  /// Build the DAE Jacobian
  void operator() (const double& time, 
		   const gridpack::math::RealVector& X, 
		   const gridpack::math::RealVector& Xdot, 
		   const double& shift, gridpack::math::RealMatrix& J)
  {
    double dt = p_daesolver->gettimestep();
    
    p_factory->setTSshift(shift);
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);
    //    p_VecMapper->mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory->setMode(XDOTVECTOBUS);

    p_VecMapper->mapToNetwork(Xdot);

    // Update ghost buses
    //    emt_network->updateBuses();

    // Evaluate the DAE Jacobian
    //    J.zero();
    p_factory->setMode(RESIDUAL_EVAL);
    p_MatMapper->mapToMatrix(J);
    J.ready();

  }

  /// Build the DAE RHS function
  void operator() (const double& time, 
		   const gridpack::math::RealVector& X, const gridpack::math::RealVector& Xdot, 
		   gridpack::math::RealVector& F)
  {
    p_factory->setTime(time);

    // Push current values in Xdot vector back into network components
    p_factory->setMode(XDOTVECTOBUS);

    p_VecMapper->mapToNetwork(Xdot);

    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);

    p_VecMapper->mapToNetwork(X);

    // Update ghost buses and branches
    emt_network->updateBuses();
    emt_network->updateBranches();

    // Evaluate the residual f(x) - xdot
    p_factory->setMode(RESIDUAL_EVAL);
    F.zero();
    p_VecMapper->mapToVector(F);
    F.ready();
    X.print(0);
    //F.print(0);
    
  }

  /// Prestep function
  void operator() (const double& time,
		   const double& timestep,
		   const gridpack::math::RealVector& X)
  {
    p_factory->setTime(time);

    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);

    p_VecMapper->mapToNetwork(X);

    // Update ghost buses and branches
    emt_network->updateBuses();
    emt_network->updateBranches();

    // Prestep
    p_factory->preStep(time,timestep);
    
  }

    /// Poststep function
  void operator() (const double& time,
		   const gridpack::math::RealVector& X)
  {
    int nsteps = p_daesolver->getstepnumber();
    if(nsteps % reuseprecon_nsteps == 0) {
      p_daesolver->reusepreconditioner(-2); // Update preconditioner at next step
    } else {
      p_daesolver->reusepreconditioner(-1); // Reuse preconditioner
    }
    p_factory->postStep(time);

    save_output(time);
  }

  // Build the residual for the nonlinear solver at tfaulton and tfaultoff
  void  operator() (const gridpack::math::RealVector& X, gridpack::math::RealVector& F)
  {
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);

    p_VecMapper->mapToNetwork(X);

    // Update ghost buses
    emt_network->updateBuses();

    // Evaluate the residual f(x) - xdot
    p_factory->setMode(FAULT_EVAL);
    p_VecMapper->mapToVector(F);
    F.ready();

  }

  // Build the Jacobian for the nonlinear solver at tfaulton or tfaultoff
  void  operator() (const gridpack::math::RealVector& X,gridpack::math::RealMatrix& J)
  {
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);

    p_VecMapper->mapToNetwork(X);

    // Update ghost buses
    emt_network->updateBuses();

    // Evaluate the fault residual Jacobian
    J.zero();
    p_factory->setMode(FAULT_EVAL);
    p_MatMapper->mapToMatrix(J);
    J.ready();

  }

  private:

  // Communicator
  gridpack::parallel::Communicator p_comm;

  // Simulation parameters
  EmtParams p_simparams;

  // Profiler
  EmtProfiler p_profiler;

  // Configuration
  gridpack::utility::Configuration *p_config;
  gridpack::utility::Configuration::CursorPtr p_configcursor;
  char p_configfile[256];


  // Set up called
  int p_isSetUp;

  // Network pointer
  boost::shared_ptr<EmtNetwork> emt_network;

  // Power flow application pointer
  gridpack::powerflow::PFAppModule *p_pfapp;
  
  // Factory
  EmtFactory *p_factory;

  // Save output function
  void save_output(const double& time);

  // Save output no input
  void save_output() {
    save_output(0.0);
  }

  // Mappers for creating vectors and matrices
  gridpack::mapper::GenVectorMap<EmtNetwork,gridpack::RealType,gridpack::math::RealVector> *p_VecMapper;
  gridpack::mapper::GenMatrixMap<EmtNetwork,gridpack::RealType,gridpack::math::RealMatrix> *p_MatMapper;

  boost::shared_ptr<gridpack::math::RealVector> p_X; // Solution vector
  boost::shared_ptr<gridpack::math::RealVector> p_R; // Residual vector
  boost::shared_ptr<gridpack::math::RealMatrix> p_J; // Jacobian matrix

  // DAE solver
  gridpack::math::RealDAESolver *p_daesolver;

  // Integration algorithm for machines
  EMTMachineIntegrationType p_emtmachineintegrationtype;

  int reuseprecon_nsteps; // Reuse preconditioner for nsteps

  void setMonitors(gridpack::utility::Configuration::CursorPtr);

  // Output requested by user?
  bool p_saveoutput;

  // Output file name
  std::string p_monitorfile;

  // File pointer
  FILE *fp_monitor;

  std::vector<BaseEMTGenModel*> monitored_gens;
  std::vector<EmtBus*> monitored_buses;

  char output_string[512]; // Output is written to this string

  /// These class needs to see inside Emt
  friend class EmtTimedFaultEvent;
  friend class EmtEventManager;
};

#endif
