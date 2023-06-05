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
     
#include <emtutils.hpp>
#include <emtnetwork.hpp>
#include <emtfactory.hpp>
#include <gridpack/math/dae_solver.hpp>

class Emt
{
public:
  
  // Let's try to avoid some typing
  
  typedef gridpack::math::DAESolver DAESolver;
  typedef DAESolver::VectorType VectorType;
  typedef DAESolver::Event Event;
  typedef DAESolver::EventPtr EventPtr;
  typedef DAESolver::EventManager EventManager;
  typedef DAESolver::EventManagerPtr EventManagerPtr;
  
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
   *
   */
  int rank() {return p_comm.rank(); };

  int size() {return p_comm.size(); };

  void setconfigurationfile(const char*);

  void setup(void);

  void initialize(void);

  void solve(void);

  void readnetworkdatafromconfig();

  /// Build the DAE Jacobian
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, 
		   const gridpack::math::Vector& Xdot, 
		   const double& shift, gridpack::math::Matrix& J)
  {
    p_factory->setTSshift(shift);
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);
    //    p_VecMapper->mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory->setMode(XDOTVECTOBUS);
    p_VecMapper->mapToBus(Xdot);

    // Update ghost buses
    //    p_network->updateBuses();

    // Evaluate the DAE Jacobian
    //    J.zero();
    p_factory->setMode(RESIDUAL_EVAL);
    p_MatMapper->mapToMatrix(J);
    J.ready();
  }

  /// Build the DAE RHS function
  void operator() (const double& time, 
		   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
		   gridpack::math::Vector& F)
  {
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);
    p_VecMapper->mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory->setMode(XDOTVECTOBUS);
    p_VecMapper->mapToBus(Xdot);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the residual f(x) - xdot
    p_factory->setMode(RESIDUAL_EVAL);
    p_VecMapper->mapToVector(F);
    F.ready();
    /*    printf("F.print():\n");
    F.print();
    exit(0);
    */
  }

  // Build the residual for the nonlinear solver at tfaulton and tfaultoff
  void  operator() (const gridpack::math::Vector& X, gridpack::math::Vector& F)
  {
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);
    p_VecMapper->mapToBus(X);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the residual f(x) - xdot
    p_factory->setMode(FAULT_EVAL);
    p_VecMapper->mapToVector(F);
    F.ready();
    //    F.print();
  }

  // Build the Jacobian for the nonlinear solver at tfaulton or tfaultoff
  void  operator() (const gridpack::math::Vector& X,gridpack::math::Matrix& J)
  {
    // Push current values in X vector back into network components
    p_factory->setMode(XVECTOBUS);
    p_VecMapper->mapToBus(X);

    // Update ghost buses
    p_network->updateBuses();

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
  boost::shared_ptr<EmtNetwork> p_network;
  
  // Factory
  EmtFactory *p_factory;

  // Mappers for creating vectors and matrices
  gridpack::mapper::BusVectorMap<EmtNetwork> *p_VecMapper;
  gridpack::mapper::FullMatrixMap<EmtNetwork> *p_MatMapper;

  boost::shared_ptr<gridpack::math::Vector> p_X; // Solution vector
  boost::shared_ptr<gridpack::math::Vector> p_R; // Residual vector
  boost::shared_ptr<gridpack::math::Matrix> p_J; // Jacobian matrix

  // DAE solver
  gridpack::math::DAESolver *p_daesolver;

  // Nonlinear solver for handling discontinuities
  gridpack::math::NonlinearSolver *p_nlsolver;

  /// Does the network need a resolve
  bool p_resolve;

  /// These class needs to see inside Emt
  friend class EmtTimedFaultEvent;
  friend class EmtEventManager;
};




#endif
