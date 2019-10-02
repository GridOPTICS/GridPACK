/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsim.hpp
 * @author Shrirang Abhyankar
 * @date   Feb 06, 2019
 *
 * @brief  Header file for the dynamic simulation application
 *
 *
 */
// -------------------------------------------------------------

#ifndef _dsim_h_
#define _dsim_h_
     
#include <dsimutils.hpp>
#include <dsimnetwork.hpp>
#include <dsimfactory.hpp>
#include <gridpack/math/dae_solver.hpp>

class DSim
{
   public:
  /**
     Basic Constructor
  */
  DSim(void); 

  /**
   * Basic constructor with commmunicator argument
   * @param comm communicator that application object is restricted to
   */
  DSim(gridpack::parallel::Communicator comm);

  /**
   * Basic destructor
   */
  ~DSim(void);

  /**
   *
   */
  int rank() {p_comm.rank(); };

  int size() {p_comm.size(); };

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
    p_VecMapper->mapToBus(X);

    // Push current values in Xdot vector back into network components
    p_factory->setMode(XDOTVECTOBUS);
    p_VecMapper->mapToBus(Xdot);

    // Update ghost buses
    p_network->updateBuses();

    // Evaluate the DAE Jacobian
    J.zero();
    p_factory->setMode(RESIDUAL_EVAL);
    p_MatMapper->mapToMatrix(J);
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
    //    F.print();
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
    p_factory->setMode(FAULT_EVAL);
    p_MatMapper->mapToMatrix(J);
    
  }

  private:

  // Communicator
  gridpack::parallel::Communicator p_comm;

  // Simulation parameters
  DSimParams p_simparams;

  // Profiler
  DSimProfiler p_profiler;

  // Configuration
  gridpack::utility::Configuration *p_config;
  gridpack::utility::Configuration::CursorPtr p_configcursor;
  char p_configfile[256];


  // Set up called
  int p_isSetUp;

  // Network pointer
  boost::shared_ptr<DSimNetwork> p_network;
  
  // Factory
  DSimFactory *p_factory;

  // Mappers for creating vectors and matrices
  gridpack::mapper::BusVectorMap<DSimNetwork> *p_VecMapper;
  gridpack::mapper::FullMatrixMap<DSimNetwork> *p_MatMapper;

  boost::shared_ptr<gridpack::math::Vector> p_X; // Solution vector
  boost::shared_ptr<gridpack::math::Vector> p_R; // Residual vector
  boost::shared_ptr<gridpack::math::Matrix> p_J; // Jacobian matrix

  // DAE solver
  gridpack::math::DAESolver *p_daesolver;

  // Nonlinear solver for handling discontinuities
  gridpack::math::NonlinearSolver *p_nlsolver;
};



#endif
