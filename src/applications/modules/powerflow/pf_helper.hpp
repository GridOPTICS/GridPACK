// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   pf_helper.cpp
 * @author William A. Perkins
 * @date   2014-11-25 07:12:27 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <ga++.h>

#include "pf_factory_module.hpp"


namespace gridpack {
namespace powerflow {

/**
 * A helper functor for the powerflow solver.
 * This is a utility functor that provides functions that build the
 * Jacobian and RHS from the network.  
 * 
 */
struct PFSolverHelper 
  : private utility::Uncopyable
{

  // The powerflow factory 
  boost::shared_ptr<PFFactoryModule> p_factory;

  // The powerflow network controlled by ::p_factory
  boost::shared_ptr<PFNetwork> p_network;

  // A place to build/store the Jacobian
  boost::shared_ptr<math::RealMatrix> J;

  // The network state estimate from previous solver iteration
  /**
   * See ::update() for why this is necessary.
   * 
   */
  boost::shared_ptr<math::RealVector> Xold;

  /**
   * The current network state estimate.
   * This vector provides a space for the nonlinear solve to store the
   * current solution estimate.  It should be filled with the initial
   * condition, handed to the nonlinear solver, and not changed
   * afterward.
   * 
   */
  boost::shared_ptr<math::RealVector> X;

  /**
   * The difference between the current and previous estimate.
   * See ::update() for why this is necessary.
   * 
   */
  boost::shared_ptr<math::RealVector> Xdelta;

  /** 
   * Constructor
   * The current network state is gathered from the network.
   * @param factory powerflow factory 
   * @param network network controlled by @c factory
   * @return 
   */
  PFSolverHelper(boost::shared_ptr<PFFactoryModule> factory,
      boost::shared_ptr<PFNetwork> network)
    : p_factory(factory), p_network(network), Xold(), Xdelta()
  {
    p_factory->setMode(State);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    Xold = vMap.mapToRealVector();
    // Xold->print();
    X.reset(Xold->clone());
    Xdelta.reset(Xold->clone());
    Xdelta->zero();
    p_factory->setMode(Jacobian);
    mapper::FullMatrixMap<PFNetwork> jMap(p_network);
    J = jMap.mapToRealMatrix();
  }
  
  /** 
   * Push the current estimated state back onto the network.
   * The network state (voltage, phase) is updated with the current
   * estimate from the nonlinear solver.
   * 
   * FIXME: The problem here is that ...::mapToBus expects the
   * @e CHANGE (old - current)? in state variables. IMHO, there should
   * be a way to set the state variables directly.  The solver should
   * be responsible for making the @e entire estimate.
   *
   * @param Xcur current state estimate from the solver
   */
  void
  update(const math::RealVector& Xcur)
  {
    
    Xdelta->equate(Xcur);
    Xdelta->scale(-1.0);
    Xdelta->add(*Xold);
    double snorm(Xdelta->norm2());
    if (Xdelta->processor_rank() == 0) {
      std::cout << "PFSolverHelper::update(): solution residual: " << snorm << std::endl;
    }

    // Xdelta->print();
    p_factory->setMode(RHS);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    vMap.mapToBus(Xdelta);
    Xold->equate(Xcur);
    
    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    p_network->updateBuses();
    
  }
  
  /** 
   * Build the Jacobian Matrix.
   * This is called by the nonlinear solver each iteration to build
   * the Jacobian from the current network state.  
   *
   * @param Xcur current state estimate
   * @param J Jacobian 
   */
  void
  operator() (const math::RealVector& Xcur, math::RealMatrix& theJ)
  {
    // In both the Netwon-Raphson and PETSc nonlinear solver (some
    // methods) implementations, the RHS function builder is called
    // before this, so we may be able to count on the current solution
    // being on the netork when here.

    // X.print();
    // update(Xcur);
    
    // Set to build Jacobian
    printf("Got to operator J\n");
    p_factory->setMode(Jacobian);
    mapper::FullMatrixMap<PFNetwork> jMap(p_network);
    
    // build the Jacobian
    jMap.mapToRealMatrix(theJ);
  }
  
  /** 
   * Build the RHS function vector.
   * This is called by the nonlinear solver each iteration to build
   * the RHS vector from the current network state.  This is also
   * responsible for updating the network state with the current
   * solution estimate, which assumes that it is called only once per
   * solver iteration.  This may not be true if certain methods are
   * used or if a finite difference Jacobian is computed by the
   * solver.
   * 
   * @param Xcur current state estimate
   * @param PQ computed RHS vector
   */
  void
  operator() (const math::RealVector& Xcur, math::RealVector& PQ)
  {
    // In both the Netwon-Raphson and PETSc nonlinear solver
    // implementations, this is called before the Jacobian builder, so
    // we may only need to map the solution back on to the network here.
    
    // X.print();
    update(Xcur);
    
    // set to build RHS vector
    p_factory->setMode(RHS);
    mapper::BusVectorMap<PFNetwork> vMap(p_network);
    
    // build the RHS vector
    vMap.mapToRealVector(PQ);
    printf("norm of PQ: %f\n",PQ.norm2());
  }
};
}
}
