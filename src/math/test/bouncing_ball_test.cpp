// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   bouncing_ball_test.cpp
 * @author Perkins
 * @date   2019-11-26 10:00:45 d3g096
 * 
 * @brief  DAESolver test based on PETSc TS example 40. 
 * 
 * 
 */
// -------------------------------------------------------------
// Created November 21, 2019 by Perkins
// Last Change: 2018-07-24 09:27:04 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <gridpack/math/dae_solver.hpp>

namespace gp = gridpack;
namespace gpu = gp::utility;
namespace gpp = gp::parallel;
namespace gpm = gp::math;

// -------------------------------------------------------------
//  class BouncingBallProblem
// -------------------------------------------------------------
template <typename TestType>
class BouncingBallProblem
  : public gpp::Distributed,
    private gpu::Uncopyable
{
public:

  typedef gpm::DAESolverT<TestType> SolverType;
  typedef typename SolverType::VectorType VectorType;
  typedef typename SolverType::MatrixType MatrixType;

  class BounceEvent;            // forward

  /// Default constructor.
  BouncingBallProblem(const gpp::Communicator& comm)
    : gpp::Distributed(comm),
      gpu::Uncopyable()
  {}

  /// Destructor
  ~BouncingBallProblem(void)
  {}

  /// Get the DAE event description
  typename SolverType::EventPtr events(void) const
  {
    typename SolverType::EventPtr result(new BounceEvent());
    return result;
  }
  
  /// Build a Jacobian
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   const double& shift, MatrixType& J)
  {
    
  }
  
  /// Build the RHS function
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   VectorType& F)
  {
  }
  
  /// solve the problem
  void solve(gpu::Configuration::CursorPtr conf)
  {
    typename SolverType::JacobianBuilder jbuilder = boost::ref(*this);
    typename SolverType::FunctionBuilder fbuilder = boost::ref(*this);
    typename SolverType::EventManagerPtr eman;
    typename SolverType::EventPtr bevent(new BounceEvent());
    SolverType solver(this->communicator(), p_size, jbuilder, fbuilder, eman);
    solver.configure(conf);
    eman->add(bevent);
  }

protected:

  /// The (local) problem size
  static const int p_size;

  /// The time limit
  double p_tmax;

  /// The reporting time step
  double p_tstep;

  /// Initial height of ball
  double p_h0;

public:
  
  // -------------------------------------------------------------
  //  class BounceEvent
  // -------------------------------------------------------------
  class BounceEvent
    : public SolverType::Event
  {
  public:

    /// Default constructor.
    BounceEvent(void)
      : SolverType::Event(2),
        p_maxBounce(5), p_numBounce(0)
    {
      this->p_dir[0] = gpm::CrossZeroNegative;
      this->p_dir[1] = gpm::CrossZeroNegative;
      this->p_term[0] = false;
      this->p_term[1] = true;
      this->p_stateIndex = 0;
    }

    /// Destructor
    ~BounceEvent(void)
    {}

  protected:

    /// Maximum number of bounces
    const int p_maxBounce;

    /// Number of bouces
    int p_numBounce;

    /// update and return event values, given the state (specialized)
    void p_update(const double& t, TestType *state)
    {
      this->p_current[0] = state[0];
      this->p_current[1] = p_maxBounce - p_numBounce;
    }
    
    /// handle any triggered events  (specialized)
    void p_handle(const bool *triggered, const double& t, TestType *state)
    {
      if (triggered[0]) {
        state[0] = 0.0;
        state[1] = -0.9*state[1];
        p_numBounce++;
      }
    }
  };
};

template <>
const int BouncingBallProblem<double>::p_size(2);

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gpp::Environment env(argc, argv);
  gpp::Communicator world;

  if (world.size() > 1) {
    std::cerr << "This program is serial" << std::endl;
    exit(3);
  }
  
  gpp::Communicator self = world.self();
  boost::scoped_ptr<gpu::Configuration> 
    config(gpu::Configuration::configuration());
  gpu::Configuration::CursorPtr cursor =
    config->getCursor("BouncingBall");
  config->open("bouncing_ball_test.xml", self);
  BouncingBallProblem<double> problem(self);
  problem.solve(cursor);
  
  return 0;
}
