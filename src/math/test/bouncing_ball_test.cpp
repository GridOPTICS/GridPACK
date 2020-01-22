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
 * @date   2020-01-22 09:00:32 d3g096
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
#include <gridpack/math/math.hpp>
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
      gpu::Uncopyable(),
      p_tmax(30.0),
      p_tstep(0.1),
      p_h0(0.0),
      p_v0(20.0)
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
    int lo, hi;
    X.localIndexRange(lo, hi);
    J.setElement(lo+0, lo+0, shift);
    J.setElement(lo+0, lo+1, -1.0);
    J.setElement(lo+1, lo+0, 0.0);
    J.setElement(lo+1, lo+1, shift);
    J.ready();
    
  }
  
  /// Build the RHS function
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   VectorType& F)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    std::vector<TestType> x(this->p_size), xdot(this->p_size);
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] - x[1]);
    F.setElement(lo+1, xdot[1] + 9.8);
    F.ready();
  }
  
  /// solve the problem
  void solve(gpu::Configuration::CursorPtr conf)
  {
    typename SolverType::JacobianBuilder jbuilder = boost::ref(*this);
    typename SolverType::FunctionBuilder fbuilder = boost::ref(*this);
    typename SolverType::EventManagerPtr eman(new typename SolverType::EventManager());
    eman->add(this->events());

    SolverType solver(this->communicator(), p_size, jbuilder, fbuilder, eman);
    solver.configure(conf);

    boost::scoped_ptr<VectorType> x(new VectorType(this->communicator(), this->p_size));
    int lo, hi;
    x->localIndexRange(lo, hi);

    x->setElement(lo+0, p_h0);
    x->setElement(lo+1, p_v0);
    x->ready();
    x->print();

    double t0(0.0), t(t0);
    solver.initialize(t0, p_tstep, *x);

    
    for (int i = 1; t <= (p_tmax - 0.5*p_tstep); ++i) {
      t = t0 + p_tstep*i;
      int mxstep(10000);
      double tout(t);
      solver.solve(tout, mxstep);
      // std::cout << "Time requested = " << t << ", "
      //           << "actual time = " << tout << ", "
      //           << "Steps = " << mxstep << std::endl;
      // x->print();
      if (solver.terminated()) break;
    }
    
    
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

  /// Initiall ball velocity
  double p_v0;

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
      // std::cout << "p_update: t = " << t
      //           << ", state[0] = " << state[0]
      //           << ", state[1] = " << state[1]
      //           << std::endl;
      this->p_current[0] = state[0];
      this->p_current[1] = p_maxBounce - p_numBounce;
    }
    
    /// handle any triggered events  (specialized)
    void p_handle(const bool *triggered, const double& t, TestType *state)
    {
      if (triggered[0]) {
        p_numBounce++;
        // std::cout << "p_handle, bounce #" << p_numBounce << ": "
        //           << "t = " << t << ", "
        //           << "state[0] = " << state[0] << ", "
        //           << "state[1] = " << state[1]
        //           << std::endl;
        state[0] = 0.0;
        state[1] = -0.9*state[1];
        std::cout << "Bounce #" << p_numBounce
                  << " at t = " << t << std::endl;
      }
      if (triggered[1]) {
        std::cout << "Maximum number of bounces reached at t = "
                  << t << std::endl;
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
  gp::Environment env(argc, argv);
  gpp::Communicator world;

  if (world.size() > 1) {
    std::cerr << "This program is serial" << std::endl;
    exit(3);
  }
  gridpack::math::Initialize(&argc,&argv);
  
  gpp::Communicator self = world.self();
  boost::scoped_ptr<gpu::Configuration> 
    config(gpu::Configuration::configuration());
  config->open("bouncing_ball_test.xml", self);
  gpu::Configuration::CursorPtr cursor =
    config->getCursor("BouncingBall");
  BouncingBallProblem<double> problem(self);
  problem.solve(cursor);

  gridpack::math::Finalize();
  
  return 0;
}
