// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   dae_solver_test.cpp
 * @author William A. Perkins
 * @date   2019-12-04 12:40:57 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/assert.hpp>
#include "dae_solver.hpp"

#include "test_main.cpp"

#ifdef TEST_REAL

typedef gridpack::RealType TestType;

#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE((y), (x), delta);

#define TEST_VALUE(r, i) (r)

#else 

typedef gridpack::ComplexType TestType;
#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE(real(y), real(x), delta); \
  BOOST_CHECK_CLOSE( abs(y), abs(x), delta);

#define TEST_VALUE(r, i) TestType(r,i)

#endif


// -------------------------------------------------------------
//  class Problem
// -------------------------------------------------------------
class Problem 
  : private gridpack::utility::Uncopyable
{
public:

  typedef gridpack::math::DAESolverT<TestType> TheSolverType;
  typedef TheSolverType::VectorType VectorType;
  typedef TheSolverType::MatrixType MatrixType;

  /// Default constructor.
  Problem(const int& local_size, const double& maxtime, const double outstep)
    : p_size(local_size),
      p_maxtime(maxtime),
      p_outstep(outstep),
      p_maxsteps(1000)
  {}

  /// Destructor
  virtual ~Problem(void) 
  {}

  /// Get the problem size
  int size(void) const
  {
    return p_size;
  }

  /// Build a Jacobian
  virtual void operator() (const double& time, 
                           const VectorType& X, 
                           const VectorType& Xdot, 
                           const double& shift, MatrixType& J) = 0;

  /// Build the RHS function
  virtual void operator() (const double& time, 
                           const VectorType& X, const VectorType& Xdot, 
                           VectorType& F) = 0;

  /// Report the time before the timestep
  static void reportPreTime(const double time)
  {
    std::cerr << "Pre  time step: the time is now " << time << std::endl;
  }

  /// Report the time after the timestep
  static void reportPostTime(const double time)
  {
    std::cerr << "Post time step: the time is now " << time << std::endl;
  }
  

  /// Get the initial solution
  virtual VectorType *initial(const gridpack::parallel::Communicator& comm) = 0;
  

  /// Solve the problem
  void solve(const gridpack::parallel::Communicator& comm,
             gridpack::utility::Configuration::CursorPtr conf)
  {
    TheSolverType::JacobianBuilder jbuilder = boost::ref(*this);
    TheSolverType::FunctionBuilder fbuilder = boost::ref(*this);
    TheSolverType::StepFunction sfunc;
    TheSolverType::EventManagerPtr eman;
    
    TheSolverType solver(comm, p_size, jbuilder, fbuilder, eman);
    solver.configure(conf);

    sfunc = &reportPreTime;
    solver.preStep(sfunc);
    sfunc = &reportPostTime;
    solver.postStep(sfunc);

    std::auto_ptr<VectorType> x(initial(comm));
    
    double t0(0.0), t(t0);
    solver.initialize(t0, 0.001, *x);

    for (int i = 1; t <= p_maxtime; ++i) {
      t = t0 + p_outstep*i;
      int mxstep(p_maxsteps);
      solver.solve(t, mxstep);
      std::cout << "Time = " << t << ", Steps = " << mxstep << std::endl;
      x->print();
    }
  }

protected:

  /// The problem size
  int p_size;

  /// The maximum/end time
  double p_maxtime;

  /// The output time step
  double p_outstep;

  /// The number of steps allowed/taken (per ::outstep)
  int p_maxsteps;

};

// -------------------------------------------------------------
// class RoberProblem
// -------------------------------------------------------------
class RoberProblem 
  : public Problem
{
public:

  /// Default constructor.
  RoberProblem(void) 
    : Problem(3, 1.0E+11, 1.0e+09)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   const double& shift, MatrixType& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size());
    X.getElementRange(lo, hi, &x[0]);
    J.setElement(lo+0, lo+0, shift + 0.04);
    J.setElement(lo+0, lo+1, -1.0e+04*x[2]);
    J.setElement(lo+0, lo+2, -1.0e+04*x[1]);
    J.setElement(lo+1, lo+0, -0.04);
    J.setElement(lo+1, lo+1, shift + 1.0e+04*x[2] + 3.0e+07*2.0*x[1]);
    J.setElement(lo+1, lo+2, 1.0e+04*x[1]);
    J.setElement(lo+2, lo+0, 0.0);
    J.setElement(lo+2, lo+1, -3.0e+07*2.0*x[1]);
    J.setElement(lo+2, lo+2, shift);
    J.ready();
  }

  /// Build the RHS vector
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   VectorType& F)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] + 0.04*x[0] - 1.0E+04*x[1]*x[2]);
    F.setElement(lo+1, xdot[1] - 0.04*x[0] + 1.0E+04*x[1]*x[2] + 3.0E+07*x[1]*x[1]);
    F.setElement(lo+2, xdot[2] - 3.0E+07*x[1]*x[1]);
    F.ready();
  }

  /// Get the initial solution
  /// 
  VectorType *initial(const gridpack::parallel::Communicator& comm)
  {
    VectorType *result = 
      new VectorType(comm, this->size());
    result->zero();
    int lo, hi;
    result->localIndexRange(lo, hi);
    result->setElement(lo+0, 1.0);
    result->ready();
    return result;
  }

};

// -------------------------------------------------------------
//  class OregoProblem
// -------------------------------------------------------------
class OregoProblem 
  : public Problem
{
public:

  /// Default constructor.
  OregoProblem()
    : Problem(3, 360.0, 30.0)
  {}

  /// Destructor
  ~OregoProblem(void)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   const double& shift, MatrixType& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size());
    X.getElementRange(lo, hi, &x[0]);
    J.setElement(lo+0, lo+0, shift - 77.27*((1. - 8.375e-6*x[0] - x[1]) - 8.375e-6*x[0]));
    J.setElement(lo+0, lo+1, -77.27*(1. - x[0]));
    J.setElement(lo+0, lo+2, 0.0);
    J.setElement(lo+1, lo+0, 1./77.27*x[1]);
    J.setElement(lo+1, lo+1, shift + 1./77.27*(1. + x[0]));
    J.setElement(lo+1, lo+2, -1./77.27);
    J.setElement(lo+2, lo+0, -0.161);
    J.setElement(lo+2, lo+1, 0.0);
    J.setElement(lo+2, lo+2, shift);
    J.ready();
  }

  /// Build the RHS vector
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   VectorType& F)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] - 77.27*(x[1] + x[0]*(1. - 8.375e-6*x[0] - x[1])));
    F.setElement(lo+1, xdot[1] - 1/77.27*(x[2] - (1. + x[0])*x[1]));
    F.setElement(lo+2, xdot[2] - 0.161*(x[0] - x[2]));
    F.ready();
  }

  /// Get the initial solution
  /// 
  VectorType *initial(const gridpack::parallel::Communicator& comm)
  {
    VectorType *result = 
      new VectorType(comm, this->size());
    result->zero();
    int lo, hi;
    result->localIndexRange(lo, hi);
    result->setElement(lo+0, 1.0);
    result->setElement(lo+1, 2.0);
    result->setElement(lo+2, 3.0);
    result->ready();
    return result;
  }

};

// -------------------------------------------------------------
//  class CEProblem
// -------------------------------------------------------------
class CEProblem 
  : public Problem
{
public:

  /// Default constructor.
  CEProblem()
    : Problem(1, 10.0, 0.5),
      p_lambda(10)
  {}

  /// Destructor
  ~CEProblem(void)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const VectorType& X, const VectorType& Xdot, 
                   const double& shift, MatrixType& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size());
    X.getElementRange(lo, hi, &x[0]);
    J.setElement(lo+0, lo+0, shift + p_lambda);
    J.ready();
  }

  /// Build the RHS vector
  void operator() (const double& t, 
                   const VectorType& X, const VectorType& Xdot, 
                   VectorType& F)
  {
    double l(p_lambda);
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<TestType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] + l*(x[0] - cos(t)));
    F.ready();
  }

  /// Get the initial solution
  /// 
  VectorType *initial(const gridpack::parallel::Communicator& comm)
  {
    VectorType *result = 
      new VectorType(comm, this->size());
    result->zero();
    int lo, hi;
    double t(0.0);
    double l(p_lambda);
    result->localIndexRange(lo, hi);
    result->setElement(lo+0, l/(l*l+1)*(l*cos(t)+sin(t)) - l*l/(l*l+1)*exp(-l*t));
    result->ready();
    return result;
  }

protected:

  double p_lambda;

};



BOOST_AUTO_TEST_SUITE(DAESolverTest)

BOOST_AUTO_TEST_CASE( Rober )
{
  gridpack::parallel::Communicator world;

  std::auto_ptr<Problem> p(new RoberProblem());

  p->solve(world, test_config);
}

BOOST_AUTO_TEST_CASE( Orego )
{
  gridpack::parallel::Communicator world;

  std::auto_ptr<Problem> p(new OregoProblem());

  p->solve(world, test_config);

}

BOOST_AUTO_TEST_CASE( CE )
{
  gridpack::parallel::Communicator world;

  std::auto_ptr<Problem> p(new CEProblem());

  p->solve(world, test_config);

}


BOOST_AUTO_TEST_SUITE_END()



