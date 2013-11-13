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
 * @date   2013-11-13 14:33:21 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/assert.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/parallel/parallel.hpp"
#include "gridpack/utilities/exception.hpp"
#include "math.hpp"
#include "dae_solver.hpp"

/// The configuration used for these tests
static gridpack::utility::Configuration::CursorPtr test_config;

// -------------------------------------------------------------
//  class Problem
// -------------------------------------------------------------
class Problem 
  : private gridpack::utility::Uncopyable
{
public:

  /// Default constructor.
  Problem(const int& local_size, const double& maxtime)
    : p_size(local_size),
      p_maxtime(maxtime),
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
                           const gridpack::math::Vector& X, 
                           const gridpack::math::Vector& Xdot, 
                           const double& shift, gridpack::math::Matrix& J) = 0;

  /// Build the RHS function
  virtual void operator() (const double& time, 
                           const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                           gridpack::math::Vector& F) = 0;

  /// Get the initial solution
  virtual gridpack::math::Vector *initial(const gridpack::parallel::Communicator& comm) = 0;
  

  /// Solve the problem
  void solve(const gridpack::parallel::Communicator& comm,
             gridpack::utility::Configuration::CursorPtr conf)
  {
    gridpack::math::DAEJacobianBuilder jbuilder = boost::ref(*this);
    gridpack::math::DAEFunctionBuilder fbuilder = boost::ref(*this);
    
    gridpack::math::DAESolver solver(comm, p_size, jbuilder, fbuilder);
    solver.configure(conf);

    std::auto_ptr<gridpack::math::Vector> x(initial(comm));

    solver.solve(0.0, 0.001, p_maxtime, p_maxsteps, *x);
  }

protected:

  /// The problem size
  int p_size;

  /// The maximum/end time
  double p_maxtime;

  /// The number of steps allowed/taken
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
    : Problem(3, 1.0E+11)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   const double& shift, gridpack::math::Matrix& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size());
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
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   gridpack::math::Vector& F)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] + 0.04*x[0] - 1.0E+04*x[1]*x[2]);
    F.setElement(lo+1, xdot[1] - 0.04*x[0] + 1.0E+04*x[1]*x[2] + 3.0E+07*x[1]*x[1]);
    F.setElement(lo+2, xdot[2] - 3.0E+07*x[1]*x[1]);
    F.ready();
  }

  /// Get the initial solution
  /// 
  gridpack::math::Vector *initial(const gridpack::parallel::Communicator& comm)
  {
    gridpack::math::Vector *result = 
      new gridpack::math::Vector(comm, this->size());
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
    : Problem(3, 360.0)
  {}

  /// Destructor
  ~OregoProblem(void)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   const double& shift, gridpack::math::Matrix& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size());
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
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   gridpack::math::Vector& F)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] - 77.27*(x[1] + x[0]*(1. - 8.375e-6*x[0] - x[1])));
    F.setElement(lo+1, xdot[1] - 1/77.27*(x[2] - (1. + x[0])*x[1]));
    F.setElement(lo+2, xdot[2] - 0.161*(x[0] - x[2]));
    F.ready();
  }

  /// Get the initial solution
  /// 
  gridpack::math::Vector *initial(const gridpack::parallel::Communicator& comm)
  {
    gridpack::math::Vector *result = 
      new gridpack::math::Vector(comm, this->size());
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
    : Problem(1, 10.0),
      p_lambda(10)
  {}

  /// Destructor
  ~CEProblem(void)
  {}

  /// Build a Jacobian
  void operator() (const double& time, 
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   const double& shift, gridpack::math::Matrix& J)
  {
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size());
    X.getElementRange(lo, hi, &x[0]);
    J.setElement(lo+0, lo+0, shift + p_lambda);
    J.ready();
  }

  /// Build the RHS vector
  void operator() (const double& t, 
                   const gridpack::math::Vector& X, const gridpack::math::Vector& Xdot, 
                   gridpack::math::Vector& F)
  {
    double l(p_lambda);
    int lo, hi;
    X.localIndexRange(lo, hi);
    BOOST_ASSERT((hi-lo) == this->size());
    std::vector<gridpack::ComplexType> x(this->size()), xdot(this->size());
    X.getElementRange(lo, hi, &x[0]);
    Xdot.getElementRange(lo, hi, &xdot[0]);
    F.setElement(lo+0, xdot[0] + l*(x[0] - cos(t)));
    F.ready();
  }

  /// Get the initial solution
  /// 
  gridpack::math::Vector *initial(const gridpack::parallel::Communicator& comm)
  {
    gridpack::math::Vector *result = 
      new gridpack::math::Vector(comm, this->size());
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

  // p->solve(world, test_config);
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



// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);

  boost::mpi::communicator world;
  boost::scoped_ptr<gridpack::utility::Configuration>
    config(gridpack::utility::Configuration::configuration());
  
  config->enableLogging();
  config->open("gridpack.xml", world);

  test_config = config->getCursor("GridPACK.MathTests");

  gridpack::math::Initialize();

  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
}
