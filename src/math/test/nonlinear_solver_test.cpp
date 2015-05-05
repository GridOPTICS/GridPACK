// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_test.cpp
 * @author William A. Perkins
 * @date   2015-03-27 08:24:09 d3g096
 * 
 * @brief  Unit tests for NonlinearSolver
 * 
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/format.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/parallel/parallel.hpp"
#include "gridpack/utilities/exception.hpp"
#include "math.hpp"
#include "nonlinear_solver.hpp"
#include "newton_raphson_solver.hpp"

/// The configuration used for these tests
static gridpack::utility::Configuration::CursorPtr test_config;

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

typedef gridpack::math::NonlinearSolverT<TestType> TheNonlinearSolver;
typedef gridpack::math::NewtonRaphsonSolverT<TestType> TheNewtonRaphsonSolver;
typedef typename TheNonlinearSolver::VectorType VectorType;
typedef typename TheNonlinearSolver::MatrixType MatrixType;

BOOST_AUTO_TEST_SUITE(NonlinearSolverTest)


// -------------------------------------------------------------
// In this test, a very small nonlinear system is solved, too small to
// be parallel. So, each process solves it separately.
// -------------------------------------------------------------

struct build_tiny_jacobian_1
{
  void operator() (const VectorType& X, MatrixType& J) const
  {
    TestType x, y;
    X.getElement(0, x);
    X.getElement(1, y);
    J.setElement(0, 0, 2.0*x-2.0);
    J.setElement(0, 1, -1);
    J.setElement(1, 0, 2.0*x);
    J.setElement(1, 1, 8.0*y);
    J.ready();
    // J.print();
  }
};

struct build_tiny_function_1
{
  void operator() (const VectorType& X, VectorType& F) const
  {
    TestType x, y;
    X.getElement(0, x);
    X.getElement(1, y);
    F.setElement(0, x*x - 2.0*x - y + 0.5);
    F.setElement(1, x*x + 4.0*y*y - 4.0);
    F.ready();
    // F.print();
  }
};

BOOST_AUTO_TEST_CASE( tiny_serial_1 )
{
  boost::mpi::communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  TheNonlinearSolver::JacobianBuilder j = build_tiny_jacobian_1();
  TheNonlinearSolver::FunctionBuilder f = build_tiny_function_1();

  TheNonlinearSolver solver(self, 2, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  VectorType X(self, 2);
  X.setElement(0, 2.00);
  X.setElement(1, 0.25);
  X.ready();
  solver.solve(X);
  BOOST_TEST_MESSAGE("tiny_serial_1 results:");
  X.print();

  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);

  TEST_VALUE_CLOSE(x, static_cast<TestType>(1.900677), 1.0e-04);
  TEST_VALUE_CLOSE(y, static_cast<TestType>(0.3112186), 1.0e-04);

}

BOOST_AUTO_TEST_CASE( tiny_nr_serial_1 )
{
  boost::mpi::communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  TheNewtonRaphsonSolver::JacobianBuilder j = build_tiny_jacobian_1();
  TheNewtonRaphsonSolver::FunctionBuilder f = build_tiny_function_1();

  TheNewtonRaphsonSolver solver(self, 2, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  // check to see if correct values came from the configuration
  BOOST_CHECK_CLOSE(solver.tolerance(), 1.0e-10, 1.0e-04);
  BOOST_CHECK_EQUAL(solver.maximumIterations(), 100);

  VectorType X(self, 2);
  X.setElement(0, 2.00);
  X.setElement(1, 0.25);
  X.ready();
  solver.solve(X);
  BOOST_TEST_MESSAGE("tiny_serial_1 results (newton-raphson):");
  X.print();

  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);

  TEST_VALUE_CLOSE(x, static_cast<TestType>(1.900677), 1.0e-04);
  TEST_VALUE_CLOSE(y, static_cast<TestType>(0.3112186), 1.0e-04);

}

// -------------------------------------------------------------
// Another tiny test. This one is example 1 (not hard) from the PETSc
// SNES examples.
// -------------------------------------------------------------

void 
build_tiny_jacobian_2(const VectorType& X, MatrixType& J)
{
  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);
  J.setElement(0, 0, 2.0*x + y);
  J.setElement(0, 1, x);
  J.setElement(1, 0, y);
  J.setElement(1, 1, x + 2.0*y);
  J.ready();
  // J.print();
}

void
build_tiny_function_2(const VectorType& X, VectorType& F)
{
  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);
  F.setElement(0, x*x + x*y - 3.0);
  F.setElement(1, x*y + y*y - 6.0);
  F.ready();
  // F.print();
}

BOOST_AUTO_TEST_CASE( tiny_serial_2 )
{
  boost::mpi::communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  TheNonlinearSolver::JacobianBuilder j = &build_tiny_jacobian_2;
  TheNonlinearSolver::FunctionBuilder f = &build_tiny_function_2;

  TheNonlinearSolver solver(self, 2, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  VectorType X(self, 2);
  X.setElement(0, 2.00);
  X.setElement(1, 3.00);
  X.ready();
  solver.solve(X);

  BOOST_TEST_MESSAGE("tiny_serial_2 results:");
  X.print();

  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);

  TEST_VALUE_CLOSE(x, static_cast<TestType>(1.0), 1.0e-04);
  TEST_VALUE_CLOSE(y, static_cast<TestType>(2.0), 1.0e-04);
}

BOOST_AUTO_TEST_CASE( tiny_nr_serial_2 )
{
  boost::mpi::communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  TheNewtonRaphsonSolver::JacobianBuilder j = &build_tiny_jacobian_2;
  TheNewtonRaphsonSolver::FunctionBuilder f = &build_tiny_function_2;

  TheNewtonRaphsonSolver solver(self, 2, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  // check to see if correct values came from the configuration
  BOOST_CHECK_CLOSE(solver.tolerance(), 1.0e-10, 1.0e-04);
  BOOST_CHECK_EQUAL(solver.maximumIterations(), 100);

  VectorType X(self, 2);
  X.setElement(0, 2.00);
  X.setElement(1, 3.00);
  X.ready();
  solver.solve(X);

  BOOST_TEST_MESSAGE("tiny_serial_2 results:");
  X.print();

  TestType x, y;
  X.getElement(0, x);
  X.getElement(1, y);

  TEST_VALUE_CLOSE(x, static_cast<TestType>(1.0), 1.0e-04);
  TEST_VALUE_CLOSE(y, static_cast<TestType>(2.0), 1.0e-04);
}

// -------------------------------------------------------------
// A larger test.  This is example 2 from the PETSc SNES examples
// -------------------------------------------------------------
struct build_thing
{
  void operator() (const VectorType& X, MatrixType& J) const
  {
    int n(X.size());
    TestType d(static_cast<double>(n - 1));
    TestType h(1.0/d);
    d *= d;

    int lo, hi;
    X.localIndexRange(lo, hi);

    for (int row = lo; row < hi; ++row) {
      if (row == 0 || row == n - 1) {
        J.setElement(row, row, 1.0);
      } else {
        int i[3] = { row,     row, row };
        int j[3] = { row - 1, row, row + 1 };
        TestType x;
        X.getElement(row, x);
        TestType A[3] = 
          { d, -2.0*d + 2.0*x, d};
        J.setElements(3, i, j, A);
      }
    }
    J.ready();
    std::cout << "build_jacobian_2 called" << std::endl;
    // J.print();
  }

  void operator() (const VectorType& X, VectorType& F) const
  {
    int n(X.size());
    TestType d(static_cast<double>(n - 1));
    TestType h(1.0/d);

    d *= d;

    int lo, hi;
    X.localIndexRange(lo, hi);
    
    std::vector<TestType> x(n);
    X.getAllElements(&x[0]);

    for (int row = lo; row < hi; ++row) {
      TestType f;
      int i(row);
      if (row == 0) {
        f = x[i];
      } else if (row == n - 1) {
        f = x[i] - 1.0;
      } else {
        TestType xp, g;
        xp = static_cast<double>(row)*h;
        g = 6.0*xp + pow(xp+1.0e-12, 6);
        f = d*(x[i-1] - 2.0*x[i] + x[i+1]) + x[i]*x[i] - g;
      }
      F.setElement(row, f);
    }
    F.ready();
    // F.print();
    std::cout << "build_function_2 called" << std::endl;
  }
};

BOOST_AUTO_TEST_CASE( example2 )
{
  boost::mpi::communicator world;
  int local_size(4);

  // Make sure local ownership specifications work
  if (world.size() > 1) {
    if (world.rank() == 0) {
      local_size -= 1;
    } else if (world.rank() == world.size() - 1) {
      local_size += 1;
    }
  }

  build_thing thing;

  TheNonlinearSolver::JacobianBuilder j = thing;
  TheNonlinearSolver::FunctionBuilder f = thing;

  TheNonlinearSolver solver(world, local_size, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  VectorType X(world, local_size);
  X.fill(0.5);
  X.ready();
  solver.solve(X);
  X.print();
}

BOOST_AUTO_TEST_CASE( example2_nr )
{
  boost::mpi::communicator world;
  int local_size(4);

  // Make sure local ownership specifications work
  if (world.size() > 1) {
    if (world.rank() == 0) {
      local_size -= 1;
    } else if (world.rank() == world.size() - 1) {
      local_size += 1;
    }
  }

  build_thing thing;

  TheNewtonRaphsonSolver::JacobianBuilder j = thing;
  TheNewtonRaphsonSolver::FunctionBuilder f = thing;

  TheNewtonRaphsonSolver solver(world, local_size, j, f);

  BOOST_REQUIRE(test_config);
  solver.configure(test_config);

  // check to see if correct values came from the configuration
  BOOST_CHECK_CLOSE(solver.tolerance(), 1.0e-10, 1.0e-04);
  BOOST_CHECK_EQUAL(solver.maximumIterations(), 100);

  VectorType X(world, local_size);
  X.fill(0.5);
  X.ready();
  solver.solve(X);
  X.print();
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
  gridpack::parallel::Communicator world;

  boost::scoped_ptr<gridpack::utility::Configuration>
    config(gridpack::utility::Configuration::configuration());
  
  config->enableLogging();
  config->open("gridpack.xml", world);

  test_config = config->getCursor("GridPACK.MathTests");

  gridpack::math::Initialize();

  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
  return result;
}
