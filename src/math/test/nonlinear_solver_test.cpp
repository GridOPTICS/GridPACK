/**
 * @file   nonlinear_solver_test.cpp
 * @author William A. Perkins
 * @date   2013-08-12 14:21:01 d3g096
 * 
 * @brief  Unit tests for NonlinearSolver
 * 
 * 
 */

#include <mpi.h>
#include <iostream>
#include <boost/format.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/parallel/parallel.hpp"
#include "gridpack/utilities/exception.hpp"
#include "math.hpp"
#include "nonlinear_solver.hpp"



BOOST_AUTO_TEST_SUITE(NonlinearSolver)


// -------------------------------------------------------------
// In this test, a very small nonlinear system is solved, too small to
// be parallel. So, each process solves it separately.
// -------------------------------------------------------------

struct build_tiny_jacobian_1
{
  void operator() (const gridpack::math::Vector& X, gridpack::math::Matrix& J) const
  {
    gridpack::ComplexType x, y;
    X.get_element(0, x);
    X.get_element(1, y);
    J.set_element(0, 0, 2.0*x-2.0);
    J.set_element(0, 1, -1);
    J.set_element(1, 0, 2.0*x);
    J.set_element(1, 1, 8.0*y);
    J.ready();
    J.print();
  }
};

struct build_tiny_function_1
{
  void operator() (const gridpack::math::Vector& X, gridpack::math::Vector& F) const
  {
    gridpack::ComplexType x, y;
    X.get_element(0, x);
    X.get_element(1, y);
    F.set_element(0, x*x - 2.0*x - y + 0.5);
    F.set_element(1, x*x + 4.0*y*y - 4.0);
    F.ready();
    F.print();
  }
};

BOOST_AUTO_TEST_CASE( tiny_serial_1 )
{
  boost::mpi::communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  gridpack::math::JacobianBuilder j = build_tiny_jacobian_1();
  gridpack::math::FunctionBuilder f = build_tiny_function_1();

  gridpack::math::NonlinearSolver solver(self, 2, j, f);
  gridpack::math::Vector X(self, 2);
  X.set_element(0, 2.00);
  X.set_element(1, 0.25);
  X.ready();
  solver.solve(X);
  X.print();

  gridpack::ComplexType x, y;
  X.get_element(0, x);
  X.get_element(1, y);

  BOOST_CHECK_CLOSE(real(x), 1.900677, 1.0e-04);
  BOOST_CHECK_CLOSE(real(y), 0.3112186, 1.0e-04);
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
  gridpack::math::Initialize();
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
}
