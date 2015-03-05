// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   complex_linear_solver_test.cpp
 * @author William A. Perkins
 * @date   2015-03-05 10:06:10 d3g096
 * 
 * @brief  
 * 
 * 
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
#include "linear_solver.hpp"

/// The configuration used for these tests
static gridpack::utility::Configuration::CursorPtr test_config;


BOOST_AUTO_TEST_SUITE(ComplexLinearSolverTest)

// ------------------------------------------------------------- This
// test solves the Helmboltz equation on a unit square. It's taken
// directly from PETSc KSP example 11
// (http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex11.c.html)
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Helmboltz)
{
  using namespace gridpack;

  const int n(6), dim(n*n);

  parallel::Communicator world;
  int local_size(dim/world.size());

  if (world.rank() == world.size() - 1) {
    local_size = dim - local_size*(world.size() - 1);
  }

  boost::scoped_ptr<math::ComplexMatrix> 
    A(new math::ComplexMatrix(world, local_size, local_size, 12));
  std::cerr << world.rank() << " of " << world.size() 
            << ": local_size = " << local_size 
            << std::endl;

  const RealType h2(1.0/static_cast<RealType>((n+1)*(n+1)));
  const RealType sigma1(100.0);
  const ComplexType sigma2(0.0, 10.0);
  
  int lo, hi;
  A->localRowRange(lo, hi);
  
  for (int row = lo; row < hi; ++row) {
    ComplexType v(-1.0);
    int i(row/n), j(row - i*n);
    if (i > 0) A->addElement(row, row - n, v);
    if (i < n-1) A->addElement(row, row + n, v);
    if (j > 0) A->addElement(row, row - 1, v);
    if (j < n-1) A->addElement(row, row + 1, v);
    v = 4.0 - sigma1*h2 + sigma2*h2;
    A->addElement(row, row, v);
  }
  A->ready();

  boost::scoped_ptr<math::ComplexVector> 
    u(new math::ComplexVector(world, local_size)),
    x(u->clone());
  u->fill(0.5);
  u->ready();
  x->fill(0.0);
  x->ready();

  // rig the game, so we know the answer
  boost::scoped_ptr<math::ComplexVector>
    b(math::multiply(*A, *u));

  boost::scoped_ptr<math::LinearSolver> 
    solver(new math::LinearSolver(*A));
  BOOST_REQUIRE(test_config);
  solver->configure(test_config);
  solver->solve(*b, *x);

  // check the answer
  x->add(*u, -1.0);
  BOOST_CHECK(x->norm2() < 1.0E-04);
  
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
