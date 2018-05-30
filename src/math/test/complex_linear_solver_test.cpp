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
 * @date   2016-12-16 09:38:23 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/format.hpp>
#include "linear_solver.hpp"
#include "linear_matrix_solver.hpp"

#include "test_main.cpp"

using namespace gridpack;

static void
assemble_helmboltz(const parallel::Communicator& comm,
                   boost::shared_ptr<math::ComplexMatrix>& A, 
                   boost::shared_ptr<math::ComplexVector>& u, 
                   boost::shared_ptr<math::ComplexVector>& b)
{
  const int n(6), dim(n*n);

  int local_size(dim/comm.size());

  if (comm.rank() == comm.size() - 1) {
    local_size = dim - local_size*(comm.size() - 1);
  }
  std::cerr << comm.rank() << " of " << comm.size() 
            << ": local_size = " << local_size 
            << std::endl;

  A.reset(new math::ComplexMatrix(comm, local_size, local_size, 6));

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

  u.reset(new math::ComplexVector(comm, local_size));
  u->fill(0.5);
  u->ready();

  b.reset(u->clone());
  math::multiply(*A, *u, *b);
}


BOOST_AUTO_TEST_SUITE(ComplexLinearSolverTest)

// ------------------------------------------------------------- This
// test solves the Helmboltz equation on a unit square. It's taken
// directly from PETSc KSP example 11
// (http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex11.c.html)
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Helmboltz)
{
  parallel::Communicator world;
  boost::shared_ptr<math::ComplexMatrix> A;
  boost::shared_ptr<math::ComplexVector> u, b, x;
  assemble_helmboltz(world, A, u, b);

  x.reset(u->clone()),
  x->fill(0.0);
  x->ready();

  boost::scoped_ptr<math::LinearSolver> 
    solver(new math::LinearSolver(*A));
  BOOST_REQUIRE(test_config);
  solver->configure(test_config);
  solver->solve(*b, *x);

  // check the answer
  x->add(*u, -1.0);
  BOOST_CHECK(x->norm2() < 1.0E-04);
}

BOOST_AUTO_TEST_CASE(HelmboltzInverse)
{
  parallel::Communicator world;
  boost::shared_ptr<math::ComplexMatrix> A, I, Ainv;
  boost::shared_ptr<math::ComplexVector> u, b, x;
  assemble_helmboltz(world, A, u, b);

  I.reset(new math::ComplexMatrix(A->communicator(),
                                  A->localRows(),
                                  A->localCols(),
                                  math::Dense));
  I->identity();

  boost::scoped_ptr<math::LinearMatrixSolver> 
    solver(new math::LinearMatrixSolver(*A));

  BOOST_REQUIRE(test_config);
  solver->configure(test_config);
  Ainv.reset(solver->solve(*I));

  x.reset(math::multiply(*Ainv, *b));

  // check the answer
  x->add(*u, -1.0);
  BOOST_CHECK(x->norm2() < 1.0E-04);
}

BOOST_AUTO_TEST_SUITE_END()

