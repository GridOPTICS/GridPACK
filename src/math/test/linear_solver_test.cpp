// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver_test.cpp
 * @author William A. Perkins
 * @date   2016-12-16 09:30:03 d3g096
 * 
 * @brief  
 * 
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/format.hpp>
#include "linear_solver.hpp"
#include "linear_matrix_solver.hpp"

#include "test_main.cpp"

BOOST_AUTO_TEST_SUITE(LinearSolverTest)

// -------------------------------------------------------------
// assemble
// -------------------------------------------------------------
void
assemble(const int imax, const int jmax, 
         gridpack::math::RealMatrix& A, 
         gridpack::math::RealVector& b)
{
  static const float k = 1000;  /* conductivity, W/m/K */
  static const float t = 0.01;  /* plate thickness, m */
  static const float W = 0.3;   /* plate width, m */
  static const float H = 0.4;   /* plate height, m */

  const float dx = W/static_cast<float>(imax);
  const float dy = H/static_cast<float>(jmax);

  int ilo, ihi;
  b.localIndexRange(ilo, ihi);

  int i, j;
  float ap, aw, ae, as, an, bp;

  int iP, iN, iS, iE, iW;

  // this fills a row in A only if it's locally owned

  for (i = 0; i < imax; i++) {
    for (j = 0; j < jmax; j++) {
      iP = i*jmax + j;
      if (ilo <= iP && iP < ihi) {
        iE = (i+1)*jmax + j;
        iW = (i-1)*jmax + j;
        iN = i*jmax + (j+1);
        iS = i*jmax + (j-1);

        bp = 0.0;
        ap = 0.0;
        if (j == 0) {             /* insulated south boundary */
          as = 0.0;
          bp += 0.0;
          ap -= 0.0;
        } else {
          as = (k/dy)*(dx*t);
        }

        if (j == jmax - 1) {      /* constant temperature (100C) north boundary */
          an = 0.0;
          bp += 2*k/dy*(dy*t)*100.0;
          ap -= -2*k/dy*(dy*t);
        } else {
          an = (k/dy)*(dx*t);
        }

        if (i == 0) {             /* constant flux (500kw/m2) west boundary */
          aw = 0.0;
          bp += 500000.0*dy*t;
          ap -= 0.0;
        } else {
          aw = (k/dx)*(dy*t);
        }

        if (i == imax - 1) {      /* insulated east boundary */
          ae = 0.0;
          bp += 0.0; 
          ap -= 0.0;
        } else {
          ae = (k/dx)*(dy*t);
        }
        
        ap += as + an + aw + ae;

        A.setElement(iP, iP, ap);

        if (an != 0.0) A.setElement(iP, iN, -an);
        if (as != 0.0) A.setElement(iP, iS, -as);
        if (ae != 0.0) A.setElement(iP, iE, -ae);
        if (aw != 0.0) A.setElement(iP, iW, -aw);
        b.setElement(iP, bp);
      }      
    }
  }
}


// -------------------------------------------------------------
/**
 * This is a simple test of the gridpack::math::LinearSolver.  This
 * problem comes from Example 7.2 in
 * 
 * Versteeg, H.K. and W. Malalasekera, 1995. An introduction to
 * computational fluid dynamics, the finite volume method. Prentice
 * Hall.  257 pp.
 * 
 */
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( Versteeg )
{
  gridpack::parallel::Communicator world;

  static const int imax = 3*world.size();
  static const int jmax = 4*world.size();
  static const int global_size = imax*jmax;
  int local_size(global_size/world.size());

  // Make sure local ownership specifications work
  if (world.size() > 1) {
    if (world.rank() == 0) {
      local_size -= 1;
    } else if (world.rank() == world.size() - 1) {
      local_size += 1;
    }
  }
    

  std::auto_ptr<gridpack::math::RealMatrix> 
    A(new gridpack::math::RealMatrix(world, local_size, local_size, 
                                 gridpack::math::Sparse));
  std::auto_ptr<gridpack::math::RealVector>
    b(new gridpack::math::RealVector(world, local_size)),
    x(new gridpack::math::RealVector(world, local_size));

  assemble(imax, jmax, *A, *b);
  A->ready();
  b->ready();

  x->fill(0.0);
  x->ready();

  A->print();
  b->print();

  std::auto_ptr<gridpack::math::RealLinearSolver> 
    solver(new gridpack::math::RealLinearSolver(*A));

  BOOST_REQUIRE(test_config);
  solver->configure(test_config);
  solver->solve(*b, *x);

  std::auto_ptr<gridpack::math::RealVector>
    res(multiply(*A, *x));
  res->add(*b, -1.0);

  double l1norm(res->norm1());
  double l2norm(res->norm2());

  if (world.rank() == 0) {
    std::cout << "Residual L1 Norm = " << l1norm << std::endl;
    std::cout << "Residual L2 Norm = " << l2norm << std::endl;
  }

  BOOST_CHECK(l1norm < 1.0e-05);
  BOOST_CHECK(l2norm < 1.0e-05);

  // solve it again

  x->fill(0.0);
  solver->resolve(*b, *x);
  multiply(*A, *x, *res);
  res->add(*b, -1.0);

  if (world.rank() == 0) {
    std::cout << "Residual L1 Norm = " << l1norm << std::endl;
    std::cout << "Residual L2 Norm = " << l2norm << std::endl;
  }

  BOOST_CHECK(l1norm < 1.0e-05);
  BOOST_CHECK(l2norm < 1.0e-05);

  l1norm = res->norm1();
  l2norm = res->norm2();

  for (int p = 0; p < world.size(); ++p) {
    if (p == world.rank()) {
      int ilo, ihi;
      x->localIndexRange(ilo, ihi);
      for (int iP = ilo; iP < ihi; ++iP) {
        gridpack::RealType val, r;
        x->getElement(iP, val);
        res->getElement(iP, r);
        int i = iP/jmax;
        int j = iP - i*jmax;
        
        std::cout << boost::str(boost::format("%8d%8d%8d%12.6f%12.3e") %
                                iP % i % j % val % r)
                  << std::endl;
        std::cout.flush();
      }
    }
    world.barrier();
  }
}

// FIXME
BOOST_AUTO_TEST_CASE ( VersteegInverse )
{
  gridpack::parallel::Communicator world;

  static const int imax = 3*world.size();
  static const int jmax = 4*world.size();
  static const int global_size = imax*jmax;
  int local_size(global_size/world.size());

  // Make sure local ownership specifications work
  if (world.size() > 1) {
    if (world.rank() == 0) {
      local_size -= 1;
    } else if (world.rank() == world.size() - 1) {
      local_size += 1;
    }
  }
    

  std::auto_ptr<gridpack::math::RealMatrix> 
    A(new gridpack::math::RealMatrix(world, local_size, local_size, 
                                 gridpack::math::Sparse)),
    I(new gridpack::math::RealMatrix(world, local_size, local_size, 
                                 gridpack::math::Sparse));
  I->identity();

  std::auto_ptr<gridpack::math::RealVector>
    b(new gridpack::math::RealVector(world, local_size));

  assemble(imax, jmax, *A, *b);
  A->ready();
  b->ready();

  std::auto_ptr<gridpack::math::RealLinearSolver> 
    solver(new gridpack::math::RealLinearSolver(*A));

  BOOST_REQUIRE(test_config);
  solver->configurationKey("LinearMatrixSolver");
  solver->configure(test_config);

  std::auto_ptr<gridpack::math::RealMatrix> 
    Ainv(solver->solve(*I));
  std::auto_ptr<gridpack::math::RealVector>
    x(multiply(*Ainv, *b));

  std::auto_ptr<gridpack::math::RealVector>
    res(multiply(*A, *x));
  res->add(*b, -1.0);
  res->print();

  // Ainv->print();

  double l1norm(res->norm1());
  double l2norm(res->norm2());

  if (world.rank() == 0) {
    std::cout << "Residual L1 Norm = " << l1norm << std::endl;
    std::cout << "Residual L2 Norm = " << l2norm << std::endl;
  }

  BOOST_CHECK(l1norm < 1.0e-05);
  BOOST_CHECK(l2norm < 1.0e-05);
  
}


// -------------------------------------------------------------
/// Test matrix inversion with LinearMatrixSolver
/**
 * Just like VersteegInverse, this solves the Versteeg heat transfer
 * problem by inverting the Matrix with LinearMatrixSolver.  This
 * would be a stupid way to solve the problem in real life.
 * 
 */
// -------------------------------------------------------------
// FIXME
BOOST_AUTO_TEST_CASE ( VersteegMatrixInverse )
{
  gridpack::parallel::Communicator world;

  static const int imax = 3*world.size();
  static const int jmax = 4*world.size();
  static const int global_size = imax*jmax;
  int local_size(global_size/world.size());

  // Make sure local ownership specifications work
  if (world.size() > 1) {
    if (world.rank() == 0) {
      local_size -= 1;
    } else if (world.rank() == world.size() - 1) {
      local_size += 1;
    }
  }
    

  boost::scoped_ptr<gridpack::math::RealMatrix> 
    A(new gridpack::math::RealMatrix(world, local_size, local_size, 
                                 gridpack::math::Sparse)),
    I(new gridpack::math::RealMatrix(world, local_size, local_size, 
                                 gridpack::math::Dense));
  I->identity();

  boost::scoped_ptr<gridpack::math::RealVector>
    b(new gridpack::math::RealVector(world, local_size));

  assemble(imax, jmax, *A, *b);
  A->ready();
  b->ready();

  boost::scoped_ptr<gridpack::math::RealLinearMatrixSolver> 
    solver(new gridpack::math::RealLinearMatrixSolver(*A));

  BOOST_REQUIRE(test_config);
  solver->configurationKey("LinearMatrixSolver");
  solver->configure(test_config);

  std::auto_ptr<gridpack::math::RealMatrix> 
    Ainv(solver->solve(*I));
  std::auto_ptr<gridpack::math::RealVector>
    x(multiply(*Ainv, *b));

  std::auto_ptr<gridpack::math::RealVector>
    res(multiply(*A, *x));
  res->add(*b, -1.0);

  // Ainv->print();

  double l1norm(res->norm1());
  double l2norm(res->norm2());

  if (world.rank() == 0) {
    std::cout << "Residual L1 Norm = " << l1norm << std::endl;
    std::cout << "Residual L2 Norm = " << l2norm << std::endl;
  }

  BOOST_CHECK(l1norm < 1.0e-05);
  BOOST_CHECK(l2norm < 1.0e-05);
  
}

BOOST_AUTO_TEST_SUITE_END()


