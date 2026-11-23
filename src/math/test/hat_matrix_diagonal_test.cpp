// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hat_matrix_diagonal_test.cpp
 * @brief  Unit test for hatMatrixDiagonal()
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <vector>
#include <boost/mpi/collectives.hpp>

#include "math.hpp"
#include "matrix.hpp"
#include "hat_matrix_diagonal.hpp"

#include "test_main.cpp"

BOOST_AUTO_TEST_SUITE(HatMatrixDiagonalTest)

BOOST_AUTO_TEST_CASE( TraceAndMethods )
{
  gridpack::parallel::Communicator world;
  const int local_rows(6), local_cols(4);
  int m(local_rows*world.size()), n(local_cols*world.size());

  // Two entries per row; the first n rows are diagonally dominant
  gridpack::math::Matrix H(world, local_rows, local_cols, gridpack::math::Sparse);
  int lo, hi;
  H.localRowRange(lo, hi);
  std::vector<double> w;
  for (int i = lo; i < hi; ++i) {
    H.setElement(i, i % n, gridpack::ComplexType(1.0 + 0.1*i, 0.0));
    H.setElement(i, (i + 1) % n, gridpack::ComplexType(-0.5, 0.0));
    w.push_back(1.0 + (i % 3));
  }
  H.ready();

  std::vector<double> d1, d2;
  gridpack::math::HatMatrixDiagonalInfo i1, i2;
  gridpack::math::hatMatrixDiagonal(H, w, d1, i1, gridpack::math::SelectedInversion);
  gridpack::math::hatMatrixDiagonal(H, w, d2, i2, gridpack::math::BlockedSolves);
  BOOST_CHECK_EQUAL(i2.method, "blocked solves");
  BOOST_CHECK_EQUAL(static_cast<int>(d1.size()), hi - lo);
  BOOST_CHECK_EQUAL(static_cast<int>(d2.size()), hi - lo);

  // The trace of the hat matrix is the number of columns
  double t1(0.0), t2(0.0), g1, g2;
  for (int i = 0; i < hi - lo; ++i) {
    BOOST_CHECK(d1[i] > 0.0 && d1[i] < 1.0 + 1.0e-12);
    BOOST_CHECK_SMALL(d1[i] - d2[i], 1.0e-10);
    t1 += d1[i];
    t2 += d2[i];
  }
  boost::mpi::all_reduce(world, t1, g1, std::plus<double>());
  boost::mpi::all_reduce(world, t2, g2, std::plus<double>());
  BOOST_CHECK_CLOSE(g1, static_cast<double>(n), 1.0e-8);
  BOOST_CHECK_CLOSE(g2, static_cast<double>(n), 1.0e-8);
  if (world.rank() == 0) {
    std::cout << "rows " << m << ", method: " << i1.method
              << (i1.note.empty() ? "" : " (" + i1.note + ")") << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
