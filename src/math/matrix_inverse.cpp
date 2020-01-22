// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   matrix_inverse.cpp
 * @author William A. Perkins
 * @date   2015-08-18 14:14:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 25, 2014 by William A. Perkins
// -------------------------------------------------------------


#include <iostream>
#include <gridpack/parallel/parallel.hpp>
#include <gridpack/environment/environment.hpp>
#include <gridpack/timer/coarse_timer.hpp>
#include <gridpack/utilities/exception.hpp>
#include <gridpack/math/math.hpp>

using namespace gridpack;

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  Environment env(argc, argv);
  parallel::Communicator world;

  std::string cinput("input.xml");
  if (argc > 1) {
    cinput = argv[1];
  }
  boost::scoped_ptr<utility::Configuration> 
    config(utility::Configuration::configuration());
  config->enableLogging(&std::cout);
  if (!config->open(cinput, world)) {
    std::cerr << argv[0] << ": error: cannot open configuration " 
              << "\"" << cinput << "\""
              << std::endl;
    return 3;
  }
  utility::Configuration::CursorPtr 
    cursor(config->getCursor("MatrixInverse"));
  std::string inmatrix;
  bool dodirect;
  inmatrix = cursor->get("Matrix", "not-a-file");
  dodirect = cursor->get("Direct", true);

  utility::CoarseTimer 
    *timer(utility::CoarseTimer::instance());
  int t_setup(timer->createCategory("Setup"));
  int t_solve(timer->createCategory("Solve"));

  try {
    timer->start(t_setup);
    boost::scoped_ptr<math::Matrix> A, I, Ainv;
    A.reset(math::matrixLoadBinary<math::Matrix::TheType, math::Matrix::IdxType>(world, 
                                                                                 inmatrix.c_str()));
    I.reset(new math::Matrix(A->communicator(), 
                                       A->localRows(), A->localCols(), 
                                       math::Dense));
    I->identity();

    timer->stop(t_setup);

    if (dodirect) {
      boost::scoped_ptr<math::LinearMatrixSolver> 
        solver(new math::LinearMatrixSolver(*A));
      solver->configure(cursor);
      timer->start(t_solve);
      Ainv.reset(solver->solve(*I));
      timer->stop(t_solve);
    } else {
      boost::scoped_ptr<math::LinearSolver> 
        solver(new math::LinearSolver(*A));
      solver->configure(cursor);
      timer->start(t_solve);
      Ainv.reset(solver->solve(*I));
      timer->stop(t_solve);
    }
  } catch (const Exception& e) {
    std::cerr << argv[0] << ": error: " << e.what() << std::endl;
    return 2;
  }
  
  timer->dump();
  return 0;
}
