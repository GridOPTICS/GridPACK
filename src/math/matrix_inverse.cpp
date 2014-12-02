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
 * @date   2014-11-26 11:22:17 d3g096
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
#include <gridpack/timer/coarse_timer.hpp>
#include <gridpack/utilities/exception.hpp>
#include <gridpack/math/math.hpp>

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  gridpack::math::Initialize();

  std::string cinput("input.xml");
  if (argc > 1) {
    cinput = argv[1];
  }
  boost::scoped_ptr<gridpack::utility::Configuration> 
    config(gridpack::utility::Configuration::configuration());
  config->enableLogging(&std::cout);
  if (!config->open(cinput, world)) {
    std::cerr << argv[0] << ": error: cannot open configuration " 
              << "\"" << cinput << "\""
              << std::endl;
    return 3;
  }
  gridpack::utility::Configuration::CursorPtr 
    cursor(config->getCursor("MatrixInverse"));
  std::string inmatrix;
  bool dodirect;
  inmatrix = cursor->get("Matrix", "not-a-file");
  dodirect = cursor->get("Direct", true);

  gridpack::utility::CoarseTimer 
    *timer(gridpack::utility::CoarseTimer::instance());
  int t_setup(timer->createCategory("Setup"));
  int t_solve(timer->createCategory("Solve"));

  try {
    timer->start(t_setup);
    boost::scoped_ptr<gridpack::math::Matrix> A, I, Ainv;
    A.reset(gridpack::math::matrixLoadBinary(world, inmatrix.c_str()));
    I.reset(new gridpack::math::Matrix(A->communicator(), 
                                       A->localRows(), A->localCols(), 
                                       gridpack::math::Matrix::Dense));
    I->identity();

    timer->stop(t_setup);

    if (dodirect) {
      boost::scoped_ptr<gridpack::math::LinearMatrixSolver> 
        solver(new gridpack::math::LinearMatrixSolver(*A));
      solver->configure(cursor);
      timer->start(t_solve);
      Ainv.reset(solver->solve(*I));
      timer->stop(t_solve);
    } else {
      boost::scoped_ptr<gridpack::math::LinearSolver> 
        solver(new gridpack::math::LinearSolver(*A));
      solver->configure(cursor);
      timer->start(t_solve);
      Ainv.reset(solver->solve(*I));
      timer->stop(t_solve);
    }
  } catch (const gridpack::Exception& e) {
    std::cerr << argv[0] << ": error: " << e.what() << std::endl;
    return 2;
  }
  
  gridpack::math::Finalize();
  timer->dump();
  return 0;
}
