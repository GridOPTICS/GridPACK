// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   small_matrix_solve.cpp
 * @author William A. Perkins
 * @date   2013-10-23 10:23:50 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>

#include "gridpack/parallel/parallel.hpp"
#include "gridpack/utilities/exception.hpp"
#include "math.hpp"
#include "linear_solver.hpp"

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  boost::mpi::communicator self = world.split(world.rank());

  gridpack::math::Initialize();

  boost::scoped_ptr<gridpack::utility::Configuration> 
    config(gridpack::utility::Configuration::configuration());
  
  config->open("small_matrix_solve.xml", world);
  
  config.reset(config->getCursor("SmallMatrixSolve"));

  BOOST_ASSERT(config);

  boost::scoped_ptr<gridpack::math::Matrix> 
    A(new gridpack::math::Matrix(self, 9, 9)),
    B(new gridpack::math::Matrix(self, 9, 3)),
    Xorig(new gridpack::math::Matrix(self, 9, 3, gridpack::math::Matrix::Dense));

  A->setElement(0, 0, gridpack::ComplexType(0, -33.8085));
  A->setElement(0, 3, gridpack::ComplexType(0,  17.3611));
  A->setElement(1, 1, gridpack::ComplexType(0, -24.3472));
  A->setElement(1, 7, gridpack::ComplexType(0,  16));
  A->setElement(2, 2, gridpack::ComplexType(0, -22.5806));
  A->setElement(2, 5, gridpack::ComplexType(0,  17.0648));
  A->setElement(3, 0, gridpack::ComplexType(0,  17.3611));
  A->setElement(3, 3, gridpack::ComplexType(0, -39.9954));
  A->setElement(3, 4, gridpack::ComplexType(0,  10.8696));
  A->setElement(3, 8, gridpack::ComplexType(0,  11.7647));
  A->setElement(4, 3, gridpack::ComplexType(0,  10.8696));
  A->setElement(4, 4, gridpack::ComplexType(0.934, -17.0633));
  A->setElement(4, 5, gridpack::ComplexType(0,  5.88235));
  A->setElement(5, 2, gridpack::ComplexType(0,  17.0648));
  A->setElement(5, 4, gridpack::ComplexType(0,  5.88235));
  A->setElement(5, 5, gridpack::ComplexType(0, -32.8678));
  A->setElement(5, 6, gridpack::ComplexType(0,  9.92063));
  A->setElement(6, 5, gridpack::ComplexType(0,  9.92063));
  A->setElement(6, 6, gridpack::ComplexType(1.03854, -24.173));
  A->setElement(6, 7, gridpack::ComplexType(0,  13.8889));
  A->setElement(7, 1, gridpack::ComplexType(0,  16));
  A->setElement(7, 6, gridpack::ComplexType(0,  13.8889));
  A->setElement(7, 7, gridpack::ComplexType(0, -36.1001));
  A->setElement(7, 8, gridpack::ComplexType(0,  6.21118));
  A->setElement(8, 3, gridpack::ComplexType(0,  11.7647));
  A->setElement(8, 7, gridpack::ComplexType(0,  6.21118));
  A->setElement(8, 8, gridpack::ComplexType(1.33901, -18.5115));
  A->ready();
  A->print();
  A->save("small_matrix_solve.mat");

  B->setElement(0, 0, gridpack::ComplexType(0,  16.4474));
  B->setElement(1, 1, gridpack::ComplexType(0,   8.34725));
  B->setElement(2, 2, gridpack::ComplexType(0,   5.51572));
  B->ready();
  B->print();

  // This is the expected answer
  Xorig->setElement(0, 0, gridpack::ComplexType(-0.802694,  0.0362709)); 
  Xorig->setElement(0, 1, gridpack::ComplexType(-0.0826026, 0.0215508));
  Xorig->setElement(0, 2, gridpack::ComplexType(-0.066821,  0.0163025));
  Xorig->setElement(1, 0, gridpack::ComplexType(-0.16276 ,  0.0424636 ));
  Xorig->setElement(1, 1, gridpack::ComplexType(-0.656183,  0.032619));
  Xorig->setElement(1, 2, gridpack::ComplexType(-0.117946,  0.0235926));
  Xorig->setElement(2, 0, gridpack::ComplexType(-0.199254 , 0.0486126));
  Xorig->setElement(2, 1, gridpack::ComplexType(-0.178494,  0.035704));
  Xorig->setElement(2, 2, gridpack::ComplexType(-0.55177,   0.0277253));
  Xorig->setElement(3, 0, gridpack::ComplexType(-0.615773 , 0.0706327));
  Xorig->setElement(3, 1, gridpack::ComplexType(-0.160858,  0.0419673));
  Xorig->setElement(3, 2, gridpack::ComplexType(-0.130125,  0.031747));
  Xorig->setElement(4, 0, gridpack::ComplexType(-0.478041 , 0.0933363));
  Xorig->setElement(4, 1, gridpack::ComplexType(-0.180994,  0.052928));
  Xorig->setElement(4, 2, gridpack::ComplexType(-0.220702,  0.0449514));
  Xorig->setElement(5, 0, gridpack::ComplexType(-0.263657 , 0.0643252)); 
  Xorig->setElement(5, 1, gridpack::ComplexType(-0.236187,  0.0472443));
  Xorig->setElement(5, 2, gridpack::ComplexType(-0.406892,  0.0366867));
  Xorig->setElement(6, 0, gridpack::ComplexType(-0.247322 , 0.0741512));
  Xorig->setElement(6, 1, gridpack::ComplexType(-0.368152,  0.0637251));
  Xorig->setElement(6, 2, gridpack::ComplexType(-0.268082,  0.0472012));
  Xorig->setElement(7, 0, gridpack::ComplexType(-0.247672 , 0.0646169));
  Xorig->setElement(7, 1, gridpack::ComplexType(-0.476813,  0.0496364));
  Xorig->setElement(7, 2, gridpack::ComplexType(-0.179478,  0.0359009));
  Xorig->setElement(8, 0, gridpack::ComplexType(-0.467188 , 0.100364));
  Xorig->setElement(8, 1, gridpack::ComplexType(-0.257734,  0.0619692));
  Xorig->setElement(8, 2, gridpack::ComplexType(-0.139857,  0.0423387));
  Xorig->ready();
  Xorig->print();

  boost::scoped_ptr<gridpack::math::LinearSolver> 
    solver(new gridpack::math::LinearSolver(*A));

  solver->configurationKey("LinearMatrixSolver");
  solver->configure(config.get());

  boost::scoped_ptr<gridpack::math::Matrix> 
    X(solver->solve(*B));

  X->scale(-1.0);
  X->add(*Xorig);
  X->print();

  std::cout << world.rank() 
            << ": Solution norm2: " << X->norm2() 
            << std::endl;

  solver.reset();
  A.reset();
  B.reset();
  Xorig.reset();
  X.reset();

  gridpack::math::Finalize();
  
  return 0;
}
