// -------------------------------------------------------------
/**
 * @file   matrixc.hpp
 * @author William A. Perkins
 * @date   2014-05-29 10:27:53 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 27, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#include <gridpack/math/matrix.hpp>

extern "C" void
matrix_initialize(gridpack::math::Matrix **m,
                  gridpack::parallel::Communicator *comm,
                  int local_rows, int local_cols, int nz_per_row)
{
  *m = new gridpack::math::Matrix(*comm, local_rows, local_cols, nz_per_row);
}

extern "C" void
matrix_finalize(gridpack::math::Matrix **m)
{
  delete *m;
  *m = NULL;
}

extern "C" int
matrix_rows(gridpack::math::Matrix *m)
{
  return m->rows();
}

extern "C" int
matrix_local_rows(gridpack::math::Matrix *m)
{
  return m->localRows();
}

extern "C" void
matrix_local_row_range(gridpack::math::Matrix *m, int *lo, int *hi)
{
  m->localRowRange(*lo, *hi);
}

extern "C" int
matrix_cols(gridpack::math::Matrix *m)
{
  return m->cols();
}

extern "C" void 
matrix_set_element(gridpack::math::Matrix *m, 
                   int i, int j, 
                   gridpack::ComplexType *x)
{
  m->setElement(i, j, *x);
}

extern "C" void 
matrix_add_element(gridpack::math::Matrix *m, 
                   int i, int j, 
                   gridpack::ComplexType *x)
{
  m->addElement(i, j, *x);
}

extern "C" void 
matrix_get_element(gridpack::math::Matrix *m, 
                   int i, int j, 
                   gridpack::ComplexType *x)
{
  m->getElement(i, j, *x);
}

extern "C" void
matrix_identity(gridpack::math::Matrix *m)
{
  m->identity();
}

extern "C" void
matrix_real(gridpack::math::Matrix *m)
{
  m->real();
}

extern "C" void
matrix_imaginary(gridpack::math::Matrix *m)
{
  m->imaginary();
}


extern "C" void
matrix_conjugate(gridpack::math::Matrix *m)
{
  m->conjugate();
}

extern "C" gridpack::ComplexType
matrix_norm2(gridpack::math::Matrix *m)
{
  return m->norm2();
}

extern "C" void
matrix_ready(gridpack::math::Matrix *m)
{
  m->ready();
}

extern "C" gridpack::math::Matrix *
matrix_clone(gridpack::math::Matrix *m)
{
  return m->clone();
}

extern "C" void
matrix_print(gridpack::math::Matrix *m)
{
  m->print();
}

extern "C" void
matrix_scale(gridpack::math::Matrix *m, 
             gridpack::ComplexType *x)
{
  m->scale(*x);
}

extern "C" void
matrix_multiply_diagonal(gridpack::math::Matrix *m, 
                         gridpack::math::Vector *x)
{
  m->multiplyDiagonal(*x);
}

extern "C" void
matrix_add(gridpack::math::Matrix *m,
           gridpack::math::Matrix *a)
{
  m->add(*a);
}


extern "C" void
matrix_zero(gridpack::math::Matrix *m)
{
  m->zero();
}





