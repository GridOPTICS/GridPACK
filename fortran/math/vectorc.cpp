// -------------------------------------------------------------
// file: vectorc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 15, 2014 by William A. Perkins
// Last Change: 2014-10-01 07:04:42 d3g096
// -------------------------------------------------------------

#include <gridpack/parallel/communicator.hpp>
#include <gridpack/math/vector.hpp>

extern "C" void
vector_initialize(gridpack::math::Vector **v, 
                  gridpack::parallel::Communicator *comm,
                  int local_length)
{
  *v = new gridpack::math::Vector(*comm, local_length);
}

extern "C" void
vector_finalize(gridpack::math::Vector **v)
{
  delete (*v);
  *v = NULL;
}

extern "C" int
vector_size(gridpack::math::Vector *v)
{
  return v->size();
}

extern "C" int
vector_local_size(gridpack::math::Vector *v)
{
  return v->localSize();
}


extern "C" void
vector_local_index_range(gridpack::math::Vector *v, int *lo, int *hi)
{
  v->localIndexRange(*lo, *hi);
}

extern "C" void
vector_set_element(gridpack::math::Vector *v, 
                   int i, 
                   gridpack::ComplexType *x)
{
  v->setElement(i, *x);
}

extern "C" void
vector_add_element(gridpack::math::Vector *v, 
                   int i, 
                   gridpack::ComplexType *x)
{
  v->addElement(i, *x);
}

extern "C" void
vector_get_element(gridpack::math::Vector *v, 
                   int i, 
                   gridpack::ComplexType *x)
{
  v->getElement(i, *x);
}

extern "C" void
vector_get_all_elements(gridpack::math::Vector *v, 
                        gridpack::ComplexType *x)
{
  // Note there is no check on length of array pointed to by x; how
  // can there be?
  v->getAllElements(x);
}

extern "C" void
vector_zero(gridpack::math::Vector *v)
{
  v->zero();
}
  
extern "C" void
vector_fill(gridpack::math::Vector *v, gridpack::ComplexType x)
{
  v->fill(x);
}
  

extern "C" gridpack::ComplexType
vector_norm1(gridpack::math::Vector *v)
{
  return v->norm1();
}

extern "C" gridpack::ComplexType
vector_norm2(gridpack::math::Vector *v)
{
  return v->norm2();
}

extern "C" gridpack::ComplexType
vector_norm_infinity(gridpack::math::Vector *v)
{
  return v->normInfinity();
}

extern "C" void
vector_ready(gridpack::math::Vector *v)
{
  v->ready();
}
  
extern "C" void
vector_print(gridpack::math::Vector *v)
{
  v->print();
}

extern "C" void
vector_scale(gridpack::math::Vector *v, gridpack::ComplexType x)
{
  v->scale(x);
}
  
extern "C" void
vector_abs(gridpack::math::Vector *v)
{
  v->abs();
}

extern "C" void
vector_real(gridpack::math::Vector *v)
{
  v->real();
}

extern "C" void
vector_imaginary(gridpack::math::Vector *v)
{
  v->imaginary();
}

extern "C" void
vector_conjugate(gridpack::math::Vector *v)
{
  v->conjugate();
}

extern "C" void
vector_exp(gridpack::math::Vector *v)
{
  v->exp();
}

extern "C" gridpack::math::Vector *
vector_clone(gridpack::math::Vector *v)
{
  return v->clone();
}

