// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_misc.hpp
 * @author William A. Perkins
 * @date   2015-07-23 09:12:20 d3g096
 * 
 * @brief Stuff necessary for PETSc math implementation but defies
 * classification.
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _petsc_misc_hpp_
#define _petsc_misc_hpp_

#include <petscsys.h>
#include <vector>
#include <algorithm>
#include "gridpack/utilities/complex.hpp"

/// Scale a complex DENSE matrix
extern PetscErrorCode sillyMatScaleComplex(Mat A, const gridpack::ComplexType& px);

// -------------------------------------------------------------
// sortPermutation
// 
// This stuff is used because ISSortPermutation is not available prior
// to PETSc 3.6.
// -------------------------------------------------------------
template <typename T, typename I = int>
struct TupleCompare
{
  typedef std::pair<T, I> tuple;
  bool operator() (const tuple& lhs, const tuple& rhs)
  {
    return lhs.first < rhs.first;
  }
};

template <typename T, typename I>
std::vector<I>
sortPermutation(const std::vector<T>& x) 
{
  size_t n(x.size());

  typedef std::pair<T, I> tuple;
  typedef std::vector<tuple> tvector;

  tvector xp(n);
  for (int i = 0; i < n; i++) {
    xp[i] = tuple(x[i], i);
  }

  std::sort(xp.begin(), xp.end(), TupleCompare<T, I>());
  
  std::vector<I> result(n);
  for (int i = 0; i < n; i++) {
    result[i] = xp[i].second;
  }
  return result;
}

#endif
