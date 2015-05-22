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
 * @file   fall_back_matrix_methods.hpp
 * @author William A. Perkins
 * @date   2015-05-22 09:44:35 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _fall_back_matrix_methods_hpp_
#define _fall_back_matrix_methods_hpp_

#include "gridpack/parallel/communicator.hpp"
#include "matrix_interface.hpp"
#include "value_transfer.hpp"

namespace gridpack {
namespace math {
namespace fallback {



template <typename T, typename I>
double 
norm2(const parallel::Communicator& comm, const BaseMatrixInterface<T, I>& A)
{
  RealType result(0.0);
  I ncols(A.cols());
  I lo, hi;
  A.localRowRange(lo, hi);

  T v;
  RealType tmp;
  ValueTransfer<T, RealType> trans(1, &v, &tmp);
  for (I i = lo; i < hi; ++i) {
    for (I j = 0; j < ncols; ++j) {
      A.getElement(i, j, v);
      v *= v;
      trans.go();
      result += tmp;
    }
  }
  comm.sum(&result, 1);
  result = sqrt(result);
  return result;
}

} // namespace fallback
} // namespace math
} // namespace gridpack


#endif
