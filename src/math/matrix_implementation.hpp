// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix_implementation.h
 * @author William A. Perkins
 * @date   2015-08-11 15:50:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 25, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _matrix_implementation_h_
#define _matrix_implementation_h_

#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/complex.hpp>
#include <gridpack/math/matrix_interface.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class MatrixImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed,
    public BaseMatrixInterface<T, I>
{
public:

  /// Default constructor.
  MatrixImplementation(const parallel::Communicator& comm)
    : utility::Uncopyable(), parallel::Distributed(comm)
  {
  }

  /// Destructor
  virtual ~MatrixImplementation(void)
  {
  }

  /// Make an exact replica of this instance
  MatrixImplementation *clone(void) const
  {
    return this->p_clone();
  }

  /// Make a sequential copy of this instance local to this processor
  MatrixImplementation *localClone(void) const
  {
    return this->p_localClone();
  }
  


protected:

  /// Make an exact replica of this instance (specialized)
  virtual MatrixImplementation *p_clone(void) const = 0;

  /// Make a sequential copy of this instance local to this processor (specialized)
  virtual MatrixImplementation *p_localClone(void) const = 0;
};

} // namespace utility
} // namespace gridpack



#endif
