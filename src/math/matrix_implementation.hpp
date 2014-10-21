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
 * @date   2014-10-21 12:52:31 d3g096
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

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------
class MatrixImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed,
    public BaseMatrixInterface<ComplexType>
{
public:

  /// Default constructor.
  MatrixImplementation(const parallel::Communicator& comm);

  /// Destructor
  virtual ~MatrixImplementation(void);

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->p_accept(visitor);
  }

  /// Allow visits by const implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

  /// Make an exact replica of this instance
  MatrixImplementation *clone(void) const
  {
    return this->p_clone();
  }


protected:

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

  /// Make an exact replica of this instance (specialized)
  virtual MatrixImplementation *p_clone(void) const = 0;
};

} // namespace utility
} // namespace gridpack



#endif
