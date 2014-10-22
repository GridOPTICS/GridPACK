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
 * @file   implementation_visitable.hpp
 * @author William A. Perkins
 * @date   2014-10-22 09:03:08 d3g096
 * 
 * @brief Declaration of the abstract ImplementationVisitable class.
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 21, 2014 by William A. Perkins
// Last Change: 2014-10-21 14:41:29 d3g096
// -------------------------------------------------------------


#ifndef _implementation_visitable_hpp_
#define _implementation_visitable_hpp_

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class ImplementationVisitable
// -------------------------------------------------------------
/// Interface for classes that accept ImplementationVisitor's
/**
 * Classes that can be visited by an ImplementationVisitor can
 * subclass from this.  It is still necessary to provide
 * specialization to recognize specific classes in the
 * ImplementationVisitor class.  
 * 
 */
class ImplementationVisitable {
protected:

  /// Protected copy constructor to avoid unwanted copies.
  ImplementationVisitable(const ImplementationVisitable& old);

public:

  /// Default constructor.
  ImplementationVisitable(void)
  {}

  /// Destructor
  ~ImplementationVisitable(void)
  {}

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

protected:

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

};

} // namespace math
} // namespace gridpack
#endif
