// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   implementation_visitor.hpp
 * @author William A. Perkins
 * @date   2015-03-05 12:51:12 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _implementation_visitor_hpp_
#define _implementation_visitor_hpp_

#include <gridpack/utilities/uncopyable.hpp>

namespace gridpack {
namespace math {

class PetscVectorWrapper;
class PetscMatrixWrapper;

// -------------------------------------------------------------
//  class ImplementationVisitor
// -------------------------------------------------------------
/**
 * This is a cyclic visitor for the various math implementation
 * classes.  It is intended to be use as a parent class for things
 * that are used to extract implementation-specific information from
 * an implementation agnostic class.
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class ImplementationVisitor : private utility::Uncopyable {
public:

  /// Default constructor.
  ImplementationVisitor(void);

  /// Destructor
  virtual ~ImplementationVisitor(void);

  /// The default visit (should just assert or do nothing)
  virtual void visit(PetscVectorWrapper&);
  virtual void visit(PetscMatrixWrapper&);
};

// -------------------------------------------------------------
//  class ConstImplementationVisitor
// -------------------------------------------------------------
/**
 * This is const version of a cyclic visitor for the various math
 * implementation classes.  It is intended to be used as a parent class
 * for things that are used to extract implementation-specific
 * information from an const implementation agnostic class.
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */
class ConstImplementationVisitor : private utility::Uncopyable {
public:

  /// Default constructor.
  ConstImplementationVisitor(void);

  /// Destructor
  virtual ~ConstImplementationVisitor(void);

  /// The default visit, const version (should just assert or do nothing)
  virtual void visit(const PetscVectorWrapper&);
  virtual void visit(const PetscMatrixWrapper&);
};

} // namespace math
} // namespace gridpack

#endif
