// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   implementation_visitor.hpp
 * @author William A. Perkins
 * @date   Mon Apr  1 09:26:28 2013
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

namespace gridpack {
namespace math {

class MatrixImplementation;
class VectorImplementation;
class LinearSolverImplementation;

class PETScVectorImplentation;
class PETScMatrixImplementation;
class PETScLinearSolverImplementation;

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

class ImplementationVisitor {
public:

  /// Default constructor.
  ImplementationVisitor(void);

  /// Destructor
  virtual ~ImplementationVisitor(void);

  /// The default visit (should just assert or do nothing)
  virtual void visit(MatrixImplementation);
  virtual void visit(VectorImplementation);
  virtual void visit(PETScVectorImplentation);
  virtual void visit(PETScMatrixImplementation);
};

} // namespace math
} // namespace gridpack

#endif
