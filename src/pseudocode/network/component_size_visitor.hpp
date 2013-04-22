// -------------------------------------------------------------
/**
 * @file   component_visitor.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief A component visitor will access the size and shape of
 * a components matrix contribution. These values are reported
 * to the calling function.
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _component_visitor_hpp_
#define _component_visitor_hpp_

namespace gridpack {
namespace network {

class MatrixInterface;


// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class ComponentVisitor {
public:

  /// Default constructor.
  ComponentVisitor(math::Matrix & matrix, int ni, int nj, int i, int j) :
      ni_(ni),
      nj_(nj),
      i_(i),
      j_(j)
  {};

  /// Destructor
  virtual ~ComponentVisitor(void){};

  virtual void getSize(MatrixInterface);
  virtual void mapData(MatrixInterface & interface) {};

  math::MatrixImplementation &  matrixImpl_;
  // TODO: are these global or local indices, I'm operating under the assumption that they are global
  int                  ni_;
  int                  nj_;
  int                  i_;
  int                  j_;
private:
};

} // namespace math
} // namespace gridpack

#endif
