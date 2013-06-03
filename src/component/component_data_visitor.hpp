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

class ComponentDataVisitor {
public:

  /// Default constructor.
  ComponentDataVisitor() {};

  /// Destructor
  virtual ~ComponentDataVisitor(void){};

  virtual void setData(MatrixInterface & interface);
  virtual void mapData(MatrixInterface & interface);

};

} // namespace math
} // namespace gridpack

#endif
