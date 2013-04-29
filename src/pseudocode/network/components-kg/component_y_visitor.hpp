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

#ifndef _component_y_visitor_hpp_
#define _component_y_visitor_hpp_

namespace gridpack {
namespace network {

#include <iostream>
class MatrixInterface;


// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class ComponentYVisitor {
public:

  /// Default constructor.
    ComponentYVisitor(void) : interface_(NULL){};
    virtual ~ComponentYVisitor(void){};

  void setInterfaceData(math::Matrix * interface){interface_ = interface;};
  MatrixInterface * getInterfaceData() const {return interface_;};
private:
  MatrixInterface   * interface_;
};

} // namespace math
} // namespace gridpack

#endif
