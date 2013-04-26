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

#ifndef _component_count_visitor_hpp_
#define _component_count_visitor_hpp_

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class ComponentCountVisitor {
public:

  /// Default constructor.
  ComponentCountVisitor(math::Matrix * matrix) : count_ (0), matrix_(matrix){};

  /// Destructor
  virtual ~ComponentCountVisitor(void){
      matrix_ = new math::Matrix(count_);
  };

  void increment(){++count_;};

private:
  int                  count_;
  math::Matrix       * matrix_;
};

} // namespace math
} // namespace gridpack

#endif
