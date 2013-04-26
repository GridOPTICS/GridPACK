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

#ifndef _component_h_visitor_hpp_
#define _component_h_visitor_hpp_

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

template <typename T>
class ComponentHVisitor {
public:

  /// Default constructor.
  ComponentHVisitor(math::Matrix * matrix)  :
      ni_(0),
      nj_(0),
      mi_(0),
      mj_(0),
      source_i_(0),
      source_j_(0),
      size_(0),
      x_(NULL),
      matrix_(matrix){};

  /// Destructor
  virtual ~ComponentHVisitor(void){};

  virtual void getComponentData(int * ni, int * nj, int * mi, int * mj, int * source_i,
          int * source_j, int * size, T * x) {
      ni_           = ni;
      nj_           = nj;
      mi_           = mi;
      mj_           = mj;
      source_i_     = source_i;
      source_j_     = source_j;
      size_         = size;
      x_            = x;
  }
  virtual void setMapData(int * ni, int * nj, int * mi, int * mj, int * source_i,
          int * source_j, int * size, T * x);
private:
  int                ni_;
  int                nj_;
  int                mi_;
  int                mj_
  int                source_i_;
  int                source_j_;
  int                size_;
  T                * x_;
  math::Matrix     * matrix_;

};

} // namespace math
} // namespace gridpack

#endif
