// -------------------------------------------------------------
/**
 * @file   component_visitor.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef matrix_interface_hpp_
#define matrix_interface_hpp_

#include <iostream>

namespace gridpack {
namespace network {


// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

template <typename T>
class MatrixInterface {
public:

  /// Default constructor.
    MatrixInterface(int ni, int nj, int mi, int mj) :
        ni_(ni),
        nj_(nj),
        mi_(mi),
        mj_(mj_),
        source_i_(0),
        source_j_(0),
        size_((nj - ni)*(mj - ni)),
        x_(new T[size_])
    {}

    MatrixInterface() :
        ni_(0),
        nj_(0),
        mi_(0),
        mj_(0),
        source_i_(0),
        source_j_(0),
        size_(0),
        x_(NULL)
    {}

    /// Destructor
  virtual ~MatrixInterface(void);
  void setSourceIndices(int i, int j) {
      source_i_       = i;
      source_j_       = j;
  }

  void setIndex(int i, int j, T value);

  void setMatrixData(T * data) {
      for (int i = 0; i < size_; i++) {
          data[i]  = x_[i];
      }
  }
protected:
private:
  int                ni_;
  int                nj_;
  int                mi_;
  int                mj_
  int                source_i_;
  int                source_j_;
  int                size_;
  T                * x_;

};  /* END MATRIX_INTERFACE */

} // namespace math
} // namespace gridpack

#endif
