// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   vector_implementation.h
 * @author William A. Perkins
 * @date   Tue Mar 26 09:59:32 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _vector_implementation_h_
#define _vector_implementation_h_


#include "gridpack/parallel/distributable.h"
#include "gridpack/utlity/uncopyable.h"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class VectorImplementation
// -------------------------------------------------------------
class VectorImplementation {
public:

  /// Default constructor.
  VectorImplementation(const parallel::Distribution& dist);

  /// Destructor
  ~VectorImplementation(void);

  /// Set an individual element
  void set_element(const int& i, const double& x)
  {
    this->set_element_(i, x);
  }

  /// Set an several elements
  void set_elements(cont int& n, const int *i, const double *x)
  {
    this->set_elements_(n, i, x);
  }

  /// Add to an individual element
  void add_element(const int& i, const double& x)
  {
    this->add_element_(i, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const double *x)
  {
    this->add_elements_(n, i, x);
  }

  /// Make all the elements zero
  void zero(void)
  {
    this->zero_();
  }

  // FIXME: more ...

protected:

  /// Set an individual element (specialized)
  virtual void set_element_(const int& i, const double& x) = 0;

  /// Set an several elements (specialized)
  virtual void set_elements_(cont int& n, const int *i, const double *x) = 0;

  /// Add to an individual element (specialized)
  virtual void add_element_(const int& i, const double& x) = 0;

  /// Add to an several elements (specialized)
  virtual void add_elements_(const int& n, const int *i, const double *x) = 0;

  /// Make all the elements zero (specialized)
  virtual void zero_(void) = 0;
};

} // namespace math
} // namespace gridpack



#endif
