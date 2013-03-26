// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   matrix_implementation.h
 * @author William A. Perkins
 * @date   Mon Mar 25 12:10:01 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 25, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _matrix_implementation_h_
#define _matrix_implementation_h_

#include "gridpack/parallel/distributable.hpp"
#include "gridpack/utlity/uncopyable.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------
class MatrixImplementation 
  : public parallel::Distributable, 
    private utility::UnCopyable 
{
public:

  /// Default constructor.
  MatrixImplementation(const parallel::Distribution& dist, 
                       const int& rows, const int& cols);

  /// Destructor
  virtual ~MatrixImplementation(void);

  /// Set an individual element
  void set_element(const int& i, const int& j, const double& x)
  {
    this->set_element_(i, j, x);
  }

  /// Set an several elements
  void set_elements(cont int& n, const int *i, const int *j, const double *x)
  {
    this->set_elements_(n, i, j, x);
  }

  /// Set all elements in a row
  void set_row(const int& nj, const int& i, const int *j, const double *x)
  {
    this->set_row_(nj, i, j, x);
  }

  /// Set all elements in a row
  void set_region(const int& ni, const int& nj, const int *i, const int *j, const double *x)
  {
    this->set_row_(ni, nj, i, j, x);
  }

  /// Add to an individual element
  void add_element(const int& i, const int& j, const double& x)
  {
    this->add_element_(i, j, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const int *j, const double *x)
  {
    this->add_elements_(n, i, j, x);
  }

  /// Add to all elements in a row
  void add_row(const int& nj, const int& i, const int *j, const double *x)
  {
    this->add_row_(nj, i, j, x);
  }

protected:

  /// Set an individual element
  virtual void set_element_(const int& i, const int& j, const double& x) = 0;

  /// Set an several element
  virtual void set_elements_(const int *i, const int *j, const double *x) = 0;

  /// Set all elements in a row
  virtual void set_row_(const int& i, const int *j, const double *x) = 0;

  /// Add to  an individual element
  virtual void add_element_(const int& i, const int& j, const double& x) = 0;

  /// Add to  an several element
  virtual void add_elements_(const int *i, const int *j, const double *x) = 0;

  /// Add to  all elements in a row
  virtual void add_row_(const int& i, const int *j, const double *x) = 0;

  
};

} // namespace utility
} // namespace gridpack



#endif
