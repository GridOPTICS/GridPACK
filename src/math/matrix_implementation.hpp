// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   matrix_implementation.h
 * @author William A. Perkins
 * @date   2013-06-12 10:27:22 d3g096
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

#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/complex.hpp>

namespace gridpack {
namespace math {

class ImplementationVisitor;
class ConstImplementationVisitor;

// -------------------------------------------------------------
//  class MatrixImplementation
// -------------------------------------------------------------
class MatrixImplementation 
  : private utility::Uncopyable,
    public parallel::Distributed
{
public:

  /// Default constructor.
  MatrixImplementation(const parallel::Communicator& comm);

  /// Destructor
  virtual ~MatrixImplementation(void);

  /// Get the global index range of the locally owned rows
  void local_row_range(int& lo, int& hi) const 
  {
    this->p_local_row_range(lo, hi);
  }

  /// Get the total number of rows in this matrix
  int rows(void) const
  {
    return this->p_rows();
  }

  /// Get the number of local rows in this matirx
  int local_rows(void) const
  {
    return this->p_local_rows();
  }

  /// Get the number of columns in this matrix
  int cols(void) const
  {
    return this->p_cols();
  }

  /// Set an individual element
  void set_element(const int& i, const int& j, const ComplexType& x)
  {
    this->p_set_element(i, j, x);
  }

  /// Set an several elements
  void set_elements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    this->p_set_elements(n, i, j, x);
  }

  // /// Set all elements in a row
  // void set_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_set_row(nj, i, j, x);
  // }

  // /// Set all elements in a row
  // void set_region(const int& ni, const int& nj, const int *i, const int *j, const ComplexType *x)
  // {
  //   this->p_set_region(ni, nj, i, j, x);
  // }

  /// Add to an individual element
  void add_element(const int& i, const int& j, const ComplexType& x)
  {
    this->p_add_element(i, j, x);
  }

  /// Add to an several elements
  void add_elements(const int& n, const int *i, const int *j, const ComplexType *x)
  {
    this->p_add_elements(n, i, j, x);
  }

  // /// Add to all elements in a row
  // void add_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_add_row(nj, i, j, x);
  // }

  /// Get an individual element
  void get_element(const int& i, const int& j, ComplexType& x) const
  {
    this->p_get_element(i, j, x);
  }

  /// Get an several elements
  void get_elements(const int& n, const int *i, const int *j, ComplexType *x) const
  {
    this->p_get_elements(n, i, j, x);
  }

  // /// Get all elements in a row
  // void get_row(const int& nj, const int& i, const int *j, const ComplexType *x)
  // {
  //   this->p_get_row(nj, i, j, x);
  // }

  // /// Get all elements in a row
  // void get_region(const int& ni, const int& nj, 
  //                 const int *i, const int *j, const ComplexType *x)
  // {
  //   this->p_get_region(ni, nj, i, j, x);
  // }

  /// Indicate the matrix is ready to use
  void ready(void)
  {
    this->p_ready();
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->p_accept(visitor);
  }

  /// Allow visits by const implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

  /// Make an exact replica of this instance
  MatrixImplementation *clone(void) const
  {
    return this->p_clone();
  }


protected:

  /// Get the global index range of the locally owned rows (specialized)
  virtual void p_local_row_range(int& lo, int& hi) const = 0;

  /// Get the total number of rows in this matrix (specialized)
  virtual int p_rows(void) const = 0;

  /// Get the number of local rows in this matirx (specialized)
  virtual int p_local_rows(void) const = 0;

  /// Get the number of columns in this matrix (specialized)
  virtual int p_cols(void) const = 0;

  /// Set an individual element
  virtual void p_set_element(const int& i, const int& j, const ComplexType& x) = 0;

  /// Set an several element
  virtual void p_set_elements(const int& n, const int *i, const int *j, 
                              const ComplexType *x) = 0;

  // /// Set all elements in a row
  // virtual void p_set_row(const int& i, const int *j, const ComplexType *x) = 0;

  // /// Set all elements in a region
  // virtual void p_set_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, const ComplexType *x) = 0;

  /// Add to  an individual element
  virtual void p_add_element(const int& i, const int& j, const ComplexType& x) = 0;

  /// Add to  an several element
  virtual void p_add_elements(const int& n, const int *i, const int *j, 
                              const ComplexType *x) = 0;

  // /// Add to  all elements in a row
  // virtual void p_add_row(const int& i, const int *j, const ComplexType *x) = 0;

  /// Get an individual element
  virtual void p_get_element(const int& i, const int& j, ComplexType& x) const = 0;

  /// Get an several element
  virtual void p_get_elements(const int& n, const int *i, const int *j, 
                              ComplexType *x) const = 0;

  // /// Get all elements in a row
  // virtual void p_get_row(const int& i, const int *j, ComplexType *x) const = 0;

  // /// Get all elements in a region
  // virtual void p_get_region(const int& ni, const int& nj, 
  //                          const int *i, const int *j, ComplexType *x) const = 0;


  /// Make this instance ready to use
  virtual void p_ready(void) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

  /// Make an exact replica of this instance (specialized)
  virtual MatrixImplementation *p_clone(void) const = 0;
};

} // namespace utility
} // namespace gridpack



#endif
