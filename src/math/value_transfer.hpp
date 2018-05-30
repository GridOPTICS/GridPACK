// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   value_transfer.hpp
 * @author William A. Perkins
 * @date   2016-07-14 13:18:45 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 20, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _value_transfer_hpp_
#define _value_transfer_hpp_

#include <boost/shared_array.hpp>
#include <boost/static_assert.hpp>

#include "numeric_type_check.hpp"

#include "gridpack/utilities/uncopyable.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class valueTransfer
// -------------------------------------------------------------
/**
 * This provides a consistent way to transfer values between RealType
 * and ComplexType arrays, being (hopefully) smart about when memory
 * is allocated. The transfer works both ways: complex <-> real. This
 * is not approprate for the case of real math with a complex-based
 * math library!
 * 
 */
template <typename FromType, typename ToType>
class ValueTransfer 
  : private utility::Uncopyable
{
public:

  /// Default constructor.
  ValueTransfer(const unsigned int& from_size, FromType* from, ToType* to = NULL)
    : utility::Uncopyable(),
      p_fromSize(from_size), p_from(from), 
      p_to()
  {
    if (to != NULL) {
      p_to.reset(to, null_deleter());
    }
  }

  /// Destructor
  ~ValueTransfer(void) 
  {  }

  /// Get the length of the "to" buffer
  unsigned int size(void) const
  {
    return p_computeToSize();
  }

  /// do whatever is necessary to set up and copy values
  void go(void)
  {
    this->p_setup();
    this->p_copy();
  }

  /// Get a pointer to the "to" buffer
  ToType *to(void)
  {
    return p_to.get();
  }

protected:

  /// Two
  static const unsigned int TWO = 2;

  /// Are the from and to types the same?
  static const bool p_isSame = boost::is_same<FromType, ToType>::value;

  /// The size of the from buffer
  const unsigned int p_fromSize;

  /// The from buffer
  FromType *p_from;

  /// The "to" buffer
  boost::shared_array<ToType> p_to;

  /// Compute the size of the "to" buffer
  virtual inline unsigned int p_computeToSize(void) const; 

  /// Do the setup  
  inline void p_setup(void);

  /// Copy the value
  virtual inline void p_copy(void);

private:

  BOOST_STATIC_ASSERT(TypeCheck<FromType>::OK::value);
  BOOST_STATIC_ASSERT(TypeCheck<ToType>::OK::value);

  /// A thing to provide a null delete operation
  struct null_deleter {
    void operator()(void const *p) const {}
  };
  
  
};


// -------------------------------------------------------------
// ValueTransfer<>::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
ValueTransfer<FromType, ToType>::p_computeToSize(void) const
{
  BOOST_STATIC_ASSERT(p_isSame);
  return p_fromSize;
}

template <>
inline unsigned int
ValueTransfer<RealType, ComplexType>::p_computeToSize(void) const
{
  return p_fromSize/TWO;
}

template <>
inline unsigned int
ValueTransfer<ComplexType, RealType>::p_computeToSize(void) const
{
  return p_fromSize*TWO;
}

// -------------------------------------------------------------
// ValueTransfer<>::p_setup
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void  
ValueTransfer<FromType, ToType>::p_setup()
{
  if (!p_to) {
    p_to.reset(new ToType[this->size()]);
  }
}

template <>
inline void  
ValueTransfer<RealType, RealType>::p_setup()
{
  if (!p_to) {
    p_to.reset(p_from, null_deleter());
  }
}

template <>
inline void  
ValueTransfer<ComplexType, ComplexType>::p_setup()
{
  if (!p_to) {
    p_to.reset(p_from, null_deleter());
  }
}


// -------------------------------------------------------------
// ValueTransfer<>::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
ValueTransfer<FromType, ToType>::p_copy(void)
{
  BOOST_STATIC_ASSERT(p_isSame);
  if (p_from != p_to.get()) 
    std::copy(p_from, p_from + p_fromSize, p_to.get());
}

template <>
inline void
ValueTransfer<RealType, ComplexType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; fidx += TWO, ++tidx) {
    p_to.get()[tidx] = ComplexType(p_from[fidx], p_from[fidx+1]);
  }
}

template <>
inline void
ValueTransfer<ComplexType, RealType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; ++fidx, tidx += TWO) {
    p_to.get()[tidx] = std::real(p_from[fidx]);
    p_to.get()[tidx+1] = std::imag(p_from[fidx]);
  }
}

// -------------------------------------------------------------
//  class ValueTransferToLibrary
// -------------------------------------------------------------
template <typename FromType, typename ToType>
class ValueTransferToLibrary 
  : public ValueTransfer<FromType, ToType>
{
public:

  /// Default constructor.
  ValueTransferToLibrary(const unsigned int& from_size, FromType* from, ToType* to = NULL)
    : ValueTransfer<FromType, ToType>(from_size, from, to)
  {}

  /// Destructor
  ~ValueTransferToLibrary(void) {}

protected:

  /// Compute the size of the "to" buffer
  inline unsigned int p_computeToSize(void) const;

  /// Copy the value
  inline void p_copy(void);
};

// -------------------------------------------------------------
// ValueTransferToLibrary<>::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
ValueTransferToLibrary<FromType, ToType>::p_computeToSize(void) const
{
  return ValueTransfer<FromType, ToType>::p_computeToSize();
}

template <>
inline unsigned int
ValueTransferToLibrary<RealType, ComplexType>::p_computeToSize(void) const
{
  return p_fromSize;
}


// -------------------------------------------------------------
// ValueTransferToLibrary<>::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
ValueTransferToLibrary<FromType, ToType>::p_copy(void)
{
  ValueTransfer<FromType, ToType>::p_copy();
}


template <>
inline void
ValueTransferToLibrary<RealType, ComplexType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; ++fidx, ++tidx) {
    p_to.get()[tidx] = p_from[fidx];
  }
}


// -------------------------------------------------------------
//  class ValueTransferFromLibrary
// -------------------------------------------------------------
template <typename FromType, typename ToType>
class ValueTransferFromLibrary 
  : public ValueTransfer<FromType, ToType>
{
public:

  /// Default constructor.
  ValueTransferFromLibrary(const unsigned int& from_size, FromType* from, ToType* to = NULL)
    : ValueTransfer<FromType, ToType>(from_size, from, to)
  {}

  /// Destructor
  ~ValueTransferFromLibrary(void) {}

protected:
  
  /// Compute the size of the "to" buffer
  inline unsigned int p_computeToSize(void) const;

  /// Copy the value
  inline void p_copy(void);
};

// -------------------------------------------------------------
// ValueTransferFromLibrary<>::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
ValueTransferFromLibrary<FromType, ToType>::p_computeToSize(void) const
{
  return ValueTransfer<FromType, ToType>::p_computeToSize();
}

template <>
inline unsigned int
ValueTransferFromLibrary<ComplexType, RealType>::p_computeToSize(void) const
{
  return p_fromSize;
}


// -------------------------------------------------------------
// ValueTransferFromLibrary<>::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
ValueTransferFromLibrary<FromType, ToType>::p_copy(void)
{
  ValueTransfer<FromType, ToType>::p_copy();
}


template <>
inline void
ValueTransferFromLibrary<ComplexType, RealType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (fidx = 0; fidx < p_fromSize; ++fidx, ++tidx) {
    p_to.get()[tidx] = std::real(p_from[fidx]);
  }
}


// -------------------------------------------------------------
//  class MatrixValueTransferToLibrary
// -------------------------------------------------------------
template <typename FromType, typename ToType>
class MatrixValueTransferToLibrary 
  : public ValueTransferToLibrary<FromType, ToType>
{
public:

  /// Default constructor.
  MatrixValueTransferToLibrary(const unsigned int& from_size, 
                               FromType* from, ToType* to = NULL)
    : ValueTransferToLibrary<FromType, ToType>(from_size, from, to)
  {}


  /// Destructor
  ~MatrixValueTransferToLibrary(void) {}

protected:

  /// Compute the size of the "to" buffer
  inline unsigned int p_computeToSize(void) const;

  /// Copy the value
  inline void p_copy(void);

};

// -------------------------------------------------------------
// MatrixValueTransferToLibrary::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
MatrixValueTransferToLibrary<FromType, ToType>::p_computeToSize(void) const
{
  return ValueTransferToLibrary<FromType, ToType>::p_computeToSize();
}

template <>
inline unsigned int
MatrixValueTransferToLibrary<ComplexType, RealType>::p_computeToSize(void) const
{
  return p_fromSize*TWO*TWO;
}

// -------------------------------------------------------------
// MatrixValueTransferToLibrary::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
MatrixValueTransferToLibrary<FromType, ToType>::p_copy(void)
{
  ValueTransferToLibrary<FromType, ToType>::p_copy();
}

template <>
inline void
MatrixValueTransferToLibrary<ComplexType, RealType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (; fidx < p_fromSize; ++fidx) {
    ComplexType x(p_from[fidx]);
    p_to.get()[tidx] = std::real(x); ++tidx;
    p_to.get()[tidx] = -std::imag(x); ++tidx;
    p_to.get()[tidx] = std::imag(x); ++tidx;
    p_to.get()[tidx] = std::real(x); ++tidx;
  }
}

// -------------------------------------------------------------
//  class MatrixValueTransferFromLibrary
// -------------------------------------------------------------
template <typename FromType, typename ToType>
class MatrixValueTransferFromLibrary 
  : public ValueTransferFromLibrary<FromType, ToType>
{
public:

  /// Default constructor.
  MatrixValueTransferFromLibrary(const unsigned int& from_size, 
                               FromType* from, ToType* to = NULL)
    : ValueTransferFromLibrary<FromType, ToType>(from_size, from, to)
  {}


  /// Destructor
  ~MatrixValueTransferFromLibrary(void) {}

protected:

  /// Compute the size of the "to" buffer
  inline unsigned int p_computeToSize(void) const;

  /// Copy the value
  inline void p_copy(void);

};

// -------------------------------------------------------------
// MatrixValueTransferFromLibrary::p_computeToSize
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline unsigned int
MatrixValueTransferFromLibrary<FromType, ToType>::p_computeToSize(void) const
{
  return ValueTransferFromLibrary<FromType, ToType>::p_computeToSize();
}

template <>
inline unsigned int
MatrixValueTransferFromLibrary<RealType, ComplexType>::p_computeToSize(void) const
{
  return p_fromSize/TWO/TWO;
}

// -------------------------------------------------------------
// MatrixValueTransferFromLibrary::p_copy
// -------------------------------------------------------------
template <typename FromType, typename ToType>
inline void
MatrixValueTransferFromLibrary<FromType, ToType>::p_copy(void)
{
  ValueTransferFromLibrary<FromType, ToType>::p_copy();
}

template <>
inline void
MatrixValueTransferFromLibrary<RealType, ComplexType>::p_copy(void)
{
  unsigned int fidx(0), tidx(0);
  for (; fidx < p_fromSize; fidx += TWO*TWO, ++tidx) {
    ComplexType x(p_from[fidx], -p_from[fidx+1]);
    p_to.get()[tidx] = x;
  }
}

// -------------------------------------------------------------
// storage_size
// -------------------------------------------------------------
/// Get the number of library elements used to represent a scalar
/** 
 * complex -> complex -> 1
 * real -> real -> 1
 * real -> complex -> 1
 * complex -> real -> 2
 * 
 * @return 
 */
template <typename ScalarType, typename LibraryType>
struct storage_size
  : boost::mpl::if_<
      typename boost::mpl::and_< 
        boost::is_same<ScalarType, ComplexType>,
        boost::is_same<LibraryType, RealType> >::type,
      boost::mpl::int_<2>,
    boost::mpl::int_<1> 
  >::type
{};






} // namespace math
} // namespace gridpack



#endif
