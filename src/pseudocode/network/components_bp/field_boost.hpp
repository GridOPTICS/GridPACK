// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   field_boost.hpp
 * @author William A. Perkins
 * @date   Tue Apr 16 09:54:52 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 16, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _field_boost_hpp_
#define _field_boost_hpp_

#include <cassert>
#include <vector>
#include <boost/variant.hpp>

// -------------------------------------------------------------
//  class BaseField
// -------------------------------------------------------------
class BaseField {
public:

  /// The available types of field values 
  enum AvailableFieldTypes {
    BOOL_TYPE,
    INT_TYPE,
    DOUBLE_TYPE,
    STRING_TYPE
  };

  /// Construct w/ known type and size
  BaseField(const AvailableFieldTypes& thetype, const size_t& size)
    : p_type(thetype), p_values(size)
  {}

  /// Destructor
  ~BaseField(void)
  {}
  
  /// Get the type of this instance
  const AvailableFieldTypes& type(void) const
  {
    return p_type;
  }

  /// Insert a value at a specific index
  template <typename T>
  void insert(const int& index, const T& value)
  {
    p_right_type<T>();
    p_values[index] = value;
  }

  /// Get a value at a specific index
  template <typename T>
  T get(const int& index) const
  {
    p_right_type<T>();
    return boost::get<T>(p_values[index]);
  }

  /// Get the size of this field
  size_t size(void) const
  {
    return p_values.size();
  }

protected:

  /// The actual field type (like a C union - sort of)
  /**
   * Must match ::AvailableFieldTypes. 
   * Should be meta-programmed. 
   * 
   */
  typedef boost::variant<bool, int, double, std::string> p_FieldType;

  /// A thing to hold a field
  typedef std::vector<p_FieldType> p_VectorType;

  /// The type of this instance
  /**
   * This may seem redundant, but boost::variant only checks variant
   * types when they are accessed not when they are set.
   * 
   */
  const AvailableFieldTypes p_type;

  /// The field values contained by this instance
  p_VectorType p_values;

  /// A way to make sure this has has the correct type
  template <typename T> void p_right_type(void) const;
  
};

// -------------------------------------------------------------
// BaseField::im_right_type instantiation
// 
// If the user of a BaseField tries to extract an element that is not
// the correct type, it is considered a programming error, even though
// it's a run time decision, and assert is called.
//
// Should be metaprogrammed
// -------------------------------------------------------------

/// Generic instantiation (fails for unknown types)
template <typename T>
void 
BaseField::p_right_type(void) const
{
  assert(false);                // for most types
}

/// Partial instantiation for bool
template <>
void 
BaseField::p_right_type<bool>(void) const
{
  assert(p_type == BOOL_TYPE);
}

/// Partial instantiation for int
template <>
void 
BaseField::p_right_type<int>(void) const
{
  assert(p_type == INT_TYPE);
}

/// Partial instantiation for double
template <>
void 
BaseField::p_right_type<double>(void) const
{
  assert(p_type == DOUBLE_TYPE);
}

/// Partial instantiation for string
template <>
void 
BaseField::p_right_type<std::string>(void) const
{
  assert(p_type == STRING_TYPE);
}

#endif
