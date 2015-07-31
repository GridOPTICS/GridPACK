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
 * @file   variable.hpp
 * @author William A. Perkins
 * @date   2015-07-30 15:26:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _variable_hpp_
#define _variable_hpp_


#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <gridpack/utilities/named.hpp>
#include <gridpack/utilities/uncopyable.hpp>

namespace gridpack {
namespace optimization {

// forward declarations
class Variable;

template <typename T> class BoundedVariableT;
typedef BoundedVariableT<int> IntegerVariable;
typedef BoundedVariableT<double> RealVariable;

class BinaryVariable;

// -------------------------------------------------------------
//  class VariableVisitor
// -------------------------------------------------------------
class VariableVisitor 
  : private utility::Uncopyable
{
public:

  /// Default constructor.
  VariableVisitor(void);

  /// Destructor
  virtual ~VariableVisitor(void);

  virtual void visit(const Variable& var);
  virtual void visit(const RealVariable& var);
  virtual void visit(const IntegerVariable& var);
  virtual void visit(const BinaryVariable& var);
};



// -------------------------------------------------------------
//  class Variable
// -------------------------------------------------------------
class Variable 
  : public utility::Named,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  Variable(void);

  /// Destructor
  virtual ~Variable(void);

  /// Allow visits from visitors
  void accept(VariableVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

protected:

  static int p_nextID;

  int p_id;                     /**< The (local) unique identifer */

  /// Allow visitor from variable visitors (specialized)
  virtual void p_accept(VariableVisitor& visitor) const
  {
    visitor.visit(*this);
  }


private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<utility::Named>(*this);
    ar & p_id;
  }

};

typedef boost::shared_ptr<Variable> VariablePtr;

// -------------------------------------------------------------
//  class BoundedVariableT
// -------------------------------------------------------------
/// Represents a variable that may or may not have bounds
template <typename T>
class BoundedVariableT 
  : public Variable
{
public:

  /// Construct with an initial value
  BoundedVariableT(const T& value)
    : Variable(),
      p_initial(value), 
      p_lowBound(p_veryLowValue),
      p_highBound(p_veryHighValue)
  {}

  /// Construct with an initial value and bounds
  BoundedVariableT(const T& value, const T& lo, const T& hi)
    : Variable(),
      p_initial(value), 
      p_lowBound(lo),
      p_highBound(hi)
  {}

  /// Destructor
  ~BoundedVariableT(void)
  {}

protected:

  static const T p_veryLowValue; /**< lower bound for unbounded variables */
  static const T p_veryHighValue; /**< upper bound for unbounded variables */

  T p_initial;                  /**< initial value */
  T p_lowBound;                 /**< lower bound */
  T p_highBound;                /**< upper bound */

  /// Constructor used for serialization
  BoundedVariableT(void)
    : Variable(),
      p_initial(), 
      p_lowBound(p_veryLowValue),
      p_highBound(p_veryHighValue)
  {}

  /// Allow visitor from variable visitors (specialized)
  void p_accept(VariableVisitor& visitor) const
  {
    visitor.visit(*this);
  }

private:

  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<Variable>(*this);
    ar & p_initial & p_lowBound & p_highBound;
  }

};



// -------------------------------------------------------------
//  class BinaryVariable
// -------------------------------------------------------------
/// Represents a variable that can be only 1 or 0
class BinaryVariable 
  : public BoundedVariableT<int>
{
public:

  /// Default constructor.
  BinaryVariable(const int& i)
    : BoundedVariableT<int>(i)
  {
    BoundedVariableT<int>::p_highBound = 1;
    BoundedVariableT<int>::p_lowBound = 0;
  }

  /// Destructor
  ~BinaryVariable(void)
  {}

protected:

  /// Allow visitor from variable visitors (specialized)
  void p_accept(VariableVisitor& visitor) const
  {
    visitor.visit(*this);
  }

private:

  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<IntegerVariable>(*this);
  }
};

} // namespace optimization
} // namespace gridpack

#endif
