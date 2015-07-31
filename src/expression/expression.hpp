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
 * @file   expression.hpp
 * @author William A. Perkins
 * @date   2015-07-31 13:43:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _expression_hpp_
#define _expression_hpp_

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>

#include "variable.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class Expression
// -------------------------------------------------------------
class Expression {
public:

  /// Default constructor.
  Expression(void) {}

  /// Destructor
  virtual ~Expression(void) {}

  /// Do whatever 
  void evaluate(void) const
  {
    this->p_evaluate();
  }

protected:

  virtual void p_evaluate(void) const = 0; 

};

typedef boost::shared_ptr<Expression> ExpressionPtr;

// -------------------------------------------------------------
//  class ConstantExpression
// -------------------------------------------------------------
template <typename T>
class ConstantExpression 
  : public Expression
{
public:

  /// Default constructor.
  ConstantExpression(const T& value)
    : Expression(), p_value(value)
  {}

  /// Destructor
  ~ConstantExpression(void)
  {}

  /// Get the parameter value
  T value(void) const
  {
    return p_value;
  }

protected:

  const T p_value;              /**< constant value */

  void p_evaluate(void) const
  {
    std::cout << p_value;
  }

  /// Constructor for serialization
  ConstantExpression(void) : p_value() {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<Expression>(*this);
    ar & const_cast<T>(p_value);
  }
};

typedef ConstantExpression<int> IntegerConstant;
typedef ConstantExpression<double> RealConstant;

// -------------------------------------------------------------
//  class VariableExpression
// -------------------------------------------------------------
class VariableExpression
  : public Expression
{
public:

  /// Default constructor.
  VariableExpression(VariablePtr v)
    : Expression(), p_var(v)
  {}

  /// Destructor
  ~VariableExpression(void) 
  {}

protected:
  
  /// The variable in question
  VariablePtr p_var;

  void p_evaluate(void) const
  {
    std::cout << p_var->name();
  }
  
private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<Expression>(*this);
    ar & p_var;
  }
};


// -------------------------------------------------------------
//  class BinaryExpression
// -------------------------------------------------------------
class BinaryExpression 
  : public Expression
{
public:

  /// Default constructor.
  BinaryExpression(ExpressionPtr lhs, ExpressionPtr rhs)
    : Expression(), p_LHS(lhs), p_RHS(rhs)
  {}

  /// Destructor
  ~BinaryExpression(void)
  {}

protected:

  /// The left hand side expression
  ExpressionPtr p_LHS;

  /// The right hand side expression
  ExpressionPtr p_RHS;

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Expression>(*this);
    ar & p_LHS & p_RHS;
  }
};


// -------------------------------------------------------------
//  class Addition
// -------------------------------------------------------------
class Addition 
  : public BinaryExpression
{
public:

  /// Default constructor.
  Addition(ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(lhs, rhs)
  {}

  /// Destructor
  ~Addition(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " + ";
    p_RHS->evaluate();
  }


private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<BinaryExpression>(*this);
  }

};

// -------------------------------------------------------------
// operator+
// -------------------------------------------------------------
ExpressionPtr operator+(ExpressionPtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr result(new Addition(lhs, rhs));
  return result;
}

template <typename T>
ExpressionPtr operator+(T lhs, ExpressionPtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  return c + rhs;
}

template <typename T>
ExpressionPtr operator+(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return lhs + c;
}

ExpressionPtr operator+(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v + rhs;
}

ExpressionPtr operator+(ExpressionPtr lhs, VariablePtr rhs)
{
  ExpressionPtr v(new VariableExpression(rhs));
  return lhs + v;
}


// -------------------------------------------------------------
//  class Multiplication
// -------------------------------------------------------------
class Multiplication 
  : public BinaryExpression
{
public:

  /// Default constructor.
  Multiplication(ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(lhs, rhs)
  {}

  /// Destructor
  ~Multiplication(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " * ";
    p_RHS->evaluate();
  }

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<BinaryExpression>(*this);
  }

};

// -------------------------------------------------------------
// operator*
// -------------------------------------------------------------
ExpressionPtr operator*(ExpressionPtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr result(new Multiplication(lhs, rhs));
  return result;
}

template <typename T>
ExpressionPtr operator*(T lhs, ExpressionPtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  return c * rhs;
}

template <typename T>
ExpressionPtr operator*(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return lhs * c;
}

ExpressionPtr operator*(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v * rhs;
}

ExpressionPtr operator*(ExpressionPtr lhs, VariablePtr rhs)
{
  ExpressionPtr v(new VariableExpression(rhs));
  return lhs * v;
}

template <typename T>
ExpressionPtr operator*(T lhs, VariablePtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  ExpressionPtr v(new VariableExpression(rhs));
  return c * v;
}

template <typename T>
ExpressionPtr operator*(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return v * c;
}



} // namespace optimization
} // namespace gridpack

#endif
