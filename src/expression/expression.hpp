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
 * @date   2015-08-28 16:45:35 d3g096
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
#include <boost/format.hpp>
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
  Expression(int prec)
    : p_precedence(prec)
  {}

  /// Copy constructor
  Expression(const Expression& old) 
    : p_precedence(old.p_precedence)
  {}

  /// Destructor
  virtual ~Expression(void) {}

  /// What is the precedence of this expression
  int precedence(void) const
  {
    return p_precedence;
  }

  /// Do whatever 
  void evaluate(void) const
  {
    std::cout << this->p_render() << std::endl;
  }

  /// Make a string representation of this instance
  std::string render(void) const
  {
    return this->p_render();
  }

protected:

  /// The precedence of this expression
  const int p_precedence;

  /// Make a string representation of this instance (specialized)
  virtual std::string p_render(void) const = 0;

  /// Constructor for serialization
  Expression(void) : p_precedence() {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & const_cast<int&>(p_precedence);
  }
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
    : Expression(0), p_value(value)
  {}

  /// Copy constructor
  ConstantExpression(const ConstantExpression& old)
    : Expression(old), p_value(old.p_value)
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

  std::string p_render(void) const;

  /// Constructor for serialization
  ConstantExpression(void) 
    : Expression(), p_value() 
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<Expression>(*this);
    ar & const_cast<T&>(p_value);
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
    : Expression(0), p_var(v)
  {}

  /// Copy constructor
  VariableExpression(const VariableExpression& old)
    : Expression(old), p_var(old.p_var)
  {}

  /// Destructor
  ~VariableExpression(void) 
  {}

  /// Get the variable name
  std::string name(void) const
  {
    return p_var->name();
  }

protected:
  
  /// The variable in question
  VariablePtr p_var;

  std::string p_render(void) const
  {
    return p_var->name();
  }
  
  /// Constructor for serialization
  VariableExpression(void) 
    : Expression(), p_var() 
  {}

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
  BinaryExpression(const int& prec, const std::string& op, 
                   ExpressionPtr lhs, ExpressionPtr rhs)
    : Expression(prec), p_operator(op), p_LHS(lhs), p_RHS(rhs)
  {}

  /// Copy constructor
  BinaryExpression(const BinaryExpression& old)
    : Expression(old), 
      p_operator(old.p_operator),
      p_LHS(old.p_LHS), p_RHS(old.p_RHS)
  {}

  /// Destructor
  ~BinaryExpression(void)
  {}

protected:

  /// The operator used for this instance
  const std::string p_operator;

  /// The left hand side expression
  ExpressionPtr p_LHS;

  /// The right hand side expression
  ExpressionPtr p_RHS;

  /// Constructor for serialization
  BinaryExpression(void) 
    : Expression()
  {}

  std::string p_render(void) const
  {
    std::string s("");
    if (p_LHS->precedence() > this->precedence()) {
      s += "( " + p_LHS->render() + ") ";
    } else {
      s += p_LHS->render();
    }
    s += " " + this->p_operator + " ";
    if (p_RHS->precedence() > this->precedence()) {
      s += "( " + p_RHS->render() + ") ";
    } else {
      s += p_RHS->render();
    }
    return s;
  }

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Expression>(*this);
    ar & const_cast<std::string&>(p_operator) & p_LHS & p_RHS;
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
    : BinaryExpression(6, "+", lhs, rhs)
  {}

  /// Copy constructor
  Addition(const Addition& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Addition(void)
  {}

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
inline
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

inline
ExpressionPtr operator+(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v + rhs;
}

inline
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
    : BinaryExpression(5, "*", lhs, rhs)
  {}

  /// Copy constructor
  Multiplication(const Multiplication& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Multiplication(void)
  {}

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
inline
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

inline
ExpressionPtr operator*(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v * rhs;
}

inline
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



// -------------------------------------------------------------
//  class Constraint
// -------------------------------------------------------------
class Constraint 
  : protected BinaryExpression
{
public:

  /// Default constructor.
  Constraint(const int& prec, const std::string& op, 
             ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(prec, op, lhs, rhs)
  {}

  /// Destructor
  ~Constraint(void)
  {}

  /// Do whatever 
  void evaluate(void) const
  {
    BinaryExpression::evaluate();
  }

  std::string render(void) const
  {
    return BinaryExpression::render();
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
//  class LessThan
// -------------------------------------------------------------
class LessThan 
  : public Constraint
{
public:

  /// Default constructor.
  LessThan(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(8, "<", lhs, rhs)
  {}

  /// Destructor
  ~LessThan(void)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Constraint>(*this);
  }
};

// -------------------------------------------------------------
//  class LessThanOrEqual
// -------------------------------------------------------------
class LessThanOrEqual 
  : public Constraint
{
public:

  /// Default constructor.
  LessThanOrEqual(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(8, "<=", lhs, rhs)
  {}

  /// Destructor
  ~LessThanOrEqual(void)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Constraint>(*this);
  }
};


// -------------------------------------------------------------
//  class GreaterThan
// -------------------------------------------------------------
class GreaterThan 
  : public Constraint
{
public:

  /// Default constructor.
  GreaterThan(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(8, ">", lhs, rhs)
  {}

  /// Destructor
  ~GreaterThan(void)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Constraint>(*this);
  }
};

// -------------------------------------------------------------
//  class GreaterThanOrEqual
// -------------------------------------------------------------
class GreaterThanOrEqual 
  : public Constraint
{
public:

  /// Default constructor.
  GreaterThanOrEqual(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(8, ">=", lhs, rhs)
  {}

  /// Destructor
  ~GreaterThanOrEqual(void)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Constraint>(*this);
  }
};

// -------------------------------------------------------------
// class Equal
// -------------------------------------------------------------
class Equal 
  : public Constraint
{
public:

  /// Default constructor.
  Equal(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(9, "==", lhs, rhs)
  {}

  /// Destructor
  ~Equal(void)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Constraint>(*this);
  }
};



typedef boost::shared_ptr<Constraint> ConstraintPtr;


// -------------------------------------------------------------
// Constraint operator<
// -------------------------------------------------------------
template <typename T>
ConstraintPtr operator<(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new LessThan(lhs, c));
  return result;
}

template <typename T>
ConstraintPtr operator<(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new LessThan(v, c));
  return result;
}

// -------------------------------------------------------------
// Constraint operator<=
// -------------------------------------------------------------
template <typename T>
ConstraintPtr operator<=(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new LessThanOrEqual(lhs, c));
  return result;
}

template <typename T>
ConstraintPtr operator<=(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new LessThanOrEqual(v, c));
  return result;
}

// -------------------------------------------------------------
// Constraint operator>
// -------------------------------------------------------------
template <typename T>
ConstraintPtr operator>(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new GreaterThan(lhs, c));
  return result;
}

template <typename T>
ConstraintPtr operator>(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new GreaterThan(v, c));
  return result;
}
// -------------------------------------------------------------
// Constraint operator>=
// -------------------------------------------------------------
template <typename T>
ConstraintPtr operator>=(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new GreaterThanOrEqual(lhs, c));
  return result;
}

template <typename T>
ConstraintPtr operator>=(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new GreaterThanOrEqual(v, c));
  return result;
}

// -------------------------------------------------------------
// Constraint operator==
// -------------------------------------------------------------
template <typename T>
ConstraintPtr operator==(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new Equal(lhs, c));
  return result;
}

template <typename T>
ConstraintPtr operator==(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  ConstraintPtr result(new Equal(v, c));
  return result;
}


} // namespace optimization
} // namespace gridpack

#endif
