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
 * @date   2015-08-28 09:03:56 d3g096
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

  /// Copy constructor
  Expression(const Expression& old) {}

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

  /// Copy constructor
  BinaryExpression(const BinaryExpression& old)
    : Expression(old), p_LHS(old.p_LHS), p_RHS(old.p_RHS)
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

  /// Copy constructor
  Addition(const Addition& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Addition(void)
  {}

protected:

  void p_evaluate(void) const
  {
    std::cout << "( ";
    p_LHS->evaluate();
    std::cout << " + ";
    p_RHS->evaluate();
    std::cout << " )";
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
    : BinaryExpression(lhs, rhs)
  {}

  /// Copy constructor
  Multiplication(const Multiplication& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Multiplication(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << "*";
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
  : public BinaryExpression
{
public:

  /// Default constructor.
  Constraint(ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(lhs, rhs)
  {}

  /// Destructor
  ~Constraint(void)
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
//  class LessThan
// -------------------------------------------------------------
class LessThan 
  : public Constraint
{
public:

  /// Default constructor.
  LessThan(ExpressionPtr lhs, ExpressionPtr rhs)
    : Constraint(lhs, rhs)
  {}

  /// Destructor
  ~LessThan(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " < ";
    p_RHS->evaluate();
  }


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
    : Constraint(lhs, rhs)
  {}

  /// Destructor
  ~LessThanOrEqual(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " <= ";
    p_RHS->evaluate();
  }

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
    : Constraint(lhs, rhs)
  {}

  /// Destructor
  ~GreaterThan(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " > ";
    p_RHS->evaluate();
  }

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
    : Constraint(lhs, rhs)
  {}

  /// Destructor
  ~GreaterThanOrEqual(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " >= ";
    p_RHS->evaluate();
  }

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
    : Constraint(lhs, rhs)
  {}

  /// Destructor
  ~Equal(void)
  {}

protected:

  void p_evaluate(void) const
  {
    p_LHS->evaluate();
    std::cout << " == ";
    p_RHS->evaluate();
  }

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
