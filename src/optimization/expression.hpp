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
 * @date   2015-11-23 11:46:06 d3g096
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

#include <boost/serialization/export.hpp>

#include <gridpack/optimization/variable.hpp>

namespace gridpack {
namespace optimization {

// required forward declarations

template <typename T> class ConstantExpression;
typedef ConstantExpression<int> IntegerConstant;
typedef ConstantExpression<double> RealConstant;

class VariableExpression;

class UnaryExpression;
class UnaryMinus;
class UnaryPlus;

class BinaryExpression;
class Multiplication;
class Division;
class Addition;
class Subtraction;
class Exponentiation;

class Constraint;
class LessThan;
class LessThanOrEqual;
class GreaterThan;
class GreaterThanOrEqual;
class Equal;

// -------------------------------------------------------------
//  class ExpressionVisitor
// -------------------------------------------------------------
/// A cyclic visitor for the Expression class tree
class ExpressionVisitor 
  : private utility::Uncopyable
{
public:

  /// Default constructor.
  ExpressionVisitor(void);

  /// Destructor
  ~ExpressionVisitor(void);

  virtual void visit(IntegerConstant& e);
  virtual void visit(RealConstant& e);
  virtual void visit(VariableExpression& e);

  virtual void visit(UnaryExpression& e);
  virtual void visit(UnaryMinus& e);
  virtual void visit(UnaryPlus& e);

  virtual void visit(BinaryExpression& e);
  virtual void visit(Multiplication& e);
  virtual void visit(Division& e);
  virtual void visit(Addition& e);
  virtual void visit(Subtraction& e);
  virtual void visit(Exponentiation& e);

  virtual void visit(Constraint& e);
  virtual void visit(LessThan& e);
  virtual void visit(LessThanOrEqual& e);
  virtual void visit(GreaterThan& e);
  virtual void visit(GreaterThanOrEqual& e);
  virtual void visit(Equal& e);

};

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

  /// Allow visits from visitors
  void accept(ExpressionVisitor& visitor) 
  {
    this->p_accept(visitor);
  }

protected:

  /// The precedence of this expression
  const int p_precedence;

  /// Make a string representation of this instance (specialized)
  virtual std::string p_render(void) const = 0;

  virtual void p_accept(ExpressionVisitor& e) = 0;

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

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

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

  /// Get the variable
  VariablePtr var(void) const
  {
    return p_var;
  }

  /// Set the variable (be careful now!)
  void var(VariablePtr v) 
  {
    p_var = v;
  }

protected:
  
  /// The variable in question
  VariablePtr p_var;

  std::string p_render(void) const
  {
    return p_var->name();
  }
  
  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
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
//  class UnaryExpression
// -------------------------------------------------------------
class UnaryExpression 
  : public Expression
{
public:

  /// Default constructor.
  UnaryExpression(const int& prec, const std::string& op, 
                  ExpressionPtr expr)
    : Expression(prec), p_operator(op), p_expr(expr)
  {}

  /// Destructor
  ~UnaryExpression(void)
  {}

  /// Get the RHS of the expresion
  ExpressionPtr rhs()
  {
    return p_expr;
  }

protected:
  
  /// The operator used for this instance
  const std::string p_operator;

  /// The expression
  ExpressionPtr p_expr;

  std::string p_render(void) const
  {
    std::string s(this->p_operator);
    if (p_expr->precedence() > this->precedence()) {
      s += "[" + p_expr->render() + "]";
    } else {
      s += p_expr->render();
    }
    return s;
  }
  
  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  UnaryExpression(void)
    : Expression(0)
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Expression>(*this);
    ar & const_cast<std::string&>(p_operator) & p_expr;
  }
};


// -------------------------------------------------------------
//  class UnaryMinus
// -------------------------------------------------------------
class UnaryMinus 
  : public UnaryExpression
{
public:

  /// Default constructor.
  UnaryMinus(ExpressionPtr expr)
    : UnaryExpression(3, "-", expr)
  {}

  /// Protected copy constructor to avoid unwanted copies.
  UnaryMinus(const UnaryMinus& old)
    : UnaryExpression(old)
  {}

  /// Destructor
  ~UnaryMinus(void) {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  UnaryMinus(void)
    : UnaryExpression()
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<UnaryExpression>(*this);
  }
};

// -------------------------------------------------------------
//  class UnaryPlus
// -------------------------------------------------------------
class UnaryPlus 
  : public UnaryExpression
{
public:

  /// Default constructor.
  UnaryPlus(ExpressionPtr expr)
    : UnaryExpression(3, "+", expr)
  {}

  /// Protected copy constructor to avoid unwanted copies.
  UnaryPlus(const UnaryPlus& old)
    : UnaryExpression(old)
  {}

  /// Destructor
  ~UnaryPlus(void) {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  UnaryPlus(void)
    : UnaryExpression()
  {}

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<UnaryExpression>(*this);
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

  const std::string& op(void) const
  {
    return p_operator;
  }

  /// Get the LHS of the expression
  ExpressionPtr lhs(void) 
  {
    return p_LHS;
  }

  
  /// Get the RHS of the expresion
  ExpressionPtr rhs()
  {
    return p_RHS;
  }

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
    if (p_LHS) {
      if (p_LHS->precedence() > this->precedence()) {
        s += "[" + p_LHS->render() + "]";
      } else {
        s += p_LHS->render();
      }
    } else {
      s += "(empty)";
    }
    s += " " + this->p_operator + " ";
    if (p_RHS) {
      if (p_RHS->precedence() > this->precedence()) {
        s += "[" + p_RHS->render() + "]";
      } else {
        s += p_RHS->render();
      }
    } else {
      s += "(empty)";
    }
    return s;
  }

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Multiplication(void) 
    : BinaryExpression()
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

inline
ExpressionPtr operator*(VariablePtr lhs, VariablePtr rhs)
{
  ExpressionPtr lv(new VariableExpression(lhs));
  ExpressionPtr rv(new VariableExpression(rhs));
  return lv * rv;
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
//  class Division
// -------------------------------------------------------------
class Division 
  : public BinaryExpression
{
public:

  /// Default constructor.
  Division(ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(5, "/", lhs, rhs)
  {}

  /// Destructor
  ~Division(void)
  {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Division(void) 
    : BinaryExpression()
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
ExpressionPtr operator/(ExpressionPtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr result(new Division(lhs, rhs));
  return result;
}

template <typename T>
ExpressionPtr operator/(T lhs, ExpressionPtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  return c / rhs;
}

template <typename T>
ExpressionPtr operator/(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return lhs / c;
}

inline
ExpressionPtr operator/(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v / rhs;
}

inline
ExpressionPtr operator/(ExpressionPtr lhs, VariablePtr rhs)
{
  ExpressionPtr v(new VariableExpression(rhs));
  return lhs / v;
}

inline
ExpressionPtr operator/(VariablePtr lhs, VariablePtr rhs)
{
  ExpressionPtr lv(new VariableExpression(lhs));
  ExpressionPtr rv(new VariableExpression(rhs));
  return lv / rv;
}


// Do not allow a nonlinear expression
template <typename T>
ExpressionPtr operator/(T lhs, VariablePtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  ExpressionPtr v(new VariableExpression(rhs));
  return c / v;
}

template <typename T>
ExpressionPtr operator/(VariablePtr lhs, T rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return v / c;
}

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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Addition(void) 
    : BinaryExpression()
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
//  class Subtraction
// -------------------------------------------------------------
class Subtraction 
  : public BinaryExpression
{
public:

  /// Default constructor.
  Subtraction(ExpressionPtr lhs, ExpressionPtr rhs)
    : BinaryExpression(6, "-", lhs, rhs)
  {}

  /// Copy constructor
  Subtraction(const Subtraction& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Subtraction(void)
  {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Subtraction(void) 
    : BinaryExpression()
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
ExpressionPtr operator+(ExpressionPtr expr)
{
  return ExpressionPtr(new UnaryPlus(expr));
}

inline
ExpressionPtr operator+(VariablePtr var)
{
  ExpressionPtr expr(new VariableExpression(var));
  return +expr;
}
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

inline
ExpressionPtr operator+(VariablePtr lhs, VariablePtr rhs)
{
  ExpressionPtr lv(new VariableExpression(lhs));
  ExpressionPtr rv(new VariableExpression(rhs));
  return lv + rv;
}

// -------------------------------------------------------------
// operator-
// -------------------------------------------------------------
inline
ExpressionPtr operator-(ExpressionPtr expr)
{
  return ExpressionPtr(new UnaryMinus(expr));
}

inline
ExpressionPtr operator-(VariablePtr var)
{
  ExpressionPtr expr(new VariableExpression(var));
  return -expr;
}

inline
ExpressionPtr operator-(ExpressionPtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr result(new Subtraction(lhs, rhs));
  return result;
}

template <typename T>
ExpressionPtr operator-(T lhs, ExpressionPtr rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(lhs));
  return c - rhs;
}

template <typename T>
ExpressionPtr operator-(ExpressionPtr lhs, T rhs)
{
  ExpressionPtr c(new ConstantExpression<T>(rhs));
  return lhs - c;
}

inline
ExpressionPtr operator-(VariablePtr lhs, ExpressionPtr rhs)
{
  ExpressionPtr v(new VariableExpression(lhs));
  return v - rhs;
}

inline
ExpressionPtr operator-(ExpressionPtr lhs, VariablePtr rhs)
{
  ExpressionPtr v(new VariableExpression(rhs));
  return lhs - v;
}

inline
ExpressionPtr operator-(VariablePtr lhs, VariablePtr rhs)
{
  ExpressionPtr lv(new VariableExpression(lhs));
  ExpressionPtr rv(new VariableExpression(rhs));
  return lv - rv;
}

// -------------------------------------------------------------
//  class Exponentiation
// -------------------------------------------------------------
class Exponentiation 
  : public BinaryExpression
{
protected:

public:

  /// Default constructor.
  Exponentiation(ExpressionPtr lhs, int exp)
    : BinaryExpression(2, "^", lhs, ExpressionPtr(new IntegerConstant(exp)))
  {}

  /// Protected copy constructor to avoid unwanted copies.
  Exponentiation(const Exponentiation& old)
    : BinaryExpression(old)
  {}

  /// Destructor
  ~Exponentiation(void)
  {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Exponentiation(void) 
    : BinaryExpression()
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
// operator^
// -------------------------------------------------------------
inline
ExpressionPtr operator^(ExpressionPtr expr, int exp)
{
  ExpressionPtr result(new Exponentiation(expr, exp));
  return result;
}

inline
ExpressionPtr operator^(VariablePtr v, int exp)
{
  ExpressionPtr vexpr(new VariableExpression(v));
  ExpressionPtr result(new Exponentiation(vexpr, exp));
  return result;
}


// -------------------------------------------------------------
//  class Constraint
// -------------------------------------------------------------
class Constraint 
  : public BinaryExpression,
    public utility::Named
{
public:

  /// Default constructor.
  Constraint(const int& prec, const std::string& op, 
             ExpressionPtr lhs, ExpressionPtr rhs);

  /// Destructor
  ~Constraint(void)
  {}

  /// Add something to the LHS
  void addToLHS(ExpressionPtr e)
  {
    if (e) {
      if (p_LHS) {
        p_LHS = e + p_LHS;
      } else {
        p_LHS = e;
      }
    }
  }

  // /// Do whatever 
  // void evaluate(void) const
  // {
  //   BinaryExpression::evaluate();
  // }
  
  // int precedence() const
  // {
  //   return BinaryExpression::precedence();
  // }

  // std::string render(void) const
  // {
  //   return BinaryExpression::render();
  // }

  // const std::string& op(void) const
  // {
  //   return BinaryExpression::op();
  // }

  // /// Get the LHS of the expression
  // ExpressionPtr lhs(void) 
  // {
  //   return BinaryExpression::lhs();
  // }
  
  // /// Get the RHS of the expresion
  // ExpressionPtr rhs()
  // {
  //   return BinaryExpression::rhs();
  // }
  
  /// Allow visits from visitors
  void accept(ExpressionVisitor& visitor) 
  {
    BinaryExpression::accept(visitor);
  }

protected:

  static int p_nextID;

  /// Constructor for serialization
  Constraint(void);

private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<BinaryExpression>(*this);
    ar & boost::serialization::base_object<utility::Named>(*this);
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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  LessThan(void) 
    : Constraint()
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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  LessThanOrEqual(void) 
    : Constraint()
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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  GreaterThan(void) 
    : Constraint()
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

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  GreaterThanOrEqual(void) 
    : Constraint()
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
    : Constraint(9, "=", lhs, rhs)
  {}

  /// Destructor
  ~Equal(void)
  {}

protected:

  void p_accept(ExpressionVisitor& e)
  {
    e.visit(*this);
  }

  /// Constructor for serialization
  Equal(void) 
    : Constraint()
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

// -------------------------------------------------------------
//  class ExpressionChecker
// -------------------------------------------------------------
class ExpressionChecker 
  : public ExpressionVisitor
{
public:

  /// The visited expression is a constant
  bool isConstant;

  /// The visited expression is an integer constant
  bool isInteger;

  /// The visited expression is a variable
  bool isVariable;

  /// The visited expression is exponentiation
  bool isExponentiation;

  /// Default constructor.
  ExpressionChecker(void)
    : ExpressionVisitor(), 
      isConstant(false), isInteger(false), isVariable(false), isExponentiation(false)
  {}

  /// Destructor
  ~ExpressionChecker(void)
  {}

  void visit(IntegerConstant& e)
  {
    isConstant = true;
    isInteger = true;
  }
  void visit(RealConstant& e)
  {
    isConstant = true;
  }
  void visit(VariableExpression& e)
  {
    isVariable = true;
  }
  void visit(Exponentiation& e)
  {
    isExponentiation = true;
  }
};

} // namespace optimization
} // namespace gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::ConstantExpression<int>);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::ConstantExpression<double>);

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::VariableExpression);

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::UnaryMinus);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::UnaryPlus);

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Multiplication);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Division);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Addition);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Subtraction);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Exponentiation);

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Constraint);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::LessThan);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::LessThanOrEqual);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::GreaterThan);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::GreaterThanOrEqual);
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Equal);


#endif
