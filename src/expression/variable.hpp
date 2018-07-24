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
 * @date   2017-02-10 07:32:27 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _variable_hpp_
#define _variable_hpp_

#include <iosfwd>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/named.hpp>

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

  virtual void visit(Variable& var);
  virtual void visit(RealVariable& var);
  virtual void visit(IntegerVariable& var);
  virtual void visit(BinaryVariable& var);
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

  /// Get this variable's id
  int id(void) const
  {
    return p_id;
  }

  /// Set this variable's id
  void id(const int& id) 
  {
    p_id = id;
  }
  
  /// reset the variable's id
  void clear()
  {
    p_nextID = 0;
  }

  /// Allow visits from visitors
  void accept(VariableVisitor& visitor) 
  {
    this->p_accept(visitor);
  }

  /// Set internal no init flag. This turn off variable declaration
  /// statement in Julia code generator.
  void setNoInit(bool flag)
  {
    p_no_init = flag;
  }

  /// return state of no init flag
  bool getNoInit()
  {
    return p_no_init;
  }

protected:

  static int p_nextID;

  int p_id;                     /**< The (local) unique identifer */

  /// Allow visitor from variable visitors (specialized)
  virtual void p_accept(VariableVisitor& visitor)
  {
    visitor.visit(*this);
  }

private:
  
  friend class boost::serialization::access;
  
  /// If set to true, do not export a variable declaration statement
  /// (at least when exporting Julia code). This can be used to
  /// to create variables that are set in an auxiliary file.
  bool p_no_init;

  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<utility::Named>(*this);
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
      p_lowBound(veryLowValue),
      p_highBound(veryHighValue)
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

  /// Get the initial value
  T initial(void) const
  {
    return p_initial;
  }

  /// Set the initial value
  void initial(const T& value) 
  {
    p_initial = value;
  }

  /// Is this variable bounded
  bool bounded(void) const
  { 
    return ((p_lowBound > veryLowValue) ||
            (p_highBound < veryHighValue));
  }

  T lowerBound(void) const
  {
    return p_lowBound;
  }
  T upperBound(void) const
  {
    return p_highBound;
  }

  /// lower bound for unbounded variables
  static const T veryLowValue; 

  /// upper bound for unbounded variables
  static const T veryHighValue; 

protected:

  T p_initial;                  /**< initial value */
  T p_lowBound;                 /**< lower bound */
  T p_highBound;                /**< upper bound */

  /// Constructor used for serialization
  BoundedVariableT(void)
    : Variable(),
      p_initial(), 
      p_lowBound(veryLowValue),
      p_highBound(veryHighValue)
  {}

  /// Allow visitor from variable visitors (specialized)
  void p_accept(VariableVisitor& visitor) 
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
  void p_accept(VariableVisitor& visitor)
  {
    visitor.visit(*this);
  }

  /// Constructor for serialization
  BinaryVariable(void)
    : BoundedVariableT<int>(0)
  {
    BoundedVariableT<int>::p_highBound = 1;
    BoundedVariableT<int>::p_lowBound = 0;
  }

private:

  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar  & boost::serialization::base_object<IntegerVariable>(*this);
  }
};

// -------------------------------------------------------------
//  class SetVariableInitial
// -------------------------------------------------------------
class SetVariableInitial 
  : public VariableVisitor
{
public:

  /// Default constructor.
  SetVariableInitial(const double& v)
    : VariableVisitor(), p_value(v)
  {}

  /// Destructor
  ~SetVariableInitial(void)
  {}

  void visit(RealVariable& var) 
  {
    var.initial(p_value);
  }
  void visit(IntegerVariable& var)
  {
    var.initial(static_cast<int>(p_value));
  }

protected:
  
  double p_value;
};

// -------------------------------------------------------------
//  class GetVariableInitial
// -------------------------------------------------------------
class GetVariableInitial 
  : public VariableVisitor
{
public:

  /// Default constructor.
  GetVariableInitial(void)
    : VariableVisitor(), p_value()
  {}

  /// Destructor
  ~GetVariableInitial(void)
  {}

  void visit(RealVariable& var) 
  {
    p_value = var.initial();
  }
  void visit(IntegerVariable& var)
  {
    p_value = static_cast<double>(var.initial());
  }

  double value(void) const
  {
    return p_value;
  }

protected:
  
  double p_value;
};

// -------------------------------------------------------------
//  class VariableTable
// -------------------------------------------------------------
class VariableTable 
  : public VariableVisitor
{
public:

  /// Default constructor.
  VariableTable(std::ostream& out);

  /// Destructor
  ~VariableTable(void);

  void visit(RealVariable& var);
  void visit(IntegerVariable& var);
  void visit(BinaryVariable& var);

protected:

  std::ostream& p_out;

  bool p_first;

  void p_header(void) const;

};


// -------------------------------------------------------------
//  class VariableCounter
// -------------------------------------------------------------
class VariableCounter 
  : public VariableVisitor
{
public:

  /// Default constructor.
  VariableCounter(void);

  /// Destructor
  ~VariableCounter(void);

  int numVar;                   /**< total number of variables */
  int numReal;                  /**< number of real variables */
  int numInt;                   /**< number of integer (not binary) variables */
  int numBin;                   /**< number of binary variables */

  void visit(Variable& var);
  void visit(RealVariable& var);
  void visit(IntegerVariable& var);
  void visit(BinaryVariable& var);

};

} // namespace optimization
} // namespace gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Variable)
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::BoundedVariableT<double>)
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::BoundedVariableT<int>)
BOOST_CLASS_EXPORT_KEY(gridpack::optimization::BinaryVariable)


#endif
