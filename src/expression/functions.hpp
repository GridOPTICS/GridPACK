// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   functions.hpp
 * @author William A. Perkins
 * @date   2017-03-22 07:56:54 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 31, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _functions_hpp_
#define _functions_hpp_

#include <gridpack/expression/expression.hpp>

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class Function
// -------------------------------------------------------------
class Function 
  : public Expression
{
public:

  /// Destructor
  ~Function(void);

  /// Get the function name
  const std::string& name(void) const
  {
    return p_fname;
  }

  /// Get the function arguments
  const std::vector<ExpressionPtr>& args(void) const
  {
    return p_args;
  }

protected:

  /// Default constructor (only friends can construct).
  Function(const std::string& name, 
           const std::vector<ExpressionPtr>& args);

  /// Constructor for serialization
  Function(void); 

  /// The operator used for this instance
  std::string p_fname;

  std::vector<ExpressionPtr> p_args;

  /// Render the function and its arguments
  std::string p_render(void) const;
         
  /// Let a visitor in
  void p_accept(ExpressionVisitor& e);

  /// Is this expression empty?
  bool p_null(void) const;
  
private:
  
  friend class boost::serialization::access;
  
  template<class Archive> 
  void serialize(Archive &ar, const unsigned int)
  {
    ar & boost::serialization::base_object<Expression>(*this);
    ar & p_fname & p_args;
  }

  /// Representation of the sine function
  friend ExpressionPtr sin(ExpressionPtr e);
  
  /// Representation of the cosine function
  friend ExpressionPtr cos(ExpressionPtr e);

  /// Representation of the sine function
  friend ExpressionPtr sin(VariablePtr e);
  
  /// Representation of the cosine function
  friend ExpressionPtr cos(VariablePtr e);

};

extern ExpressionPtr sin(ExpressionPtr e);
extern ExpressionPtr cos(ExpressionPtr e);
extern ExpressionPtr sin(VariablePtr e);
extern ExpressionPtr cos(VariablePtr e);

typedef boost::shared_ptr<Function> FunctionPtr;

} // namespace optimization
} // namespace gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::optimization::Function)


#endif
