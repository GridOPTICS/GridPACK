// -------------------------------------------------------------
// file: field_boost_test.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 16, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include "field_boost.hpp"


// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  BaseField fld(BaseField::INT_TYPE, 3);

  fld.insert<int>(0, 1);
  fld.insert<int>(1, 2);
  fld.insert<int>(2, 3);

  int i, j, k;

  i = fld.get<int>(0);
  j = fld.get<int>(1);
  k = fld.get<int>(2);

  std::cout << i << ", " << j << ", " << k << std::endl;

  std::string s;

  // this will cause an assertion, if enabled, otherwise a
  // boost::bad_get exception is thrown

  s = fld.get<std::string>(2);

  return 0;
}
