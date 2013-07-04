/**
 * @file   component_serialization.cpp
 * @author William A. Perkins
 * @date   2013-07-04 09:26:38 d3g096
 * 
 * @brief  Serialization tests for various network component classes
 * 
 * 
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/mpi/collectives.hpp>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "data_collection.hpp"

#include <boost/serialization/scoped_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

typedef boost::archive::binary_iarchive inarchive;
typedef boost::archive::binary_oarchive outarchive;

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

static const double delta(1.0e-8);

// -------------------------------------------------------------
// make_a_data_collection
// -------------------------------------------------------------
gridpack::component::DataCollection *
make_a_data_collection(char *key, const int& v, const char *cval)
{
  gridpack::component::DataCollection *result = 
    new gridpack::component::DataCollection();

  bool bval(true);
  int ival(v);
  long lval(v);
  float fval(v);
  double dval(v);
  gridpack::ComplexType cplxval(v, v);
  
  std::string sval(cval);

  result->addValue(key, bval);
  result->addValue(key, ival);
  result->addValue(key, lval);
  result->addValue(key, fval);
  result->addValue(key, dval);
  result->addValue(key, cval);
  result->addValue(key, cplxval);

  return result;
  
}

// -------------------------------------------------------------
// check_data_collection
// -------------------------------------------------------------
void
check_data_collection(char *key,
                      gridpack::component::DataCollection& dcin, 
                      gridpack::component::DataCollection& dcout)
{
  bool bvalin, bvalout;
  int ivalin, ivalout;
  long lvalin, lvalout;
  float fvalin, fvalout;
  double dvalin, dvalout;
  gridpack::ComplexType cplxin, cplxout;
  std::string svalin, svalout;

  dcin.getValue(key, &bvalin);
  dcout.getValue(key, &bvalout);
  BOOST_CHECK_EQUAL(bvalin, bvalout);

  dcin.getValue(key, &ivalin);
  dcout.getValue(key, &ivalout);
  BOOST_CHECK_EQUAL(ivalin, ivalout);

  dcin.getValue(key, &lvalin);
  dcout.getValue(key, &lvalout);
  BOOST_CHECK_EQUAL(lvalin, lvalout);

  dcin.getValue(key, &fvalin);
  dcout.getValue(key, &fvalout);
  BOOST_CHECK_CLOSE(fvalin, fvalout, delta);

  dcin.getValue(key, &dvalin);
  dcout.getValue(key, &dvalout);
  BOOST_CHECK_CLOSE(dvalin, dvalout, delta);

  dcin.getValue(key, &cplxin);
  dcout.getValue(key, &cplxout);
  BOOST_CHECK_CLOSE(abs(cplxin), abs(cplxout), delta);
  BOOST_CHECK_CLOSE(real(cplxin), real(cplxout), delta);

  dcin.getValue(key, &svalin);
  dcout.getValue(key, &svalout);
  BOOST_CHECK_EQUAL(svalin, svalout);

}


BOOST_AUTO_TEST_SUITE(ComponentSerialization)

BOOST_AUTO_TEST_CASE( DataCollection_bin )
{
  char key[] = "key name";
  boost::scoped_ptr<gridpack::component::DataCollection> 
    dcin(make_a_data_collection(key, 14, "A string value")),
    dcout;

  std::stringstream obuf;
  { 
    outarchive oa(obuf);
    oa << dcin;
  }

  { 
    inarchive ia(obuf);
    ia >> dcout;
  }

  check_data_collection(key, *dcin, *dcout);
}

BOOST_AUTO_TEST_CASE( DataCollection_mpi )
{
  boost::mpi::communicator world;
  char key[] = "key name";

  // each process makes an instance

  typedef boost::shared_ptr<gridpack::component::DataCollection> DCPtr;

  std::string s(boost::lexical_cast<std::string>(world.rank()));
  DCPtr dcin(make_a_data_collection(key, world.rank(), s.c_str()));

  // all instances are gathered to the root process

  std::vector<DCPtr> dcvector(world.size());
  if (world.rank() == 0) {
    gather(world, dcin, dcvector, 0);
  } else {
    gather(world, dcin, 0);
  }

  // all instances are scattered back to the original process

  DCPtr dcout;
  if (world.rank() == 0) {
    scatter(world, dcvector, dcout, 0);
  } else {
    scatter(world, dcout, 0);
  }

  // check to make sure it did not change

  check_data_collection(key, *dcin, *dcout);
}

BOOST_AUTO_TEST_SUITE_END()


// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}
