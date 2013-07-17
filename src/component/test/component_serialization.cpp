/**
 * @file   component_serialization.cpp
 * @author William A. Perkins
 * @date   2013-07-11 14:43:59 d3g096
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
#include <boost/lexical_cast.hpp>

#include <boost/serialization/scoped_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

typedef boost::archive::binary_iarchive inarchive;
typedef boost::archive::binary_oarchive outarchive;

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "data_collection.hpp"
#include "base_component.hpp"

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

// -------------------------------------------------------------
// gather_scatter
// -------------------------------------------------------------
/// Test MPI serialization using using gather and scatter

/** 
 * Each process creates a thing.  All things are gathered to the zero
 * process then scattered back to each process.  
 * 
 * Thing must serializable and copyable.
 * 
 * @param thing_in original thing
 * @param thing_out thing after gather/scatter
 */
template <typename Thing> void
gather_scatter(const boost::mpi::communicator& comm,
               const Thing& thing_in, Thing& thing_out)
{
  std::vector<Thing> thing_vector(comm.size());
  if (comm.rank() == 0) {
    gather(comm, thing_in, thing_vector, 0);
  } else {
    gather(comm, thing_in, 0);
  }

  if (comm.rank() == 0) {
    scatter(comm, thing_vector, thing_out, 0);
  } else {
    scatter(comm, thing_out, 0);
  }
}

// -------------------------------------------------------------
//  class BogusBus
// -------------------------------------------------------------
class BogusBus 
  : public gridpack::component::BaseBusComponent
{
public:

  /// Default constructor.
  explicit BogusBus(const int& id)
    : gridpack::component::BaseBusComponent(), p_name("BogusBus#")
  {
    this->setMatVecIndex(id);
    p_name += boost::lexical_cast<std::string>(id);
  }

  /// Destructor
  ~BogusBus(void) {}

  /// Get the name
  const std::string& name(void) const
  {
    return p_name;
  }

protected:

  /// A name
  std::string p_name;

private:

  /// Constructor for serialization only
  BogusBus(void) 
    : gridpack::component::BaseBusComponent(), p_name()
  {
  }

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseBusComponent>(*this)
       & p_name;
  }

};

BOOST_CLASS_EXPORT(BogusBus);

// -------------------------------------------------------------
//  class BogusBranch
// -------------------------------------------------------------
class BogusBranch 
  : public gridpack::component::BaseBranchComponent
{
public:

  /// Default constructor.
  explicit BogusBranch(const int& id)
    : gridpack::component::BaseBranchComponent(), p_name("BogusBranch#")
  {
    this->setMatVecIndices(id, 0);
    p_name += boost::lexical_cast<std::string>(id);
  }

  /// Destructor
  ~BogusBranch(void) {}

  /// Get the name
  const std::string& name(void) const
  {
    return p_name;
  }

protected:

  /// A name
  std::string p_name;

private:

  /// Constructor for serialization only
  BogusBranch(void) 
    : gridpack::component::BaseBranchComponent(), p_name()
  {
  }

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseBranchComponent>(*this)
       & p_name;
  }

};

BOOST_CLASS_EXPORT(BogusBranch);

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
  DCPtr dcout;

  gather_scatter(world, dcin, dcout);

  // check to make sure it did not change

  check_data_collection(key, *dcin, *dcout);
}

BOOST_AUTO_TEST_CASE ( Component_bin )
{
  static int the_id(1);
  boost::scoped_ptr<gridpack::component::BaseComponent> 
    busin(new BogusBus(the_id)), 
    busout,
    branchin(new BogusBranch(the_id)),
    branchout;


  std::stringstream obuf;
  { 
    outarchive oa(obuf);
    oa << busin << branchin;
  }

  { 
    inarchive ia(obuf);
    ia >> busout >> branchout;
  }
  int inid;
  int outid;
  int junk;

  BogusBus 
    *bogusin = dynamic_cast<BogusBus *>(busin.get()),
    *bogusout = dynamic_cast<BogusBus *>(busout.get());
  
  BOOST_REQUIRE(bogusin != NULL);
  BOOST_REQUIRE(bogusout != NULL);
  
  bogusin->getMatVecIndex(&inid);
  bogusout->getMatVecIndex(&outid);
  
  BOOST_CHECK_EQUAL(inid, the_id);
  BOOST_CHECK_EQUAL(inid, outid);
  
  BOOST_REQUIRE(!(bogusin->name()).empty());
  BOOST_CHECK_EQUAL(bogusin->name(), bogusout->name());

  BogusBranch
    *brogusin = dynamic_cast<BogusBranch *>(branchin.get()),
    *brogusout = dynamic_cast<BogusBranch *>(branchout.get());

  BOOST_REQUIRE(brogusin != NULL);
  BOOST_REQUIRE(brogusout != NULL);
  
  brogusin->getMatVecIndices(&inid, &junk);
  brogusout->getMatVecIndices(&outid, &junk);
  
  BOOST_CHECK_EQUAL(inid, the_id);
  BOOST_CHECK_EQUAL(inid, outid);
  
  BOOST_REQUIRE(!(brogusin->name()).empty());
  BOOST_CHECK_EQUAL(brogusin->name(), brogusout->name());
}

BOOST_AUTO_TEST_CASE ( Component_mpi )
{

  boost::mpi::communicator world;

  typedef boost::shared_ptr<gridpack::component::BaseComponent> CompPtr;
  CompPtr compin, compout;
  int inid, outid, junk;

  // check BogusBus 

  compin.reset(new BogusBus(world.rank()));
  compout.reset();

  gather_scatter(world, compin, compout);

  BogusBus 
    *bogusin = dynamic_cast<BogusBus *>(compin.get()),
    *bogusout = dynamic_cast<BogusBus *>(compout.get());
  
  BOOST_REQUIRE(bogusin != NULL);
  BOOST_REQUIRE(bogusout != NULL);
  
  bogusin->getMatVecIndex(&inid);
  bogusout->getMatVecIndex(&outid);
  
  BOOST_CHECK_EQUAL(inid, world.rank());
  BOOST_CHECK_EQUAL(inid, outid);
  
  BOOST_REQUIRE(!(bogusin->name()).empty());
  BOOST_CHECK_EQUAL(bogusin->name(), bogusout->name());

  // check BogusBranch

  compin.reset(new BogusBranch(world.rank()));
  compout.reset();

  gather_scatter(world, compin, compout);

  BogusBranch
    *brogusin = dynamic_cast<BogusBranch *>(compin.get()),
    *brogusout = dynamic_cast<BogusBranch *>(compout.get());

  BOOST_REQUIRE(brogusin != NULL);
  BOOST_REQUIRE(brogusout != NULL);
  
  brogusin->getMatVecIndices(&inid, &junk);
  brogusout->getMatVecIndices(&outid, &junk);
  
  BOOST_CHECK_EQUAL(inid, world.rank());
  BOOST_CHECK_EQUAL(inid, outid);
  
  BOOST_REQUIRE(!(brogusin->name()).empty());
  BOOST_CHECK_EQUAL(brogusin->name(), brogusout->name());

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
