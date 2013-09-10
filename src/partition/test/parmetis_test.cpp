/**
 * @file   parmetis_test.cpp
 * @author William A. Perkins
 * @date   2013-09-10 14:28:45 d3g096
 * 
 * @brief  Unit tests of ParMETIS-specific code
 * 
 * @test
 */

#include <ga++.h>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "simple_adjacency.hpp"
#include "parmetis/parmetis_graph_wrapper.hpp"
#include "gridpack/parallel/printit.hpp"


BOOST_AUTO_TEST_SUITE( ParMETISTest )

BOOST_AUTO_TEST_CASE( graph_wrapper )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(5*world.size());

  using gridpack::network::AdjacencyList;
  using gridpack::network::ParMETISGraphWrapper;

  std::auto_ptr<AdjacencyList> 
    adlist(simple_adjacency_list(world, global_nodes));

  std::auto_ptr<ParMETISGraphWrapper>
    wrapper(new ParMETISGraphWrapper(*adlist));

  std::vector<idx_t> vtxdist;
  std::vector<idx_t> xadj;
  std::vector<idx_t> adjncy;

  wrapper->get_csr_local(vtxdist, xadj, adjncy);

  printit(world, vtxdist, "vtxdist:");
  printit(world, xadj, "xadj:");
  printit(world, adjncy, "adjncy:");

  // Not really sure what to check here

  BOOST_CHECK(true);

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
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  GA::Terminate();
}

