/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   partition_test.cpp
 * @author William A. Perkins
 * @date   2013-09-10 14:23:17 d3g096
 * 
 * @brief  Unit test suite for various partition classes.
 * 
 * @test
 */
// -------------------------------------------------------------

#include <ctime>
#include <iostream>
#include <iterator>
#include <ga++.h>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "simple_adjacency.hpp"
#include "graph_partitioner.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "simple_adjacency.hpp"

BOOST_AUTO_TEST_SUITE( PartitionTest )


// -------------------------------------------------------------
// adjacency unit test
// -------------------------------------------------------------
/// 
/** 
 * @test
 * 
 * 
 */
BOOST_AUTO_TEST_CASE( adjacency )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(4*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::AdjacencyList;

  std::auto_ptr<AdjacencyList> 
    adlist(simple_adjacency_list(world, global_nodes));

  for (int p = 0; p < world.size(); ++p) {
    if (world.rank() == p) {
      AdjacencyList::IndexVector nbr;
      std::cout << p << ": Results" << std::endl;
      for (size_t i = 0; i < adlist->nodes(); ++i) {
        adlist->node_neighbors(i, nbr);
        std::cout << p << ": node " << adlist->node_index(i) << ": ";
        std::copy(nbr.begin(), nbr.end(),
                  std::ostream_iterator<AdjacencyList::Index>(std::cout, ","));
        std::cout << std::endl;
      }
    }
    world.barrier();
  }

  for (size_t i = 0; i < adlist->nodes(); ++i) {
    AdjacencyList::Index node(adlist->node_index(i));
    AdjacencyList::IndexVector nbr;
    adlist->node_neighbors(i, nbr);

    // all nodes should connect to at least one other, but no more than two
    BOOST_CHECK(nbr.size() >= 1 && nbr.size() <= 2);

    // node n should connect to node n-1, except for 0
    AdjacencyList::IndexVector::iterator p;
    if (node > 0) {
      p = std::find(nbr.begin(), nbr.end(), node - 1);
      BOOST_CHECK(p != nbr.end());
    }

    // node n should connect to node n+1 except for the last
    if (node < global_nodes - 1) {
      p = std::find(nbr.begin(), nbr.end(), node + 1);
      BOOST_CHECK(p != nbr.end());
    }
  }
    
}

/// 
/**
 * @test
 * 
 */
BOOST_AUTO_TEST_CASE( random_partition )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(5*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::GraphPartitioner;

  std::auto_ptr<GraphPartitioner> 
    partitioner(simple_graph_partitioner(world, global_nodes));

  partitioner->partition();

  world.barrier();

  GraphPartitioner::IndexVector node_dest, edge_dest;
  partitioner->node_destinations(node_dest);
  partitioner->edge_destinations(edge_dest);
  for (int p = 0; p < world.size(); ++p) {
    if (world.rank() == p) {
      int nnodes(partitioner->nodes());
      std::cout << p << ": nodes: ";
      for (int i = 0; i < nnodes; ++i) {
        std::cout << partitioner->node_index(i) << ",";
      }
      std::cout << std::endl;
      std::cout << p << ": parts: ";
      std::copy(node_dest.begin(), node_dest.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
      int nedges(partitioner->edges());
      std::cout << p << ": edges: ";
      for (int i = 0; i < nedges; ++i) {
        std::cout << partitioner->edge_index(i) << ",";
      }
      std::cout << std::endl;
      std::cout << p << ": parts: ";
      std::copy(edge_dest.begin(), edge_dest.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
    }
    world.barrier();
  }
}

/// Partition an extremely unbalanced graph
/**
 * @test
 * 
 * In this test, a simple, linear graph is created entirely on process
 * zero.  The graph is then partitioned evenly amongst all
 * participating processes.
 */
BOOST_AUTO_TEST_CASE( unbalanced_partition )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(5*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::GraphPartitioner;
  using gridpack::network::AdjacencyList;

  GraphPartitioner::IndexVector my_nodes;
  GraphPartitioner::IndexVector my_edges;

  if (world.rank() == 0) {
    my_nodes.resize(global_nodes);
    my_edges.resize(global_edges);
    for (GraphPartitioner::Index i = 0; i < my_nodes.size(); ++i) {
      my_nodes[i] = i;
    }
    for (GraphPartitioner::Index i = 0; i < my_edges.size(); ++i) {
      my_edges[i] = i;
    }
  }

  GraphPartitioner partitioner(world);

  for (GraphPartitioner::IndexVector::iterator n = my_nodes.begin();
       n != my_nodes.end(); ++n) {
    partitioner.add_node(*n);
  }
  for (GraphPartitioner::IndexVector::iterator e = my_edges.begin();
       e != my_edges.end(); ++e) {
    GraphPartitioner::Index from(*e), to(from+1);
    if (from < 0) from = AdjacencyList::bogus;
    if (to > global_nodes - 1) to = AdjacencyList::bogus;
    partitioner.add_edge(*e, from, to);
  }

  partitioner.partition();

  world.barrier();

  GraphPartitioner::IndexVector node_dest, edge_dest;
  partitioner.node_destinations(node_dest);
  partitioner.edge_destinations(edge_dest);
  for (int p = 0; p < world.size(); ++p) {
    if (world.rank() == p) {
      std::cout << p << ": nodes: ";
      std::copy(my_nodes.begin(), my_nodes.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
      std::cout << p << ": parts: ";
      std::copy(node_dest.begin(), node_dest.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
      std::cout << p << ": edges: ";
      std::copy(my_edges.begin(), my_edges.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
      std::cout << p << ": parts: ";
      std::copy(edge_dest.begin(), edge_dest.end(),
                std::ostream_iterator<GraphPartitioner::Index>(std::cout, ","));
      std::cout << std::endl;
    }
    world.barrier();
  }

  /// @todo FIXME: need to test some stuff here

}

BOOST_AUTO_TEST_SUITE_END()


// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
/**
 * @test
 * 
 */
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
/**
 * @test
 * 
 */
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  GA::Terminate();
}




