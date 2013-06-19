// -------------------------------------------------------------
/**
 * @file   partition_test.cpp
 * @author William A. Perkins
 * @date   2013-06-19 15:48:27 d3g096
 * 
 * @brief  Unit test suite for various partition classes.
 * 
 * 
 */
// -------------------------------------------------------------

#include <ctime>
#include <iostream>
#include <iterator>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "graph_partitioner.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen(12 /*time(NULL)*/);

// -------------------------------------------------------------
// random_vector
// I is some integer type
// -------------------------------------------------------------
template <typename I>
void 
random_scattered_vector(const boost::mpi::communicator& comm, 
                        const int& global_size,
                        std::vector<I>& local_values)
{
  int me(comm.rank());
  int nproc(comm.size());

  std::vector< std::vector<I> > toscatter(nproc);

  if (me == 0) {
    for (I i = 0; i < global_size; ++i) {
      boost::random::uniform_int_distribution<> dist(0, nproc-1);
      int p(dist(gen));
      toscatter[p].push_back(i);
    }
  }
  scatter(comm, toscatter, local_values, 0);
}

BOOST_AUTO_TEST_SUITE( Partition )


// -------------------------------------------------------------
// adjacency
//
// In this test, a simple circular network is imagined:
// 
//  0   1   2   n-2 n-1
//  +---+---+-//-+---+
//    0   1       e-1
// 
// where n is the number of nodes, e is the number of edges. It
// follows that e = n - 1.  Each edge's e connects to nodes e-1 and
// e. In the end, each node n's adjacency 2: n-1 and n+1.  Except on the ends. 
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( adjacency )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(4*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::AdjacencyList;

  AdjacencyList::IndexVector my_nodes;
  AdjacencyList::IndexVector my_edges;
  random_scattered_vector(world, global_nodes, my_nodes);
  random_scattered_vector(world, global_nodes, my_edges);


  AdjacencyList adlist(world);
  for (AdjacencyList::IndexVector::iterator n = my_nodes.begin();
       n != my_nodes.end(); ++n) {
    adlist.add_node(*n);
  }
  for (AdjacencyList::IndexVector::iterator e = my_edges.begin();
       e != my_edges.end(); ++e) {
    AdjacencyList::Index from(*e), to(from+1);
    if (from < 0) from = AdjacencyList::bogus;
    if (to > global_nodes - 1) to = AdjacencyList::bogus;
    adlist.add_edge(*e, from, to);
  }
  adlist.ready();

  for (int p = 0; p < world.size(); ++p) {
    if (world.rank() == p) {
      std::cout << p << ": Starting with " << std::endl;
      std::cout << p << ": node indexes: ";
      std::copy(my_nodes.begin(), my_nodes.end(),
                std::ostream_iterator<AdjacencyList::Index>(std::cout, ","));
      std::cout << std::endl;
      std::cout << p << ": edge indexes: ";
      std::copy(my_edges.begin(), my_edges.end(),
                std::ostream_iterator<AdjacencyList::Index>(std::cout, ","));
      std::cout << std::endl;
      
      AdjacencyList::IndexVector nbr;
      std::cout << p << ": Results" << std::endl;
      for (size_t i = 0; i < adlist.nodes(); ++i) {
        adlist.node_neighbors(i, nbr);
        std::cout << p << ": node " << adlist.node_index(i) << ": ";
        std::copy(nbr.begin(), nbr.end(),
                  std::ostream_iterator<AdjacencyList::Index>(std::cout, ","));
        std::cout << std::endl;
      }
    }
    world.barrier();
  }

  for (size_t i = 0; i < adlist.nodes(); ++i) {
    AdjacencyList::Index node(adlist.node_index(i));
    AdjacencyList::IndexVector nbr;
    adlist.node_neighbors(i, nbr);

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

BOOST_AUTO_TEST_CASE( random_partition )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(5*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::GraphPartitioner;
  using gridpack::network::AdjacencyList;

  GraphPartitioner::IndexVector my_nodes;
  GraphPartitioner::IndexVector my_edges;
  random_scattered_vector(world, global_nodes, my_nodes);
  random_scattered_vector(world, global_nodes, my_edges);

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

  GraphPartitioner::IndexVector node_dest;
  partitioner.node_destinations(node_dest);
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
    }
    world.barrier();
  }
}



#if 0

// This does not work with ParMETIS

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

  GraphPartitioner::IndexVector node_dest;
  partitioner.node_destinations(node_dest);
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
    }
    world.barrier();
  }
}

#endif

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




