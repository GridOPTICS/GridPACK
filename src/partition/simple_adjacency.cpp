/**
 * @file   simple_adjacency.cpp
 * @author William A. Perkins
 * @date   2013-07-25 09:35:05 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include <ctime>
#include <iostream>
#include <iterator>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include "gridpack/parallel/parallel.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen(12*time(NULL));

#include "graph_partitioner.hpp"

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
    for (int i = 0; i < global_size; ++i) {
      boost::random::uniform_int_distribution<> dist(0, nproc-1);
      int p(dist(gen));
      toscatter[p].push_back(i);
    }
  }
  scatter(comm, toscatter, local_values, 0);
}

// -------------------------------------------------------------
// simple_adjacency
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
/// Make a simple AdjacencyList instance
/** 
 * This makes a simple AdjacencyList instance that represents a graph
 * like this
 * 
 *  0   1   2   n-2 n-1
 *  +---+---+-//-+---+
 *    0   1       e-1
 * 
 * where n is the number of nodes, e is the number of edges. It
 * follows that e = n - 1.  Each edge's e connects to nodes e-1 and
 * e. In the end, each node n's adjacency 2: n-1 and n+1.  Except on
 * the ends.
 * 
 * 
 * @param global_nodes total number of nodes
 * @param global_edges total number of edges
 * 
 * @return adjacency list instance
 */
gridpack::network::AdjacencyList *
simple_adjacency_list(const gridpack::parallel::Communicator& comm,
                      const int& global_nodes)
{
  using gridpack::network::AdjacencyList;

  AdjacencyList::IndexVector my_nodes;
  AdjacencyList::IndexVector my_edges;
  random_scattered_vector(comm, global_nodes, my_nodes);
  random_scattered_vector(comm, global_nodes-1, my_edges);


  AdjacencyList *adlist(new AdjacencyList(comm));
  for (AdjacencyList::IndexVector::iterator n = my_nodes.begin();
       n != my_nodes.end(); ++n) {
    adlist->add_node(*n);
  }
  for (AdjacencyList::IndexVector::iterator e = my_edges.begin();
       e != my_edges.end(); ++e) {
    AdjacencyList::Index from(*e), to(from+1);
    if (from < 0) from = AdjacencyList::bogus;
    if (to > global_nodes - 1) to = AdjacencyList::bogus;
    adlist->add_edge(*e, from, to);
  }
  adlist->ready();
  return adlist;
}

// -------------------------------------------------------------
// simple_graph_partitioner
// -------------------------------------------------------------
gridpack::network::GraphPartitioner *
simple_graph_partitioner(const gridpack::parallel::Communicator& comm,
                         const int& global_nodes)
{
  using gridpack::network::AdjacencyList;
  using gridpack::network::GraphPartitioner;

  GraphPartitioner::IndexVector my_nodes;
  GraphPartitioner::IndexVector my_edges;
  random_scattered_vector(comm, global_nodes, my_nodes);
  random_scattered_vector(comm, global_nodes-1, my_edges);

  GraphPartitioner *partitioner = 
    new GraphPartitioner(comm);

  for (GraphPartitioner::IndexVector::iterator n = my_nodes.begin();
       n != my_nodes.end(); ++n) {
    partitioner->add_node(*n);
  }
  for (GraphPartitioner::IndexVector::iterator e = my_edges.begin();
       e != my_edges.end(); ++e) {
    GraphPartitioner::Index from(*e), to(from+1);
    if (from < 0) from = AdjacencyList::bogus;
    if (to > global_nodes - 1) to = AdjacencyList::bogus;
    partitioner->add_edge(*e, from, to);
  }
  return partitioner;
}
