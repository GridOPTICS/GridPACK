// -------------------------------------------------------------
/**
 * @file   adjacency_list.cpp
 * @author William A. Perkins
 * @date   2013-06-18 12:04:11 d3g096
 * 
 * @brief  Implementation of AdjacencyList
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 17, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <boost/assert.hpp>
#include <boost/mpi/collectives.hpp>
#include "adjacency_list.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class AdjacencyList
// -------------------------------------------------------------

// -------------------------------------------------------------
// AdjacencyList static members
// -------------------------------------------------------------
const AdjacencyList::Index AdjacencyList::bogus(-1);

// -------------------------------------------------------------
// AdjacencyList:: constructors / destructor
// -------------------------------------------------------------
AdjacencyList::AdjacencyList(const parallel::Communicator& comm)
  : parallel::Distributed(comm),
    utility::Uncopyable(),
    p_nodes(), p_edges(), p_adjacency()
{
  // empty
}

AdjacencyList::AdjacencyList(const parallel::Communicator& comm,
                             const int& local_nodes, const int& local_edges)
  : parallel::Distributed(comm),
    utility::Uncopyable(),
    p_nodes(), p_edges(), p_adjacency()
{
  p_nodes.reserve(local_nodes);
  p_edges.reserve(local_edges);
  p_adjacency.reserve(local_nodes);
}

AdjacencyList::~AdjacencyList(void)
{
  // empty
}

// -------------------------------------------------------------
// AdjacencyList::node_index
// -------------------------------------------------------------
AdjacencyList::Index 
AdjacencyList::node_index(const int& local_index) const
{
  BOOST_ASSERT(local_index < this->nodes());
  return p_nodes[local_index];
}

// -------------------------------------------------------------
// AdjacencyList::edge_index
// -------------------------------------------------------------
AdjacencyList::Index 
AdjacencyList::edge_index(const int& local_index) const
{
  BOOST_ASSERT(local_index < this->edges());
  return p_edges[local_index].index;
}



// -------------------------------------------------------------
// AdjacencyList::ready
// -------------------------------------------------------------
void
AdjacencyList::ready(void)
{
  int me(this->processor_rank());
  int nproc(this->processor_size());

  p_adjacency.clear();
  p_adjacency.resize(p_nodes.size());

  IndexVector current_indexes;
  IndexVector connected_indexes;

  for (int p = 0; p < nproc; ++p) {

    // broadcast the node indexes owned by process p to all processes,
    // all processes work on these at once

    current_indexes.clear();
    if (me == p) {
      std::copy(p_nodes.begin(), p_nodes.end(), 
                std::back_inserter(current_indexes));
      std::cout << me << ": node indexes: ";
      std::copy(current_indexes.begin(), current_indexes.end(),
                std::ostream_iterator<Index>(std::cout, ","));
      std::cout << std::endl;
    }
    broadcast(this->communicator(), current_indexes, p);

    // loop over the process p's node index set
    
    int local_index(0);
    for (IndexVector::iterator n = current_indexes.begin(); 
         n != current_indexes.end(); ++n, ++local_index) {
      std::cout << me << ": current node index: " << *n << std::endl;
      
      // determine the local edges that refer to the current node index
 
      connected_indexes.clear();
      for (p_EdgeVector::iterator e = p_edges.begin();
           e != p_edges.end(); ++e) {
        if (*n == e->conn.first && e->conn.second != bogus) {
          connected_indexes.push_back(e->conn.second);
          std::cout << me << ": found connection: edge " << e->index
                    << " (" << e->conn.first << ", " << e->conn.second << ")"
                    << std::endl;
        }
        if (*n == e->conn.second && e->conn.first != bogus) {
          connected_indexes.push_back(e->conn.first);
          std::cout << me << ": found connection: edge " << e->index
                    << " (" << e->conn.first << ", " << e->conn.second << ")"
                    << std::endl;
        }
      }

      // gather all connections for the current node index to the
      // node's owner process, we have to gather the vectors because
      // processes will have different numbers of connections

      if (me == p) {
        size_t allsize;
        reduce(this->communicator(), 
               connected_indexes.size(), allsize, std::plus<size_t>(), p);

        std::vector<IndexVector> all_connected_indexes;
        gather(this->communicator(), 
               connected_indexes, all_connected_indexes, p);
        p_adjacency[local_index].clear();
        for (std::vector<IndexVector>::iterator k = all_connected_indexes.begin();
             k != all_connected_indexes.end(); ++k) {
          std::copy(k->begin(), k->end(), 
                    std::back_inserter(p_adjacency[local_index]));
        }
      } else {
        reduce(this->communicator(), 
               connected_indexes.size(), std::plus<size_t>(), p);
        gather(this->communicator(), connected_indexes, p);
      }
      this->communicator().barrier();
    }
    this->communicator().barrier();
  }
}

// -------------------------------------------------------------
// AdjacencyList::node_neighbors
// -------------------------------------------------------------
void
AdjacencyList::node_neighbors(const int& local_index,
                              IndexVector& global_neighbor_indexes) const
{
  BOOST_ASSERT(local_index < p_adjacency.size());
  global_neighbor_indexes.clear();
  std::copy(p_adjacency[local_index].begin(), p_adjacency[local_index].end(),
            std::back_inserter(global_neighbor_indexes));

}


} // namespace network
} // namespace gridpack
