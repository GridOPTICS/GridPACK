// -------------------------------------------------------------
/**
 * @file   graph_partitioner_implementation.cpp
 * @author William A. Perkins
 * @date   2013-07-24 14:17:55 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <algorithm>
#include <ga++.h>
#include <boost/mpi/collectives.hpp>
#include <boost/scoped_ptr.hpp>
#include "graph_partitioner_implementation.hpp"


namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class GraphPartitionerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// GraphPartitionerImplementation:: constructors / destructor
// -------------------------------------------------------------
GraphPartitionerImplementation::GraphPartitionerImplementation(const parallel::Communicator& comm)
  : parallel::Distributed(comm), utility::Uncopyable(),
    p_adjacency_list(comm), 
    p_node_destinations(),
    p_edge_destinations()
{
  // empty
}

GraphPartitionerImplementation::GraphPartitionerImplementation(const parallel::Communicator& comm,
                                                               const int& local_nodes, 
                                                               const int& local_edges)
  : parallel::Distributed(comm), utility::Uncopyable(),
    p_adjacency_list(comm, local_nodes, local_edges), 
    p_node_destinations(local_nodes),
    p_edge_destinations(local_edges)
{
  // empty
}

GraphPartitionerImplementation::~GraphPartitionerImplementation(void)
{
}

// -------------------------------------------------------------
// GraphPartitionerImplementation::node_destinations
// -------------------------------------------------------------
void
GraphPartitionerImplementation::node_destinations(GraphPartitionerImplementation::IndexVector& dest) const
{
  dest.clear();
  std::copy(p_node_destinations.begin(), p_node_destinations.end(),
            std::back_inserter(dest));
}

// -------------------------------------------------------------
// GraphPartitionerImplementation::edge_destinations
// -------------------------------------------------------------
void
GraphPartitionerImplementation::edge_destinations(GraphPartitionerImplementation::IndexVector& dest) const
{
  dest.clear();
  std::copy(p_edge_destinations.begin(), p_edge_destinations.end(),
            std::back_inserter(dest));
}

// -------------------------------------------------------------
// GraphPartitionerImplementation::ghost_node_destinations
// -------------------------------------------------------------
void
GraphPartitionerImplementation::ghost_node_destinations(GraphPartitionerImplementation::MultiIndexVector& dest) const
{
  dest.clear();
  std::copy(p_ghost_node_destinations.begin(), 
            p_ghost_node_destinations.end(),
            std::back_inserter(dest));
}

// -------------------------------------------------------------
// GraphPartitionerImplementation::ghost_edge_destinations
// -------------------------------------------------------------
void
GraphPartitionerImplementation::ghost_edge_destinations(GraphPartitionerImplementation::IndexVector& dest) const
{
  dest.clear();
  std::copy(p_ghost_edge_destinations.begin(), 
            p_ghost_edge_destinations.end(),
            std::back_inserter(dest));
}

// -------------------------------------------------------------
// GraphPartitionerImplementation::partition
// -------------------------------------------------------------
void
GraphPartitionerImplementation::partition(void)
{
  p_adjacency_list.ready();
  this->p_partition();          // fills p_node_destinations

  int maxdim(1);
  int dims[maxdim], lo[maxdim], hi[maxdim], ld[maxdim];
  ld[0] = 1;

  int locnodes(p_adjacency_list.nodes());
  int locedges(p_adjacency_list.edges());
  int allnodes;
  int alledges;

  all_reduce(this->communicator(), locnodes, allnodes, std::plus<int>());
  all_reduce(this->communicator(), locedges, alledges, std::plus<int>());

  // make two GAs, one that holds the node source and another that
  // node destination; each is indexed by global node index

  std::vector<int> nodeidx(locnodes);
  std::vector<int *> stupid(locnodes);
  for (Index n = 0; n < locnodes; ++n) {
    nodeidx[n] = p_adjacency_list.node_index(n);
    stupid[n] = &nodeidx[n];
  }

  dims[0] = allnodes;
  boost::scoped_ptr<GA::GlobalArray> 
    node_dest(new GA::GlobalArray(MT_C_INT, 1, dims, "Node Destinations Process", NULL)),
    node_src(new GA::GlobalArray(MT_C_INT, 1, dims, "Node Source Process", NULL));
  node_dest->scatter(&p_node_destinations[0], &stupid[0], locnodes);

  { 
    std::vector<int> nsrc(locnodes, this->processor_rank());
    node_src->scatter(&nsrc[0], &stupid[0], locnodes);
  }
  
  node_src->print();
  node_dest->print();

  // edges are assigned to the same partition as the lowest numbered
  // node to which it connects, which are extracted from the node
  // destination GA.

  nodeidx.resize(locedges);
  stupid.resize(locedges);
  std::vector<int> e1dest(locedges);

  for (Index e = 0; e < locedges; ++e) {
    Index n1, n2;
    p_adjacency_list.edge(e, n1, n2);
    nodeidx[e] = std::min(n1, n2);
    stupid[e] = &nodeidx[e];
  }

  node_dest->gather(&e1dest[0], &stupid[0], locedges);

  p_edge_destinations.clear();
  p_edge_destinations.reserve(locedges);
  std::copy(e1dest.begin(), e1dest.end(), 
            std::back_inserter(p_edge_destinations));

  // determine (possible) destinations for ghost edges (highest numbered node) 

  std::vector<int> e2dest(locedges);
  for (Index e = 0; e < locedges; ++e) {
    Index n1, n2;
    p_adjacency_list.edge(e, n1, n2);
    nodeidx[e] = std::max(n1, n2);
    stupid[e] = &nodeidx[e];
  }

  node_dest->gather(&e2dest[0], &stupid[0], locedges);

  p_ghost_edge_destinations.reserve(locedges);
  std::copy(e2dest.begin(), e2dest.end(), 
            std::back_inserter(p_ghost_edge_destinations));

}
} // namespace network
} // namespace gridpack
