/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   graph_partioner_implementation.hpp
 * @author William A. Perkins
 * @date   2013-07-24 14:17:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _graph_partioner_implementation_hpp_
#define _graph_partioner_implementation_hpp_

#include <gridpack/partition/adjacency_list.hpp>

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class GraphPartitionerImplementation
// -------------------------------------------------------------
class GraphPartitionerImplementation 
  : public parallel::Distributed,
    private utility::Uncopyable
{
public:

  /// The index type
  typedef AdjacencyList::Index Index;

  /// A vector of Indexes
  typedef AdjacencyList::IndexVector IndexVector;

  /// A vector of IndexVectors
  typedef std::vector<IndexVector> MultiIndexVector;

  /// Default constructor.
  GraphPartitionerImplementation(const parallel::Communicator& comm);

  /// Construct w/ known local sizes  (guesses to size containers, maybe)
  GraphPartitionerImplementation(const parallel::Communicator& comm,
                                 const int& local_nodes, const int& local_edges);

  /// Destructor
  virtual ~GraphPartitionerImplementation(void);

  /// Add the global index and original index of a local node
  void add_node(const Index& global_index, const Index& original_index)
  {
    p_adjacency_list.add_node(global_index, original_index);
  }
  
  /// Add the global index of a local edge and what it connects using the
  /// original indices of buses at either end
  void add_edge(const Index& edge_index, 
                const Index& node_index_1,
                const Index& node_index_2)
  {
    p_adjacency_list.add_edge(edge_index, node_index_1, node_index_2);
  }

  /// Get the global indices of the buses at either end of a branch
  void get_global_edge_ids(int idx, Index *node_index_1, Index *node_index_2) const
  {
    p_adjacency_list.get_global_edge_ids(idx, node_index_1, node_index_2);
  }

  /// Get the number of local nodes
  size_t nodes(void) const
  {
    return p_adjacency_list.nodes();
  }

  /// Get the global node index given a local index
  Index node_index(const int& local_index) const
  {
    return p_adjacency_list.node_index(local_index);
  }

  /// Get the number of local edges
  size_t edges(void) const
  {
    return p_adjacency_list.edges();
  }

  /// Get the global edge index given a local index
  Index edge_index(const int& local_index) const
  {
    return p_adjacency_list.edge_index(local_index);
  }

  /// Partition the graph
  void partition(void);

  /// Get the node destinations
  void node_destinations(IndexVector& dest) const;

  /// Get the edge destinations
  void edge_destinations(IndexVector& dest) const;

  /// Get the destinations of ghosted nodes
  void ghost_node_destinations(MultiIndexVector& dest) const;

  /// Get the destinations of ghosted edges
  void ghost_edge_destinations(IndexVector& dest) const;

protected:

  /// Adjacency list builder
  AdjacencyList p_adjacency_list;

  /// A list of processors where local nodes should go
  IndexVector p_node_destinations;

  /// A list of processors where local edges should go
  IndexVector p_edge_destinations;

  /// A list of processors where local nodes should go
  MultiIndexVector p_ghost_node_destinations;

  /// A list of processors where local edges should go
  IndexVector p_ghost_edge_destinations;

  /// Partition the graph (specialized)
  virtual void p_partition(void) = 0;

};


} // namespace network
} // namespace gridpack


#endif
