/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   adjacency_list.hpp
 * @author William A. Perkins
 * @date   2013-12-11 09:25:04 d3g096
 * 
 * @brief  Declaration of the AdjacencyList class
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 17, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _adjacency_list_hpp_
#define _adjacency_list_hpp_

#include <utility>
#include <vector>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>

namespace gridpack {
namespace network {


// -------------------------------------------------------------
//  class AdjacencyList
// -------------------------------------------------------------
/**
 * This class provides a way to assemble a graph adjacency list when
 * graph node and edge information is arbitrarily distributed.  For
 * example, one process may "own" a certain node, but the way that
 * node is connected to others may be on another process.
 * 
 * It just works with integral indexes.  It is assumed that the global
 * indexes start with 0, are unique, and that the largest index is N-1
 * if there are N nodes/edges.  There's probably some special term for
 * that.
 * 
 */
class AdjacencyList : 
    public parallel::Distributed,
    private utility::Uncopyable
{
public:

  /// In case we need to change the type
  typedef unsigned int Index;

  /// A thing to hold Indexes
  typedef std::vector<Index> IndexVector;

  /// The index that means unconnected
  static const Index bogus;

  /// Default constructor.
  AdjacencyList(const parallel::Communicator& comm);

  /// Construct with known local sizes (guesses to size containers, maybe)
  AdjacencyList(const parallel::Communicator& comm,
                const int& local_nodes, const int& local_edges);

  /// Destructor
  ~AdjacencyList(void);

  /// Add the global index of a local node
  void add_node(const Index& node_index)
  {
    p_nodes.push_back(node_index);
  }
  
  /// Add the global index of a local edge and what it connects
  void add_edge(const Index& edge_index, 
                const Index& node_index_1,
                const Index& node_index_2)
  {
    p_Edge tmp;
    tmp.index = edge_index;
    tmp.conn = std::make_pair<Index, Index>(node_index_1, node_index_2);
    p_edges.push_back(tmp);
  }

  /// Get the number of local nodes
  size_t nodes(void) const
  {
    return p_nodes.size();
  }

  /// Get the global node index given a local index
  Index node_index(const int& local_index) const;

  /// Get the number of local edges
  size_t edges(void) const
  {
    return p_edges.size();
  }

  /// Get the global edge index given a local index
  Index edge_index(const int& local_index) const;

  /// Get an edges connected global node indexes 
  void edge(const int& local_index, Index& node1, Index& node2) const;

  /// Indicate that the graph is complete
  void ready(void);

  /// Get the neighbors of the specified (local) node
  void node_neighbors(const int& local_index,
                      IndexVector& global_neighbor_indexes) const;

  /// Get the number of neighbors of the specified (local) node
  size_t node_neighbors(const int& local_index) const;

protected:

  typedef std::pair<Index, Index> p_NodeConnect;
  typedef std::pair<bool, bool> p_Connected;

  struct p_Edge {
    Index index;
    p_NodeConnect conn;
    p_Connected found;
    p_Edge() : index(0), conn(), found(false, false) {}
  };
  typedef std::vector<p_Edge> p_EdgeVector;

  typedef std::vector<IndexVector> p_Adjacency;

  /// The list of local node indexes
  IndexVector p_nodes;
  
  /// The list of local edges
  p_EdgeVector p_edges;

  /// The resulting adjacency for local nodes
  p_Adjacency p_adjacency;
  

};

} // namespace network
} // namespace gridpack


#endif
