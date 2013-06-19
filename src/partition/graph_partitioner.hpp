// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   graph_partitioner.hpp
 * @author William A. Perkins
 * @date   2013-06-18 12:43:19 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _graph_partitioner_hpp_
#define _graph_partitioner_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/partition/graph_partitioner_implementation.hpp>

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class GraphPartitioner
// -------------------------------------------------------------
/// A class that serves as an interface to a graph partitioning library
class GraphPartitioner : private utility::Uncopyable {
public:

  typedef GraphPartitionerImplementation::Index Index;
  typedef GraphPartitionerImplementation::IndexVector IndexVector;

  /// Default constructor.
  GraphPartitioner(const parallel::Communicator& comm);

  /// Construct w/ known local sizes  (guesses to size containers, maybe)
  GraphPartitioner(const parallel::Communicator& comm,
                   const int& local_nodes, const int& local_edges);

  /// Destructor
  ~GraphPartitioner(void);

  /// Add the global index of a local node
  void add_node(const Index& node_index)
  {
    p_impl->add_node(node_index);
  }
  
  /// Add the global index of a local edge and what it connects
  void add_edge(const Index& edge_index, 
                const Index& node_index_1,
                const Index& node_index_2)
  {
    p_impl->add_edge(edge_index, node_index_1, node_index_2);
  }

  /// Get the number of local nodes
  size_t nodes(void) const
  {
    return p_impl->nodes();
  }

  /// Get the global node index given a local index
  Index node_index(const int& local_index) const
  {
    return p_impl->node_index(local_index);
  }

  /// Get the number of local edges
  size_t edges(void) const
  {
    return p_impl->edges();
  }

  /// Get the global edge index given a local index
  Index edge_index(const int& local_index) const
  {
    return p_impl->edge_index(local_index);
  }

  /// Partition the graph
  void partition(void)
  {
    p_impl->partition();
  }

  /// Get the node destinations
  void node_destinations(IndexVector& dest) const
  {
    p_impl->node_destinations(dest);
  }

  /// Get the edge destinations
  void edge_destinations(IndexVector& dest) const
  {
    p_impl->edge_destinations(dest);
  }

protected:

  /// The actual implementation
  boost::scoped_ptr<GraphPartitionerImplementation> p_impl;
  
};


} // namespace network
} // namespace gridpack


#endif
