// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   parmetis_graph_wrapper.hpp
 * @author William A. Perkins
 * @date   2013-07-09 10:29:42 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _parmetis_graph_wrapper_hpp_
#define _parmetis_graph_wrapper_hpp_

#include <parmetis.h>
#include <vector>
#include <boost/scoped_ptr.hpp>
#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/parallel/distributed.hpp"
#include "adjacency_list.hpp"

namespace GA {
class GlobalArray;
}

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class ParMETISGraphWrapper
// -------------------------------------------------------------
/// Enapsulation of a ParMETIS format graph
/**
 * The main reason for this class is to be able to have a ParMETIS
 * representation of a graph whose distribution is different than the
 * actual, current distribution (as represented by an AdjacencyList
 * instance).
 *
 * ParMETIS requires that at least one of the graph nodes exists on
 * each processor.  This class, hopefully, will be used to overcome
 * that limitation, so, for examples, graphs that exist on only one
 * processor can be partitioned and redistributed
 * 
 * In addition, AdjacencyList works with arbitrarily ordered global
 * node indices.  ParMETIS requires special ordering of nodes. This
 * class is used to translate between the two indexing systems.  In
 * both cases, the global node indices must start with zero.
 */
class ParMETISGraphWrapper 
  : public parallel::Distributed,
    private utility::Uncopyable 
{
public:

  /// Default constructor.
  explicit ParMETISGraphWrapper(const AdjacencyList& alist);

  /// Destructor
  ~ParMETISGraphWrapper(void);

  /// Get the local part of the "Distributed CSR graph" (used by ParMETIS)
  void get_csr_local(std::vector<idx_t>& vtxdist,
                     std::vector<idx_t>& xadj,
                     std::vector<idx_t>& adjncy) const;

  /// Assign partition number for local ParMETIS graph nodes
  void set_partition(const std::vector<idx_t>& vtxdist, 
                     const std::vector<idx_t>& part);

  /// Get partition number for local AdjacencyList nodes
  void get_partition(AdjacencyList::IndexVector& part) const;

protected:

  /// The adjacency list upon which this graph is based
  const AdjacencyList& p_adjacency;

  /// The total number of nodes involved
  int p_global_nodes;

  /// The total number of edges
  int p_global_edges;

  /// A GA to hold the original node data
  /**
   * This is a 2D GA. It's used to hold several things that need to be
   * remembered about the graph nodes: global node id (j=0), initial
   * owner process(j=1), destination process (j=2)
   * 
   */
  boost::scoped_ptr<GA::GlobalArray> p_node_data;

  /// A GA to hold the "local" node id indexed with the original global id
  /**
   * This is used translate between the node indexes as provided by
   * ::p_adjacency to indexes ParMETIS can understand.
   *
   * Each processor may contribute some node indexes. These may be in
   * any order.  Say the 14th node index (13, zero-based) in the
   * global list is 3. The contents of this array at index 3 would be
   * 13.
   *
   */
  boost::scoped_ptr<GA::GlobalArray> p_local_node_id;

  /// The the low index into ::p_node_data for this process's data
  int p_node_lo;

  /// The high index into ::p_node_data for this process's data
  int p_node_hi;
  
  // The GA arrays look like the ParMETIS "Serial CSR" graph format

  /// The graph node index into p_adjncy_ga
  /**
   * The ::p_adjncy_gbl and ::p_xadj_gbl arrays are supposed to look
   * just like what is called the "Serial CSR" graph format in
   * ParMETIS.  
   * 
   */
  boost::scoped_ptr<GA::GlobalArray> p_xadj_gbl;

  /// The global node adjacency list
  /**
   * The ::p_adjncy_gbl and ::p_xadj_gbl arrays are supposed to look
   * just like what is called the "Serial CSR" graph format in
   * ParMETIS.  
   * 
   */
  boost::scoped_ptr<GA::GlobalArray> p_adjncy_gbl;

  /// The initialize routine
  void p_initialize(void);

  /// Routine to initialize the GA part of this instance
  void p_initialize_gbl(const int& gblnodes, const int& locnodes,
                        const int& gbledges, const int& locedges);

};

} // namespace networkx
} // namespace gridpack

#endif
