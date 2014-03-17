/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   parmetis_graph_wrapper.cpp
 * @author William A. Perkins
 * @date   2014-02-13 09:28:27 d3g096
 * 
 * @brief  Implementation of the ParMETISGraphWrapper class
 * 
 * 
 */
// -------------------------------------------------------------


#include <ga++.h>
#include <boost/mpi/collectives.hpp>
#include <boost/assert.hpp>
#include <boost/lambda/lambda.hpp>

namespace bl = boost::lambda;

#include "parmetis/parmetis_graph_wrapper.hpp"

// define some constants for readability 
static const int one(1);
static const int two(2);

static const int num_node_data(3);

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class ParMETISGraphWrapper
// -------------------------------------------------------------

// -------------------------------------------------------------
// ParMETISGraphWrapper:: constructors / destructor
// -------------------------------------------------------------
ParMETISGraphWrapper::ParMETISGraphWrapper(const AdjacencyList& alist)
  : parallel::Distributed(alist.communicator()),
    utility::Uncopyable(),
    p_adjacency(alist), 
    p_global_nodes(0), p_global_edges(0),
    p_node_data(), p_local_node_id(), 
    p_node_lo(-1), p_node_hi(-1), 
    p_xadj_gbl(), p_adjncy_gbl()
{
  p_initialize();
}

ParMETISGraphWrapper::~ParMETISGraphWrapper(void)
{
}

// -------------------------------------------------------------
// ParMETISGraphWrapper::p_initialize_gbl
// -------------------------------------------------------------
/** 
 * This routine sets up the Global Arrays used by this instance.  
 * 
 * @param gblnodes 
 * @param locnodes 
 * @param gbledges 
 * @param locedges 
 */
void
ParMETISGraphWrapper::p_initialize_gbl(const int& gblnodes, const int& locnodes,
                                       const int& gbledges, const int& locedges)
{
  int maxdim(two);
  int dims[maxdim], lo[maxdim], hi[maxdim], ld[maxdim];
  int tmp[maxdim];
  int sum;
  ld[0] = 1;
  ld[1] = 1;

  int theGAgroup(communicator().getGroup());
  int oldGAgroup = GA_Pgroup_get_default();
  GA_Pgroup_set_default(theGAgroup);

  std::vector<int> ndata(locnodes);

                                // make the node data array (2D)

  dims[0] = gblnodes;  dims[1] = num_node_data;
  p_node_data.reset(new GA::GlobalArray(MT_C_INT, two, dims, 
                                        "ParMETIS Wrapper Node Data", NULL));
  // p_node_data->printDistribution();
  p_node_data->zero();

  // put the global node indexes or ids

  if (locnodes > 0) {

    for (int n = 0; n < locnodes; ++n) {
      ndata[n] = p_adjacency.node_index(n);
    }

    lo[0] = p_node_lo; lo[1] = 0;
    hi[0] = p_node_hi; hi[1] = 0;
    p_node_data->put(lo, hi, &ndata[0], ld);

    // put the node owner process
  
    std::fill(ndata.begin(), ndata.end(), this->processor_rank());
    lo[0] = p_node_lo; lo[1] = 1;
    hi[1] = p_node_hi; hi[1] = 1;
    p_node_data->put(lo, hi, &ndata[0], ld);

  }

  communicator().sync();

                                // make the node id table

  dims[0] = gblnodes;
  dims[1] = two;
  p_local_node_id.reset(new GA::GlobalArray(MT_C_INT, one, dims,
                                            "Node index translator", NULL));

  // p_local_node_id->printDistribution();

  if (locnodes > 0) {
    lo[0] = p_node_lo; lo[1] = 0;
    hi[1] = p_node_hi; hi[1] = 0;
    p_node_data->get(lo, hi, &ndata[0], ld);

    std::vector<int> nidx(locnodes);
    std::vector<int> junk(ndata);
    std::vector<int*> junkidx(locnodes);
    sum = p_node_lo;
    std::vector<int>::iterator i(nidx.begin());
    std::vector<int>::iterator j(junk.begin());
    std::vector<int*>::iterator ip(junkidx.begin());

    for (; i != nidx.end(); ++i, ++j, ++ip, ++sum) {
      *i = sum;
      *ip = &(*j);
    }
    
    lo[0] = p_node_lo; lo[1] = 2;
    hi[1] = p_node_hi; hi[1] = 2;
    p_node_data->put(lo, hi, &nidx[0], ld);

    p_local_node_id->scatter(&nidx[0], &junkidx[0], nidx.size());
  }

  communicator().sync();

  // p_node_data->print();

  // p_local_node_id->print();

                                // make the adjacency list and its index

  dims[0] = gblnodes+1; 
  p_xadj_gbl.reset(new GA::GlobalArray(MT_C_INT, one, dims, 
                                       "ParMETIS Adjacency Index", NULL));
  p_xadj_gbl->zero();

  dims[0] = 2*gbledges; 
  p_adjncy_gbl.reset(new GA::GlobalArray(MT_C_INT, one, dims,
                                         "ParMETIS Adjacency List", NULL));
  p_adjncy_gbl->zero();

  std::vector<AdjacencyList::Index> nbrs;
  std::vector<int> inbrs;
  for (int p = 0; p < this->processor_size(); ++p) {
    if (p == this->processor_rank()) {
      if (locnodes > 0) {
	lo[0] = p_node_lo;
	hi[0] = p_node_lo;
	p_xadj_gbl->get(lo, hi, &tmp[0], ld);
	for (int i = 0; i < locnodes; ++i) {
	  nbrs.clear();
	  p_adjacency.node_neighbors(i, nbrs);
	  inbrs.clear();
	  std::copy(nbrs.begin(), nbrs.end(), std::back_inserter(inbrs));

	  lo[0] = tmp[0];
	  hi[0] = tmp[0] + inbrs.size() - 1;
	  if (hi[0] >= lo[0]) p_adjncy_gbl->put(lo, hi, &inbrs[0], ld);

	  int idx(p_node_lo + i + 1);
	  tmp[0] += inbrs.size();
	  lo[0] = idx;
	  hi[0] = idx;
	  p_xadj_gbl->put(lo, hi, &tmp[0], ld);

	}
      }
    }
    communicator().sync();
  }

  // p_xadj_gbl->print();
  // p_adjncy_gbl->print();

  GA_Pgroup_set_default(oldGAgroup);
}
  

// -------------------------------------------------------------
// ParMETISGraphWrapper::p_initialize
// -------------------------------------------------------------
void
ParMETISGraphWrapper::p_initialize(void)
{

  int localnodes(p_adjacency.nodes());
  p_global_nodes = 0;
  { 
    std::vector<int> nodes_by_proc;
    all_gather(this->communicator().getCommunicator(), 
               localnodes, nodes_by_proc);
    
    std::vector<int> oldnodedist(nodes_by_proc.size() + 1);
    for (size_t i = 0; i < nodes_by_proc.size(); ++i) {
      oldnodedist[i] = p_global_nodes;
      p_global_nodes += nodes_by_proc[i];
    }
    oldnodedist.back() = p_global_nodes;
    
    p_node_lo = oldnodedist[this->processor_rank()];
    p_node_hi = oldnodedist[this->processor_rank() + 1] - 1;
  }

  int localedges(p_adjacency.edges());
  p_global_edges = 0;

  all_reduce(this->communicator().getCommunicator(), 
             localedges, p_global_edges, std::plus<int>());

  p_initialize_gbl(p_global_nodes, localnodes, p_global_edges, localedges);

}

// -------------------------------------------------------------
// ParMETISGraphWrapper::get_csr_local
// -------------------------------------------------------------
void
ParMETISGraphWrapper::get_csr_local(std::vector<idx_t>& vtxdist,
                                    std::vector<idx_t>& xadj,
                                    std::vector<idx_t>& adjncy) const
{
  BOOST_ASSERT(p_node_data);
  BOOST_ASSERT(p_local_node_id);
  BOOST_ASSERT(p_xadj_gbl);
  BOOST_ASSERT(p_adjncy_gbl);
  BOOST_ASSERT(p_global_nodes > 0);
  BOOST_ASSERT(p_global_edges > 0);

  int nproc(this->processor_size());
  int me(this->processor_rank());

                                // build the node distribution vector

  vtxdist.clear();
  vtxdist.resize(nproc+1, 0);
  for (int i = 0; i < p_global_nodes; ++i) {
    int p(i % nproc);
    vtxdist[p] += 1;
  }
  int localnodes(vtxdist[me]);
  int sum(0);
  for (int p = 0; p < nproc; ++p) {
    int tmp(vtxdist[p]);
    vtxdist[p] = sum;
    sum += tmp;
  }
  vtxdist[nproc] = sum;

                                // extract adjacency index

  int maxdim(two);
  int lo[maxdim], hi[maxdim], ld[maxdim];
  ld[0] = 1;
  ld[1] = 1;

  std::vector<int> tmp(localnodes+1);
  lo[0] = vtxdist[me];
  hi[0] = vtxdist[me+1];
  p_xadj_gbl->get(lo, hi, &tmp[0], ld);

  xadj.clear();
  xadj.reserve(tmp.size());
  std::copy(tmp.begin(), tmp.end(), std::back_inserter(xadj));

  sum = xadj.front();
  std::for_each(xadj.begin(), xadj.end(), bl::_1 -= bl::var(sum));
  
  lo[0] = tmp.front();
  hi[0] = tmp.back() - 1;
  int nidxsize(hi[0] - lo[0] + 1);
  std::vector<int> nidx(nidxsize);
  p_adjncy_gbl->get(lo, hi, &nidx[0], ld);

  {  
    std::vector<int*> junkidx(nidx.size());
    std::vector<int>::iterator i(nidx.begin());
    std::vector<int*>::iterator ip(junkidx.begin());

    for (; i != nidx.end(); ++i, ++ip) {
      *ip = &(*i);
    }
    tmp.resize(nidx.size());
    p_local_node_id->gather(&tmp[0], &junkidx[0], tmp.size());
  }


  adjncy.clear();
  adjncy.reserve(tmp.size());
  std::copy(tmp.begin(), tmp.end(), std::back_inserter(adjncy));
  communicator().sync();
}

// -------------------------------------------------------------
// ParMETISGraphWrapper::set_partition
// -------------------------------------------------------------
/** 
 * The caller has used the ParMETIS graph produced by
 * ::get_csr_local() to partition the graph.  
 * 
 * @param vtxdist ParMETIS graph node distribution (from ::get_csr_local)
 * @param part ParMETIS graph node partition assignments
 */
void
ParMETISGraphWrapper::set_partition(const std::vector<idx_t>& vtxdist, 
                                    const std::vector<idx_t>& part)
{
  int me(this->processor_rank());
  int lo[2], hi[2], ld[2];
  lo[0] = vtxdist[me]; lo[1] = 2;
  hi[0] = vtxdist[me+1]-1; hi[1] = 2;
  ld[0] = 1; ld[1] = 1;
  
  BOOST_ASSERT(p_node_data);
  BOOST_ASSERT(part.size() == (hi[0] - lo[0] + 1));

  // idx_t may not be same as int

  std::vector<int> tmp(part.size());
  std::copy(part.begin(), part.end(), tmp.begin());

  p_node_data->put(lo, hi, &tmp[0], ld);
  communicator().sync();
  // p_node_data->print();
}
// -------------------------------------------------------------
// ParMETISGraphWrapper::get_partition
// -------------------------------------------------------------
void
ParMETISGraphWrapper::get_partition(AdjacencyList::IndexVector& part) const
{

  part.clear();
  if (p_node_hi >= p_node_lo) {
    int lo[2], hi[2], ld[2];
    lo[0] = p_node_lo; lo[1] = 2;
    hi[0] = p_node_hi; hi[1] = 2;
    ld[0] = 1; ld[1] = 1;
  
    std::vector<int> tmp(p_node_hi - p_node_lo + 1);
    p_node_data->get(lo, hi, &tmp[0], ld);

    part.reserve(tmp.size());
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(part));
  }
  communicator().sync();
}

} // namespace network
} // namespace gridpack
