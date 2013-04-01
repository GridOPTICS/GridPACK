// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   network.h
 * @author William A. Perkins
 * @date   Fri Mar 22 11:32:28 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 22, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _network_h_
#define _network_h_

#include <vector>
#include "gridpack/parallel/distribution.hpp"
#include "gridpack/network/component.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class NetworkInterface
// -------------------------------------------------------------
class NetworkInterface {
public:

  /// Default constructor.
  NetworkInterface(void);

  /// Destructor
  ~NetworkInterface(void);

protected:

  /// Protected copy constructor to avoid unwanted copies.
  NetworkInterface(const NetworkInterface& old);

};



// -------------------------------------------------------------
//  class Network
// -------------------------------------------------------------
/// A generic representation of a power network
/**
 * This serves as the base class for all power network
 * representations.  It is essentially a container for network
 * components.
 * 
 */
template <typename Bus_, 
          typename Branch_, 
          typename Generator_,
          typename Load>
class Network : 
    public NetworkInterface, 
    public parallel::Distributed {
public:

  /// The kind of bus in this kind of network (must be subclass of BusInterface)
  typedef Bus_ BusType;

  /// A container for buses
  typedef std::list<BusType> BusList;

  /// Iterator types for bus the bus container
  typedef std::list<BusType>::iterator BusIterator;
  typedef std::list<BusType>::const_iterator ConstBusIterator;

  /// The kind of bus in this kind of network (must be subclass of BranchInterface)
  typedef Branch_ BranchType;

  /// A container for buses
  typedef std::list<BranchType> BranchList;

  /// Iterator types for bus the bus container
  typedef std::list<BranchType>::iterator BranchIterator;
  typedef std::list<BranchType>::const_iterator ConstBranchIterator;

  // [More declarations for other component types that need to be contained]

  /// Default constructor.
  explicit Network(const parallel::distribution& dist);

  /// Destructor
  ~Network(void);

  /// Iterator for local buses
  BusIterator local_bus_begin(void);
  BusIterator local_bus_end(void);

  ConstBusIterator local_bus_begin(void) const;
  ConstBusIterator local_bus_end(void) const;

  // [ other exposed iterators ]
  
  /// A way to redistribute the network
  void redistribute(void);

  

protected:

  /// Group of buses assigned to the local process
  BusList local_buses_;

  /// Copies of busese assigned to other processes, but connected to local branches
  BusList ghost_buses_;

  /// Group of branches assigned to the local process
  BranchList local_branches_;

  /// Copies of branches assigned to other processes
  BranchList ghost_branches_;
  
  // [Similar declarations for other component lists (may only need local lists)

private:

  /// Protected copy constructor to avoid unwanted copies.
  Network(const Network& old);

};



#endif
