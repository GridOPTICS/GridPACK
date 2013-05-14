// -------------------------------------------------------------
/**
 * @file   fields.hpp
 * @author Bruce Palmer
 * @date   April 25, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _fields_h_
#define _fields_h_

#include <vector>
#include <map>
#include "gridpack/parallel/distribution.hpp"
#include "gridpack/network/base_network.hpp"
#include "smart_ptr.hpp"
// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------
namespace gridpack {
namespace network {

template <class elem>
class BaseField  {
  public:
    /**
     * Constructor
     * @param network: Network associated with field
     */
    BaseField(BaseNetwork *network);

    /**
     * Destructor
     */
    ~BaseField(void);

    /**
     * Index operator
     * @param index: index of element in field
     * @return: element at network index
     */
    elem& operator[] (int index);

    /**
     * Return size of field
     */
    int size(void);

    friend class BaseNetwork;
  private:
    /**
     * Add another element to the field
     * @param new_elem: element to be appended to field
     * @return: index of new element
     */
    int append(elem new_elem);

    /**
     * Add another element to the field. Use default element
     * @return: index of new element
     */
    int append(void);

    /**
     * Delete element from field
     * @param index: index of element to be deleted
     * @return: success or failure of delete operation
     */
    bool delete(int index);

    /**
     * Clear all elements from field
     */
    void clear(void);

    BaseNetwork *p_network;
    vector<elem> p_vector;
    elem p_type;
};

// -------------------------------------------------------------
//  class BusField:
//  This class represents fields defined on the network buses
// -------------------------------------------------------------
template <class elem>
class BusField : public BaseField<elem> {
  public:
    /**
     * Constructor
     */
    BusField(BaseNetwork *network);

    /**
     * Destructor
     */
    ~BusField(void);

    /**
     * Get the local indices of the branches that are attached
     * to the bus element indicated by index
     * @param index: index of bus element
     * @return: list of local indices of branches attached to
     *         this bus
     */
    vector<int> neighborBranches(int index);

    /**
     * Get the local indices of the buses that are connected
     * to this bus via a branch element
     * @param index: index of bus element
     * @return: list of local indices of neighboring buses
     */
    vector<int> neighborBuses(int index);

    /**
     * Determine whether or not bus is local or ghost
     * @param: index of bus
     * @return: bus is local (true) or ghost (false)
     */
    bool active(int index);

  private:
    /**
     * Add another element to the field where element is assumed to be a base
     * network component
     * @param new_elem: element to be appended to field
     * @return: index of new element
     */
    int appendComponent(smart_ptr<elem> new_elem);

  friend class BaseNetwork;
};

// -------------------------------------------------------------
//  class BranchField:
//  This class represents fields defined on the network branches
// -------------------------------------------------------------
template <class elem>
class BranchField : public BaseField<elem> {
  public:
    /**
     * Constructor
     */
    BranchField(BaseNetwork *network);

    /**
     * Destructor
     */
    ~BranchField(void);

    /**
     * Get one of the terminal buses on this branch
     * @param index: index of branch element
     * @return: index to the bus at one end of the branch
     */
    int getBus1(int index);

    /**
     * Get the other terminal bus on this branch
     * @param index: index of branch element
     * @return: index to the bus at the other end of the
     *        branch
     */
    int getBus2(int index);

    /**
     * Determine whether or not branch is local or ghost
     * @param: index of bus
     * @return: bus is local (true) or ghost (false)
     */
    bool active(int index);

  friend class BaseNetwork;
};
}  // namespace network
}  // namespace gridpack
#endif
