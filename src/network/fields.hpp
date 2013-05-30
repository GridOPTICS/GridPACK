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
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
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
#if 0
    BaseField(BaseNetwork *network)
    {
      p_network = network;
      p_type = elem;
    };
#else
    BaseField()
    {
    };
#endif

    /**
     * Destructor
     */
    ~BaseField(void)
    {
      p_vector.clear();
    };

    /**
     * Index operator
     * @param index: index of element in field
     * @return: element at network index
     */
    elem& operator[] (int index)
    {
      return p_vector[index];
    };

    /**
     * Return size of field
     */
    int size(void);

#if 0
    friend class BaseNetwork;
#endif
    /**
     * Add another element to the field
     * @param new_elem: element to be appended to field
     * @return: index of new element
     */
    int append(elem new_elem)
    {
      p_vector.push_back(new_elem);
      return p_vector.size()-1;
    };

    /**
     * Add another element to the field. Use default element
     * @return: index of new element
     */
    int append(void)
    {
      elem p;
      p_vector.push_back(p);
      return p_vector.size()-1;
    };

    /**
     * Remove last element from field
     */
    void pop_back(void)
    {
      p_vector.pop_back(); 
    }

  protected:
    /**
     * Delete element from field
     * @param index: index of element to be deleted
     * @return: success or failure of delete operation
     */
    bool deleteElement(int index)
    {
      int size = p_vector.size();
      if (index >= size || index < 0) {
        return false;
      }
      if (index == size-1) {
        p_vector.pop_back();
        return true;
      }
      vector<elem>::iterator p;
      p = p_vector.begin();
      for (i=0; i<index; i++) p++;
      p_vector.erase(p);
      return true;
    };

    /**
     * Clear all elements from field
     */
    void clear(void)
    {
      clear();
    };

#if 0
    BaseNetwork *p_network;
#endif
    std::vector<elem> p_vector;
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
#if 0
    BusField(BaseNetwork *network);
#else
    BusField()
    {
    };
#endif

    /**
     * Destructor
     */
    ~BusField(void)
    {
      this->clear();
    };

    /**
     * Get the local indices of the branches that are attached
     * to the bus element indicated by index
     * @param index: index of bus element
     * @return: list of local indices of branches attached to
     *         this bus
     */
    std::vector<int> neighborBranches(int index)
    {
      //TODO: Need some kind of implementation
    };

    /**
     * Get the local indices of the buses that are connected
     * to this bus via a branch element
     * @param index: index of bus element
     * @return: list of local indices of neighboring buses
     */
    std::vector<int> neighborBuses(int index)
    {
      //TODO: Need some kind of implementation
    };

    /**
     * Determine whether or not bus is local or ghost
     * @param: index of bus
     * @return: bus is local (true) or ghost (false)
     */
    bool active(int index)
    {
      //TODO: Need some kind of implementation
    };

#if 0
  friend class BaseNetwork;
#endif
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
#if 0
    BranchField(BaseNetwork *network);
#else
    BranchField()
    {
    };
#endif

    /**
     * Destructor
     */
    ~BranchField(void)
    {
      this->clear();
    };

    /**
     * Get one of the terminal buses on this branch
     * @param index: index of branch element
     * @return: index to the bus at one end of the branch
     */
    int getBus1(int index)
    {
      //TODO: Some kind of implementation
    };

    /**
     * Get the other terminal bus on this branch
     * @param index: index of branch element
     * @return: index to the bus at the other end of the
     *        branch
     */
    int getBus2(int index)
    {
      //TODO: Some kind of implementation
    };

    /**
     * Determine whether or not branch is local or ghost
     * @param: index of bus
     * @return: bus is local (true) or ghost (false)
     */
    bool active(int index)
    {
      //TODO: Some kind of implementation
    };

#if 0
  friend class BaseNetwork;
#endif
};
}  // namespace network
}  // namespace gridpack
#endif
