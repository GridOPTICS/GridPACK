// -------------------------------------------------------------
/**
 * @file   base_network_component.hpp
 * @author Bruce Palmer
 * @date   April 8, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_network_component_h_
#define _base_network_component_h_

#include "gridpack/parallel/distributed.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
// TODO: Might want to put MatrixIndices and VectorIndex operations into a
//       separate class since these can probably be implemented once for all
//       network components
namespace gridpack{
namespace component{

class MatVecInterface {
  public:
    /**
     * Constructor
     */
    MatVecInterface(void);

    /**
     * Destructor
     */
    virtual ~MatVecInterface(void);

    /**
     * Provide the matrix indices for the network component based
     * on location in the network topology. Return false if this
     * component does not contribute to matrix elements
     * @param idx, jdx: row and column indices of matrix
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixIndices(int *idx, int *jdx) const;

    /**
     * Return size of matrix block contributed by component
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixSize(int *isize, int *jsize) const;

    /**
     * Return the values of the matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixValues(void *values);

    /**
     * Provide the vector index for the network component based
     * on location in the network topology. Return false if this
     * component does not contribute vector elements
     * @param idx: vector index
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorIndex(int *idx) const;

    /**
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block.
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorValues(void *values);

    //TODO: May need to include routines that support moving values from vectors
    //      back into network components.

  private:

};

// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------
class BaseNetworkComponent
  : public MatVecInterface {
  public:
    /**
     * Simple constructor
     */
    BaseNetworkComponent();

    /**
     * Constructor
     * @param network: pointer to network the component is associated with
     * @param idx: local bus or branch index that network is associated with
     */
#if 0
    virtual BaseNetworkComponent(boost::smart_ptr<BaseNetwork> network, int idx);
#endif

    /**
     * Constructor without local network index
     * @param network: pointer to network the component is associated with
     */
#if 0
    virtual BaseNetworkComponent(boost::smart_ptr<BaseNetwork> network);
#endif

    /**
     * Destructor
     */
    virtual ~BaseNetworkComponent(void);

    /**
     * Set the network associated with the component
     */
#if 0
    virtual void setNetwork(boost::smart_ptr<BaseNetwork> network);
#endif

    /**
     * Set the value of the local network index the component is associated with
     * @param idx: value of local network index
     */
    virtual void setIndex(int idx);

    /**
     * Get the value of the local network index the component is associated with
     * @return: value of local network index
     */
    virtual int getIndex();

    /**
     * Return the size of the component for use in packing and
     * unpacking routines. This might not be needed, but throw
     * it in for now.
     * @return: size of network component
     */
    virtual int size(void) const;

  private:
#if 0
    stlplus::smart_ptr<BaseNetwork> p_network;
#endif
    int p_idx;

};

}    // component
}    // gridpack
#endif
