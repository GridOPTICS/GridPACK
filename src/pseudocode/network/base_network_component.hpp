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

#include "gridpack/parallel/distribution.hpp"
// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------
class BaseNetworkComponent
  : public MatVecInterface {
  public:
    /**
     * Constructor
     */
    BaseNetworkComponent(void);

    /**
     * Destructor
     */
    ~BaseNetworkComponent(void);

    /**
     * Return the size of the component for use in packing and
     * unpacking routines. This might not be needed, but throw
     * it in for now.
     * @return: size of network component
     */
    int Size(void);

  private:

};

// TODO: Might want to put MatrixIndices and VectorIndex operations into a
//       separate class since these can probably be implemented once for all
//       network components

class MatVecInterface {
  public:
    /**
     * Constructor
     */
    virtual MatVecInterface(void);

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
    virtual bool MatrixIndices(int *idx, int *jdx) const;

    /**
     * Return size of matrix block contributed by component
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool MatrixSize(int *isize, int *jsize) const;

    /**
     * Return the values of the matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool MatrixValues(void *values);

    /**
     * Provide the vector index for the network component based
     * on location in the network topology. Return false if this
     * component does not contribute vector elements
     * @param idx: vector index
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool VectorIndex(int *idx) const;

    /**
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool VectorSize(int *isize) const;

    /**
     * Return the values of the vector block.
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool VectorValues(void *values);

    //TODO: May need to include routines that support moving values from vectors
    //      back into network components.

  private:

};

#endif
