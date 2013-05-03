#ifndef _se_factory_h_
#define _se_factory_h_

#define BUSES "buses"
#define BRANCHES "branches"
#define RAW_BUS_DATA "raw_bus_data"
#define RAW_BRANCH_DATA "raw_branch_data"

class SEFactory {

public:
  /**
   * Basic constructor
   */
  SEFactory(BaseNetwork *network);

  /**
   * Basic destructor
   */
  ~SEFactory();

  /**
   * Move data from data collection objects on network into network components
   * @return: true of all necessary data fields where found
   */
  bool load(void);

  /**
   * Set state of all network components so that the admittance matrix and H
   * matrix can be built from the network components
   */
  void setState(void);

private:
  
  BaseNetwork *p_network;
};
#endif //_se_factory_h_
