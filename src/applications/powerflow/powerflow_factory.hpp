#define BUSES "buses"
#define BRANCHES "branches"
#define RAW_BUS_DATA "raw_bus_data"
#define RAW_BRANCH_DATA "raw_branch_data"

class PowerflowFactory {

public:
  /**
   * Basic constructor
   */
  PowerflowFactory(BaseNetwork *network);

  /**
   * Basic destructor
   */
  ~PowerflowFactory();

  /**
   * Move data from data collection objects on network into network components
   * @return: true of all necessary data fields where found
   */
  bool load(void);

  /**
   * Set state of all network components so that the admittance matrix can be
   * built from the network components
   */
  void setState(void);
};
