/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _no_print_h
#define _no_print_h

// Simple singleton object to keep track of printing status

namespace gridpack{

class NoPrint {
public:
  
  /**
   * Retrieve instance of the NoPrint object
   */
  static NoPrint *instance();

  /**
   * set status of print object
   * @param flag status of print object
   */
  void setStatus(bool flag);

  /**
   * return status of NoPrint object
   * @return true if no print desired
   */
  bool status();

protected:
  /**
   * Constructor
   */
  NoPrint();

  /**
   * Destructor
   */
  ~NoPrint();

private:

  static NoPrint *p_instance;

  bool                p_status;
};


}    // gridpack

#endif // _no_print_h
