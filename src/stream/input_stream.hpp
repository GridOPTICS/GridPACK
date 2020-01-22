/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _input_stream_h
#define _input_stream_h

#include <string>
#include <fstream>

#ifdef USE_GOSS
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/utilities/exception.hpp"
#include <activemq/library/ActiveMQCPP.h>
#include <decaf/lang/Thread.h>
#include <decaf/lang/Runnable.h>
#include <decaf/util/concurrent/CountDownLatch.h>
#include <decaf/lang/Integer.h>
#include <decaf/lang/Long.h>
#include <decaf/lang/System.h>
#include <activemq/core/ActiveMQConnectionFactory.h>
#include <activemq/util/Config.h>
#include <cms/Connection.h>
#include <cms/Session.h>
#include <cms/TextMessage.h>
#include <cms/BytesMessage.h>
#include <cms/MapMessage.h>
#include <cms/ExceptionListener.h>
#include <cms/MessageListener.h>
using namespace activemq::core;
using namespace decaf::util::concurrent;
using namespace decaf::util;
using namespace decaf::lang;
using namespace cms;
#endif

// Simple class to open a stream of data and read lines

namespace gridpack{
namespace stream{

class InputStream {
public:

  /**
   * Constructor
   */
  InputStream();

  /**
   * Destructor
   */
  ~InputStream();

  /**
   * Open a file
   * @param filename name of file to open
   * @return true of file is successfully opened
   */
  bool openFile(std::string file);

#ifdef USE_GOSS
  /**
   * Open GOSS channel
   * @param topic topic that that is used by server
   * @param URI channel URI
   * @param username server username
   * @param passwd server password
   */
  bool openChannel(const char *topic, const char *URI, const char *username,
      const char *passwd);
#endif

  /**
   * Close a file or other input stream
   */
  void close();

  /**
   * Get next line from streaming source
   * @param line next line from streaming source
   * @return true if another line is found
   */
  bool nextLine(std::string &line);
  
  /**
   * Report if a stream is currently open
   * @return true if stream is open
   */
  bool isOpen();

private:

  std::ifstream p_fout;

  bool p_srcFile;

  bool p_srcChannel;

  bool p_isOpen;

};


}    // stream
}    // gridpack


#endif // _input_stream_h
