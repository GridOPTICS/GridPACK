/*
 * Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   goss_utils.hpp
 * @author Bruce Palmer
 * @date   2017-05-24 14:49:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _goss_utils_h_
#define _goss_utils_h_

//#define GOSS_DEBUG

#ifdef GOSS_DEBUG
#define USE_GOSS
#endif

#ifdef USE_GOSS
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/utilities/exception.hpp"
#ifndef GOSS_DEBUG
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

namespace gridpack {
namespace goss {
// -------------------------------------------------------------
// A set of classes to support output of information using the
// GOSS communication libraries
// -------------------------------------------------------------

class GOSSUtils {

  public:

  /**
   * Retrieve instance of the GOSSUtils object
   * @return pointer to GOSSUtils object
   */
  static GOSSUtils *instance();

  /**
   * Initialize goss bus specifying all topics that will be used in this
   * application. This must be called on the world communicator.
   * @param topics list of topics that will be sent by the application
   * @param URI channel URI
   * @param username server username
   * @param passwd server password
   */
  void initGOSS(std::vector<std::string> &topics, const char *URI,
      const char* username, const char *passwd);

  /**
   * Open GOSS channel
   * @param comm communicator for whatever module is opening the channel
   * @param topic channel topic. This must be chosen from list of topics
   *              used to initialize the GOSSUtils object
   */
  void openGOSSChannel(gridpack::parallel::Communicator &comm,
      const std::string topic);

  /**
   * Close GOSS channel.
   */
  void closeGOSSChannel(gridpack::parallel::Communicator &comm);

  /**
   * Send a message over an open GOSS channel
   * @param text message to be sent
   */
  void sendGOSSMessage(std::string &text);

  /**
   * Send a message over a GOSS channel
   */
  void sendChannelMessage(const char *topic, const char *URI,
      const char *username, const char *passwd, const char *msg);

  /**
   * Shut down GOSS communication
   */
  void terminateGOSS();

  protected:

  /**
   * Simple constructor
   */
  GOSSUtils();

  /**
   * Simple destructor
   */
  ~GOSSUtils();

  private:

  /**
   * Channel parameters
   */
  std::string p_URI;
  std::string p_username;
  std::string p_passwd;
  std::string p_current_topic;

#ifndef GOSS_DEBUG
  Connection *p_connection;
  Session *p_session;
  Destination *p_destination;
  MessageProducer *p_producer;
#endif
  bool p_open;
  int p_grp;
  /**
   * Pointer to single instance of GOSSUtils object
   */
  static GOSSUtils *p_instance;
};
}
}
#endif
#endif
