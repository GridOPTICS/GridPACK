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

#ifdef USE_GOSS
#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
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

namespace gridpack {
namespace goss {
// -------------------------------------------------------------
// A set of classes to support output of information using the
// GOSS communication libraries
// -------------------------------------------------------------

class GOSSUtils {
  public:

  /**
   * Simple constructor
   */
  GOSSUtils()
  {
  }

  /**
   * Simple destructor
   */
  ~GOSSUtils()
  {
  }

  /**
   * Send a simple message over a GOSS channel
   */
  void sendChannelMessage(const char *topic, const char *URI,
      const char *username, const char *passwd, const char *msg)
  {
    Connection *connection;
    Session *session;
    Destination *destination;
    MessageProducer *producer;

    std::string brokerURI = URI;
    std::auto_ptr<ActiveMQConnectionFactory>
      connectionFactory(new ActiveMQConnectionFactory(brokerURI)) ;
    // Create a Connection
    std::string User = username;
    std::string Pass = passwd;
    connection = connectionFactory->createConnection(User, Pass);
    connection->start();
    // Create a Session
    session = connection->createSession(Session::AUTO_ACKNOWLEDGE);
    // Create the destination (Topic or Queue)
    destination = session->createTopic(topic);
    // Create a MessageProducer from the Session to the Topic
    producer = session->createProducer(destination);
    producer->setDeliveryMode(DeliveryMode::NON_PERSISTENT);

    std::auto_ptr<TextMessage> message(session->createTextMessage(msg));
    producer->send(message.get());
    // Clean up
    delete connection;
    delete session;
    delete destination;
    delete producer;
  }
};
}
}
#endif
#endif
