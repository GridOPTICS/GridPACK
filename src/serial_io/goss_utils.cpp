/*
 * Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   goss_utils.cpp
 * @author Bruce Palmer
 * @date   2017-05-24 14:49:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifdef USE_GOSS
#include "gridpack/serial_io/goss_utils.hpp"
#include "gridpack/utilities/string_utils.hpp"

gridpack::goss::GOSSUtils *gridpack::goss::GOSSUtils::p_instance = NULL;

/**
 * Retrieve instance of the GOSSUtils object
 * @return pointer to GOSSUtils object
 */
gridpack::goss::GOSSUtils* gridpack::goss::GOSSUtils::instance()
{
  if (p_instance == NULL) {
    p_instance =  new GOSSUtils();
  }
  return p_instance;
}

/**
 * Initialize goss bus specifying all topics that will be used in this
 * application. This must be called on the world communicator.
 * @param topics list of topics that will be sent by the application
 * @param URI channel URI
 * @param username server username
 * @param passwd server password
 */
void gridpack::goss::GOSSUtils::initGOSS(std::vector<std::string> &topics,
    const char *URI, const char* username, const char *passwd)
{
  char sbuf[128];
  int nvals = topics.size();
  int i;
  // Store channel parameters
  p_URI = URI;
  p_username = username;
  p_passwd = passwd;

  // Concatenate topics into a single string
  if (topics.size() == 0) {
    printf("No topics provided when initializing GOSSUtils\n");
  }
  std::string list = topics[0];
  if (nvals > 0) {
    list = topics[0];
    for (i=1; i<nvals; i++) {
      sprintf(sbuf," %s",topics[i].c_str());
      list.append(sbuf);
    }
  }

  // Send string using topic "topic/goss/gridpack"
  if (GA_Nodeid()==0) {
    char topic[128];

    Connection *connection;
    Session *session;
    Destination *destination;
    MessageProducer *producer;


    std::auto_ptr<ActiveMQConnectionFactory>
      connectionFactory(new ActiveMQConnectionFactory(p_URI)) ;
    // Create a Connection
    connection = connectionFactory->createConnection(p_username, p_passwd);
    connection->start();

    // Create a Session
    session = connection->createSession(Session::AUTO_ACKNOWLEDGE);

    // Create the destination (Topic or Queue)
    sprintf(topic,"topic/goss/gridpack/topic_list");
    //destination = session->createQueue(topic);
    printf("topic = %s\n", topic);
    destination = session->createTopic(topic);


    // Create a MessageProducer from the Session to the Topic
    producer = session->createProducer(destination);
    producer->setDeliveryMode(DeliveryMode::NON_PERSISTENT);

    // Send topics
    std::auto_ptr<TextMessage> message_dat(
        session->createTextMessage(list));
    printf("list: %s\n", list.c_str());
    producer->send(message_dat.get());
    // Close data connection
    // Send final message indicating that channel is being close
    std::string buf = "Closing channel";
    std::auto_ptr<TextMessage>
      end_message(session->createTextMessage(buf));
    producer->send(end_message.get());
    if (connection) delete connection;
    if (session) delete session;
    if (destination) delete destination;
    if (producer) delete producer;
  }
}

/**
 * Open GOSS channel
 * @param comm communicator for whatever module is opening the channel
 * @param topic channel topic. This must be chosen from list of topics
 *              used to initialize the GOSSUtils object
 */
void gridpack::goss::GOSSUtils::openGOSSChannel(gridpack::parallel::Communicator &comm,
    std::string topic)
{
  p_grp = comm.getGroup();
  if (!p_open && GA_Pgroup_nodeid(p_grp) == 0) {
    printf("Opening Channel\n");

    std::auto_ptr<ActiveMQConnectionFactory> connectionFactory(new
        ActiveMQConnectionFactory(p_URI));
    // Create a Connection
    p_connection = connectionFactory->createConnection(p_username,
        p_passwd);
    p_connection->start();

    // Create a Session
    p_session = p_connection->createSession(Session::AUTO_ACKNOWLEDGE);

    // Create the destination (Topic or Queue)
    std::string new_topic("topic.goss.gridpack.");
    //std::string new_topic("topic/goss/gridpack/");
    // Make sure there is no white space around topic
    gridpack::utility::StringUtils utils;
    utils.trim(topic);
    new_topic.append(topic);
    printf("new topic = %s\n", new_topic.c_str());
    p_destination = p_session->createTopic(new_topic);

    // Create a MessageProducer from the Session to the Topic
    p_producer = p_session->createProducer(p_destination);
    p_producer->setDeliveryMode(DeliveryMode::NON_PERSISTENT);

    p_open = true;

    char sbuf[128];
    gridpack::utility::CoarseTimer *timer =
      gridpack::utility::CoarseTimer::instance();
    sprintf(sbuf,"_goss_channel_opened %f",timer->currentTime());
    std::auto_ptr<TextMessage>
      message(p_session->createTextMessage(sbuf));
    p_producer->send(message.get());
  } else {
    if (GA_Pgroup_nodeid(p_grp) == 0)
    {
      printf("ERROR: Channel already opened\n");
    }
  }
}

/**
 * Close GOSS channel.
 */
void gridpack::goss::GOSSUtils::closeGOSSChannel(gridpack::parallel::Communicator &comm)
{
  if (GA_Pgroup_nodeid(p_grp) == 0) {
    // Send final message indicating the channel is being closed
    std::string buf = "_goss_channel_closed";
    std::auto_ptr<TextMessage> message(p_session->createTextMessage(buf));
    p_producer->send(message.get());
    if (p_connection) delete p_connection;
    if (p_session) delete p_session;
    if (p_destination) delete p_destination;
    if (p_producer) delete p_producer;
    p_open = false;
  }
}

/**
 * Send a message over an open GOSS channel
 * @param text message to be sent
 */
void gridpack::goss::GOSSUtils::sendGOSSMessage(std::string &text)
{
  if (GA_Pgroup_nodeid(p_grp) == 0) {
    if (p_open) {
      std::auto_ptr<TextMessage> message(
          p_session->createTextMessage(text));
      printf("Sending message of length %d\n",text.length());
      p_producer->send(message.get());
    } else {
      printf("No GOSS channel is open");
    }
  } 
}

/**
 * Send a message over a GOSS channel
 */
void gridpack::goss::GOSSUtils::sendChannelMessage(const char *topic, const char *URI,
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

/**
 * Shut down GOSS communication
 */
void gridpack::goss::GOSSUtils::terminateGOSS()
{
  if (GA_Nodeid()==0) {
    char topic[128];

    Connection *connection;
    Session *session;
    Destination *destination;
    MessageProducer *producer;


    std::auto_ptr<ActiveMQConnectionFactory>
      connectionFactory(new ActiveMQConnectionFactory(p_URI)) ;
    // Create a Connection
    connection = connectionFactory->createConnection(p_username, p_passwd);
    connection->start();

    // Create a Session
    session = connection->createSession(Session::AUTO_ACKNOWLEDGE);

    // Create the destination (Topic or Queue)
    sprintf(topic,"topic/goss/gridpack/close_goss");
    destination = session->createTopic(topic);

    // Create a MessageProducer from the Session to the Topic
    producer = session->createProducer(destination);
    producer->setDeliveryMode(DeliveryMode::NON_PERSISTENT);

    // Send final message indicating that channel is being close
    std::string buf = "Closing GOSS";
    std::auto_ptr<TextMessage>
      end_message(session->createTextMessage(buf));
    producer->send(end_message.get());
    if (connection) delete connection;
    if (session) delete session;
    if (destination) delete destination;
    if (producer) delete producer;
  }
}

/**
 * Simple constructor
 */
gridpack::goss::GOSSUtils::GOSSUtils()
{
  p_open = false;
  p_connection = NULL;
  p_session = NULL;
  p_destination = NULL;
  p_producer = NULL;
  p_grp = GA_Pgroup_get_world();
}

/**
 * Simple destructor
 */
gridpack::goss::GOSSUtils::~GOSSUtils()
{
}
#endif
