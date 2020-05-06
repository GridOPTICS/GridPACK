/** Copyright (c) 2013 Battelle Memorial Institute
*     Licensed under modified BSD License. A copy of this license can be found
*     in the LICENSE file in the top level directory of this distribution.
*/

/**
 * @file goss_client.cpp
 * @author Arun Veeramany
 * @brief  This is a wrapper class to publish to/subscribe from GOSS
 * @version 1.0
 * @date 2019-09-13
 * 
 * @copyright Copyright (c) 2019
 * 
 */


#include "gridpack/serial_io/goss_client.hpp"

/**
 * @brief Construct a new GOSSClient object
 *        Invokes connect with supplied connection parameters after object is constructed
 * @param URI A URI e.g. tcp://poorva.pnl.gov:61616?wireFormat=openwire
 * @param username credentials to get connected to GOSS - username
 * @param password credentials to get connected to GOSS - password
 */
gridpack::goss::GOSSClient::GOSSClient(const std::string &URI, const std::string &username, const std::string &password)
{
	connect(URI, username, password);
}

/**
 * @brief Connect to a GOSS server with supplied connection parameters
 * 
 * @param URI A URI e.g. tcp://poorva.pnl.gov:61616?wireFormat=openwire
 * @param username credentials to get connected to GOSS - username
 * @param password credentials to get connected to GOSS - password
 */
void gridpack::goss::GOSSClient::connect(const std::string &URI, const std::string &username, const std::string &password)
{
	activemq::library::ActiveMQCPP::initializeLibrary();
	ActiveMQConnectionFactory *connectionFactory = new ActiveMQConnectionFactory(URI);
	p_connection = connectionFactory->createConnection(username, password);
	delete connectionFactory;
	p_connection->start();
	p_session = p_connection->createSession(Session::AUTO_ACKNOWLEDGE);
	p_message_producer = p_session->createProducer(NULL);

	if (p_connection == NULL || p_session == NULL || p_message_producer == NULL)
		throw std::exception();
}

/**
 * @brief Destroy the GOSSClient::GOSSClient object
 *        Deletes objects and frees up memory used for a connection, session and messaging
 */
gridpack::goss::GOSSClient::~GOSSClient()
{

	try
	{
		p_message_producer->close();
		p_session->close();

		delete p_message_producer;
		delete p_session;
		p_connection->close();
		delete p_connection;

		activemq::library::ActiveMQCPP::shutdownLibrary();
	}
	catch (...)
	{
		std::cout << "Error shutting down ActiveMQ" << std::endl;
	}
}



/**
 * @brief Check if a GOSS connection is valid
 *
 * @return true if connection is valid, false otherwise
 */
bool gridpack::goss::GOSSClient::isConnectionValid()
{
	return p_connection != NULL;
}

/**
 * @brief Publish a message to GOSS Queue
 *
 * @param topic Topic to which message is published
 * @param topic_detail Details of the topic to be published
 * @param response_id An id to use while subscribing to the response
 * @return true  Attempt to publish succeeded
 * @return false Attempt to publish failed
 */
bool gridpack::goss::GOSSClient::publish(const std::string &topic, const std::string &topic_detail, const std::string &response_id)
{

	try
	{
		Destination *destination = p_session->createQueue(topic);
		TextMessage *text_message = p_session->createTextMessage(topic_detail);

		std::string response_id_rev = (response_id == "") ? "" : response_id + ".";

		Destination *destination_header = p_session->createQueue(response_id_rev + "reply.goss.request.data.file");
		text_message->setCMSReplyTo(destination_header);

		p_message_producer->send(destination, text_message);

		delete destination;
		delete destination_header;
		delete text_message;
	}
	catch (...)
	{
		return false;
	}

	return true;
}


/**
 * @brief Subscribe to a message from GOSS Queue
 *
 * @param topic Topic from which message needs to be subscribed
 * @return  subscribed string from GOSS
 */
std::string gridpack::goss::GOSSClient::subscribe(const std::string &topic)
{
	std::string response = "";

	try
	{
		Destination *destination = p_session->createQueue(topic);
		MessageConsumer* consumer = p_session->createConsumer( destination );

	        const ActiveMQBytesMessage* message = (ActiveMQBytesMessage*)consumer->receive();
		std::vector< unsigned char > p = message->getContent();
		response = std::string ( {p.begin(), p.end() } );

		delete message;
		delete consumer;
		delete destination;
	}
	catch(...)
	{
		std::cerr <<  "Error receving message from GOSS" << std::endl;
	}

	return response;
}

/**
 * @brief Extract substring between two strings
 *
 * @param str String from which to extract the substring
 * @param STARTDELIMETER Left delimiting string
 * @param STOPDELIMITER Right delimiting string
 * @return  Substring between the supplied delimeters
 */
std::string gridpack::goss::GOSSClient::getStringBetweenStrings(std::string str, std::string STARTDELIMITER, std::string  STOPDELIMITER)
{
   unsigned first = str.find(STARTDELIMITER);
   unsigned last = str.find(STOPDELIMITER);
   std::string strNew = str.substr (first+STARTDELIMITER.size(),last-first-STOPDELIMITER.size());
   return strNew;
}

/**
 * @brief Subscribe to a file as string from GOSS
 *
 * @param topic Topic from which file needs to be subscribed
 * @return  subscribed file as a string from GOSS. Returned JSON+XML format is malformed, 
 *          so using a substring extract for now.
 */
std::string gridpack::goss::GOSSClient::subscribeFile(const std::string &topic)
{
	std::string wrappedfile = subscribe(topic);
	std::string filecontents = getStringBetweenStrings(wrappedfile, "file_data=","},\"response" );
	return filecontents;
}


/**
 * @brief Subscribe to a file as vector of strings from GOSS
 *
 * @param topic Topic from which file needs to be subscribed
 * @return  subscribed file as a string from GOSS. Returned JSON+XML format is malformed,
 *          so using a substring extract for now.
 */
std::vector<std::string> gridpack::goss::GOSSClient::subscribeFileAsVector(const std::string &topic)
{
        std::string wrappedfile = subscribe(topic);
        std::string filecontents = getStringBetweenStrings(wrappedfile, "file_data=","},\"response" );

        std::vector <std::string> v;
        std::string line;
        std::istringstream sin(filecontents);
        while(getline(sin,line)){
            v.push_back(line);
        }

        return v;
}

