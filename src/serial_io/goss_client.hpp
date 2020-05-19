/** Copyright (c) 2013 Battelle Memorial Institute
*     Licensed under modified BSD License. A copy of this license can be found
*     in the LICENSE file in the top level directory of this distribution.
*/

/**
 * @file goss_client.hpp
 * @author Arun Veeramany
 * @brief Header for wrapper around publishing to / subscribing from GOSS
 * @version 1.0
 * @date 2019-09-13
 *
 * @copyright Copyright (c) 2019
 *
 */


#ifndef _goss_client_h_
#define _goss_client_h_

#ifdef GOSS_DEBUG
#define USE_GOSS
#endif

#include<vector>
#include <sstream>

#include "activemq/library/ActiveMQCPP.h"
#include <decaf/lang/Thread.h>
#include <decaf/lang/Runnable.h>
#include <decaf/util/concurrent/CountDownLatch.h>
#include <decaf/lang/Integer.h>
#include <decaf/lang/Long.h>
#include <decaf/lang/System.h>
#include "activemq/core/ActiveMQConnectionFactory.h"
#include "activemq/util/Config.h"
#include <cms/Connection.h>
#include <cms/Session.h>
#include <cms/TextMessage.h>
#include <cms/BytesMessage.h>
#include <cms/MapMessage.h>
#include <cms/ExceptionListener.h>
#include <cms/MessageListener.h>
#include "activemq/commands/ActiveMQBytesMessage.h"

using namespace activemq::core;
using namespace activemq::commands;
using namespace decaf::util::concurrent;
using namespace decaf::util;
using namespace decaf::lang;
using namespace cms;


//gridpack - goss integration namespaces
namespace gridpack {
namespace goss {


/**
 * @brief Wrapper class to publish to / subscribe from GOSS
 * 
 */
class GOSSClient
{

protected:
	//These are freed up in the destructor
	Connection *p_connection = NULL;
	Session *p_session = NULL;
	MessageProducer *p_message_producer = NULL;

public:
	//Constructors and Destructor
	GOSSClient(const std::string &URI, const std::string &username, const std::string &password);
	GOSSClient() {}
	~GOSSClient();

protected:
	//Utility member functions
	std::string getStringBetweenStrings(std::string str, std::string STARTDELIMITER, std::string  STOPDELIMITER);

public:
	//Members to get connected, publish and subscribe
	void connect(const std::string &URI, const std::string &username, const std::string &password);
	bool publish(const std::string &topic, const std::string &topic_detail, const std::string &response_id="");
	bool isConnectionValid();
	std::string subscribe(const std::string &topic);
	std::string subscribeFile(const std::string &topic);
        std::vector<std::string> subscribeFileAsVector(const std::string &topic);


};

}
}

#endif
