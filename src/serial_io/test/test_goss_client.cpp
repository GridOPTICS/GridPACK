/**
 * @file test_goss_client.cpp
 * @author your name (arun.veeramany@pnnl.gov)
 * @brief 
 * @version 0.1
 * @date 2019-09-13
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "gridpack/serial_io/goss_client.hpp"

using namespace gridpack::goss;


int main()
{

   std::cout << "This is goss client test." << std::endl;

   std::string username = "system";
   std::string password = "manager";
   std::string URI = "tcp://constance.pnl.gov:61616?wireFormat=openwire"; //openwire (61616), stomp (61614)

   //std::string URI = "tcp://poorva.pnl.gov:61616?wireFormat=openwire"; //openwire (61616), stomp (61614)
   //std::string URI = "tcp://gridpack2:61616?wireFormat=openwire"; //openwire (61616)U, stomp (61614)

   std::cout << "Attempting to connect to " << URI << std::endl;


   //Step 1: Get authenticated
   GOSSClient client(URI, username, password);
   std::cout << "GOSS authentication is complete. Connected to " << URI << std::endl;

   //Step 2: Publish
   std::string topic_list = "topic/goss/gridpack/topic_list";
   std::string topics = "simulation1 simulation2";
   client.publish(topic_list, topics);

   std::string topic2 = "topic.goss.gridpack.simulation1";
   std::string topic_detail2 = "1 2 3 4 5";
   client.publish(topic2, topic_detail2);


   std::string topic3 = "topic.goss.gridpack.simulation3";
   std::string topic_detail3 = "1 2 3 4 5 6 7 8";
   client.publish(topic3, topic_detail3);


   //Step 3: Subscribe
   client.publish("goss.request.data.file", "{ \"simulation_id\": \"temp1234\",    \"file_path\": \"ecp_problem3a/rts_contingencies.xml\"}", "temp1234"  );
   std::vector<std::string> file_contents = client.subscribeFileAsVector("temp1234.reply.goss.request.data.file");
 
   for (std::vector<std::string>::iterator line = file_contents.begin(); line != file_contents.end(); ++line)
	std::cout << *line << std::endl;

   return 0;
}
