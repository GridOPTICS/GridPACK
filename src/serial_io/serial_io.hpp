/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   serial_io.hpp
 * @author Bruce Palmer
 * @author Arun Veeramany (GOSS)
 * @date   2018-03-16 07:20:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _serial_io_h_
#define _serial_io_h_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/utilities/exception.hpp"
#ifdef USE_GOSS
#include "gridpack/serial_io/goss_client.hpp"
#endif

namespace gridpack {
namespace serial_io {

// -------------------------------------------------------------
// A set of classes to support output of information from buses
// and branches to standard output. Each bus or branch is
// responsible for creating a string that can be written to
// standard out. These modules then organize these sequentially
// and write them from process 0
// -------------------------------------------------------------

template <class _network>
class SerialBusIO {
  public:

  /**
   * Simple constructor
   * @param max_str_len the maximum string length written out by any bus
   * @param network the network for which output is desired
   */
  SerialBusIO(int max_str_len,
              boost::shared_ptr<_network> network)
  {
    p_GAgrp = network->communicator().getGroup();
    p_GA_type = NGA_Register_type(max_str_len);
    p_network = network;
    p_size = max_str_len;
  
    // Find total number of buses in network and create GAs for moving strings
    // around
    int nbus = p_network->totalBuses();
    int one = 1;
    p_stringGA = GA_Create_handle();
    GA_Set_data(p_stringGA,one,&nbus,p_GA_type);
    GA_Set_pgroup(p_stringGA, p_GAgrp);
    GA_Allocate(p_stringGA);
    p_maskGA = GA_Create_handle();
    GA_Set_data(p_maskGA,one,&nbus,C_INT);
    GA_Set_pgroup(p_maskGA, p_GAgrp);
    GA_Allocate(p_maskGA);
//#ifdef USE_GOSS
//    p_goss = NULL;
//    p_channel = false;
//#endif
  }

#ifdef USE_GOSS
  /* Connect to GOSS
   * @param URI e.g. tcp://gridpack2:61616?wireFormat=openwire
   * @param user username to connect to GOSS
   * @param password  password to connect to GOSS
   */
  void connectToGOSS(std::string URI, std::string user, std::string password, std::string topic)
  {
    m_client.connect(URI, user, password); 
    m_topic =  topic;

   ///Just an example for testing
   m_client.publish("goss.request.data.file", "{ \"simulation_id\": \"temp1234\",    \"file_path\": \"ecp_problem3a/rts_contingencies.xml\"}", "temp1234"  );
   std::string file_contents = m_client.subscribeFile("temp1234.reply.goss.request.data.file");
   std::istringstream file_stream(file_contents);
   std::string line;
   while (std::getline(file_stream, line))
        std::cout << line << std::endl;
  }
#endif

  /**
   * Simple Destructor
   */
  ~SerialBusIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    this->close();
  }

  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      this->close();
      p_fout.reset(new std::ofstream);
      p_fout->open(filename);
    }
  }

  /**
   * return IO stream
   * @return IO stream to file
   */
  boost::shared_ptr<std::ofstream> getStream()
  {
    return p_fout;
  }

  /**
   * Set IO stream to point to existing file
   * @param IO stream to file
   */
  void setStream(boost::shared_ptr<std::ofstream> stream)
  {
    p_fout = stream;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_fout) {
        if (p_fout->is_open()) p_fout->close();
      }
    }
    p_fout.reset();
  }

  /**
   * Write output from buses to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_fout) {
      write(*p_fout, signal);
    } else {
      std::cout << "Choosing GOSS for the output stream" << std::endl;
      write(std::cout, signal);
    }
  }

#ifdef USE_GOSS
 /**
   * Send a topic list to GOSS before any IO
   */
  
  void sendTopicList(const char *topics)
  {

    if (m_client.isConnectionValid() && GA_Pgroup_nodeid(p_GAgrp)==0) 
    {
	    std::string topic = "topic/goss/gridpack/topic_list";
	    std::string topic_list(topic);
	    m_client.publish(topic_list, topics);   //topics are space delimited
            //std::cout << "GOSS Acknowledgement: " << m_client.subscribe("topic.goss.gridpack.acknowledge") << std::endl;
    }

  }

  /**
   * Open a channel for IO. This assumes that the application has already
   * specified a complete list of topics using GOSSInit
   * @param topic tag used in publish-subscribe that will identify messages
   * from this program
   * @param URI string for address of broker
   * @param username account name for server recieve messages
   * @param passwd password for server
   */
 /* void openChannel(const char *topic)
  {
    if (!p_goss) p_goss = gridpack::goss::GOSSUtils::instance();
    if (!p_channel && GA_Pgroup_nodeid(p_GAgrp)==0) {
      gridpack::parallel::Communicator comm = p_network->communicator();
      p_goss->openGOSSChannel(comm,topic);
      p_channel = true;
    } else {
      if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
        printf("ERROR: Channel already opened\n");
      }
    }
  }*/

  /**
   * Close IO channel - obsolete
   */
  void closeChannel()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      //gridpack::parallel::Communicator comm = p_network->communicator();
      //p_goss->closeGOSSChannel(comm);
      //p_channel = false;
    }
  }
#endif

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param str character string containing the header
   */
  void header(const char *str)
  {
    if (p_fout) 
    {
      *p_fout << str;
    } 
#ifdef USE_GOSS
    else if (m_client.isConnectionValid()) 
    {
      std::cout << "Choosing GOSS for the output stream" << std::endl;
      m_client.publish(m_topic, str);
    }
#endif
    else
      std::cout << str;


/*
    if (p_fout) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
*/
  }

  /**
   * This is a function that can use the machinery that has been set up in the
   * serialBusIO class to move data to the head node
   * @param vector of data at the head node (this vector will be empty on all
   * other nodes)
   * @signal character string that can be used to specify which data is being
   * requested
   */
  template <class _data_type> void gatherData(std::vector<_data_type>
      &data_vector, const char* signal = NULL)
  {
    int nBus = p_network->numBuses();
    if (sizeof(_data_type) > p_size) {
      char buf[256];
      sprintf(buf,"SerialBusIO::gatherData: data_type size inconsistent"
          " with allocated size: data: %ld allocated: %d\n",
          sizeof(_data_type),p_size);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    _data_type data;
    int dsize = sizeof(_data_type);
    int nwrites = 0;
    int i;
    int one = 1;
    GA_Zero(p_maskGA);

    // Count up total strings being written from this processor
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->getDataItem(&data,signal)) {
        nwrites++;
      }
    }

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBus; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBus(i) &&
            p_network->getBus(i)->getDataItem(ptr,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBusIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr ++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    data_vector.clear();
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              memcpy(&data,ptr,dsize);
              data_vector.push_back(data);
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

#ifdef USE_GOSS
  /**
   * Dump the contents of the channel buffer. Allows users more control over
   * messsages to the GOSS server
   */
  void dumpChannel()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      //printf("Sending message of length %d\n",p_channel_buf.length());
      //p_goss->sendGOSSMessage(p_channel_buf);
    }
  }
#endif

  /**
   * Write output from buses to a vector of strings
   * @param signal an optional character string used to control contents of
   *                output
   * @return vector of strings containing output strings
   */
  std::vector<std::string> writeStrings(const char *signal = NULL)
  {
    int nBus = p_network->numBuses();
    char *string;
    int nwrites = 0;
    int i;
    int one = 1;
    std::vector<std::string> ret;
    string = (char*)malloc(p_size*sizeof(char));

    // Count up total strings being written from this processor
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->serialWrite(string,p_size,signal)) {
        nwrites++;
      }
    }
    free(string);

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf = NULL;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBus; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBus(i) &&
            p_network->getBus(i)->serialWrite(ptr,p_size,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBusIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              std::string tmp = ptr;
              ret.push_back(tmp);
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    return ret;
  }


  protected:

  /**
   * Write output from buses to standard out
   * @param out stream object for output
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(std::ostream & out, const char *signal = NULL)
  {
    int nBus = p_network->numBuses();
    char *string;
    int nwrites = 0;
    int i;
    int one = 1;
    string = (char*)malloc(p_size*sizeof(char));

    // Count up total strings being written from this processor
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->serialWrite(string,p_size,signal)) {
        nwrites++;
      }
    }
    free(string);

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf = NULL;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBus; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBus(i) &&
            p_network->getBus(i)->serialWrite(ptr,p_size,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBusIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
#ifndef USE_GOSS
              out << ptr;
#else
              if (m_client.isConnectionValid()) {
		std::cout << "publishing " << m_topic << " to GOSS" << std::endl;
		m_client.publish(m_topic,ptr);

              /*if (p_channel) {
                p_channel_buf.append(ptr);*/
		
              } else {
		std::cout <<  "goss p_connection is null, writing to stdout" << std::endl;
                out << ptr;
              }
#endif
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param out stream object for output
   * @param str character string containing the header
   */
/*  void header(std::ostream & out, const char *str)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
#ifndef USE_GOSS
      out << str;
#else
      if (p_channel) {
        p_channel_buf.append(str);
      } else {
        out << str;
      }
#endif
    }
  }*/

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
#ifdef USE_GOSS
    gridpack::goss::GOSSClient m_client;
    std::string m_topic;

/*    gridpack::goss::GOSSUtils *p_goss;
    std::string p_channel_buf;
    bool p_channel;*/
#endif
};

template <class _network>
class SerialBranchIO {
  public:
  /**
   * Simple constructor
   * @param max_str_len the maximum string length written out by any branch
   * @param network the network for which output is desired
   */
  SerialBranchIO(int max_str_len,
                 boost::shared_ptr<_network> network)
  {
    p_GAgrp = network->communicator().getGroup();
    p_GA_type = NGA_Register_type(max_str_len);
    p_network = network;
    p_size = max_str_len;

    // Find total number of branches in network and create GAs for moving strings
    // around
    int nbranch = p_network->totalBranches();
    int one = 1;
    p_stringGA = GA_Create_handle();
    GA_Set_data(p_stringGA,one,&nbranch,p_GA_type);
    GA_Set_pgroup(p_stringGA, p_GAgrp);
    GA_Allocate(p_stringGA);
    p_maskGA = GA_Create_handle();
    GA_Set_data(p_maskGA,one,&nbranch,C_INT);
    GA_Set_pgroup(p_maskGA, p_GAgrp);
    GA_Allocate(p_maskGA);
    p_fout.reset();
  }

  /**
   * Simple Destructor
   */
  ~SerialBranchIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    this->close();
  }

#ifdef GOSS
  /* Connect to GOSS
   * @param URI e.g. tcp://gridpack2:61616?wireFormat=openwire
   * @param user username to connect to GOSS
   * @param password  password to connect to GOSS
   */	
  void connectToGOSS(std::string URI, std::string user, std::string password)
  {
    m_client.connect(URI, user, password);
  }


  void sendTopicList(const char *topics)
  {

    if (m_client.isConnectionValid() && GA_Pgroup_nodeid(p_GAgrp)==0)
    {
            std::string topic = "topic/goss/gridpack/topic_list";
            std::string topic_list(topic);
            m_client.publish(topic_list, topics);   //topics are space delimited
            //std::cout << "GOSS Acknowledgement: " << m_client.subscribe("topic.goss.gridpack.acknowledge") << std::endl;
    }

  }
#endif
  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      this->close();
      p_fout.reset(new std::ofstream);
      p_fout->open(filename);
    }
  }

  /**
   * return IO stream
   * @return IO stream to file
   */
  boost::shared_ptr<std::ofstream> getStream()
  {
    return p_fout;
  }

  /**
   * Set IO stream to point to existing file
   * @param IO stream to file
   */
  void setStream(boost::shared_ptr<std::ofstream> stream)
  {
    p_fout = stream;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_fout) {
        if (p_fout->is_open()) p_fout->close();
      }
    }
    p_fout.reset();
  }

  /**
   * Write output from branches to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_fout) {
      write(*p_fout, signal);
    } else {
      write(std::cout, signal);
    }
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param str character string containing the header
   */
  /*void header(const char *str)
  {
    if (p_fout) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
  }*/

  /**
   * This is a function that can use the machinery that has been set up in the
   * serialBranchIO class to move data to the head node
   * @param vector of data at the head node (this vector will be empty on all
   * other nodes)
   * @signal character string that can be used to specify which data is being
   * requested
   */
  template <class _data_type> void gatherData(std::vector<_data_type>
      &data_vector, const char* signal = NULL)
  {
    int nBranch = p_network->numBranches();
    if (sizeof(_data_type) > p_size) {
      char buf[256];
      sprintf(buf,"SerialBranchIO::gatherData: data_type size inconsistent"
          " with allocated size: data: %ld allocated: %d\n",
          sizeof(_data_type),p_size);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    _data_type data;
    int dsize = sizeof(_data_type);
    int nwrites = 0;
    int i;
    int one = 1;
    GA_Zero(p_maskGA);

    // Count up total strings being written from this processor
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->getDataItem(&data,signal)) nwrites++;
    }

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;

    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBranch; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBranch(i) &&
            p_network->getBranch(i)->getDataItem(ptr,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBranchIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    data_vector.clear();
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              memcpy(&data,ptr,dsize);
              data_vector.push_back(data);
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

  /**
   * Write output from branches to a vector of strings
   * @param signal an optional character string used to control contents of
   *                output
   * @return a vector of strings contain output from branches
   */
  std::vector<std::string> writeStrings(const char *signal = NULL)
  {
    int nBranch = p_network->numBranches();
    char *string;
    string = new char[p_size];
    int nwrites = 0;
    int i;
    int one = 1;
    std::vector<std::string> ret;

    // Count up total strings being written from this processor
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->serialWrite(string,p_size,signal)) nwrites++;
    }
    delete [] string;

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBranch; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBranch(i) &&
            p_network->getBranch(i)->serialWrite(ptr,p_size,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBranchIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              std::string tmp = ptr;
              ret.push_back(tmp);
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    return ret;
  }
  protected:

  /**
   * Write output from branches to standard out
   * @param out stream object for output
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(std::ostream & out, const char *signal = NULL)
  {
    int nBranch = p_network->numBranches();
    char *string;
    string = new char[p_size];
    int nwrites = 0;
    int i;
    int one = 1;

    // Count up total strings being written from this processor
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->serialWrite(string,p_size,signal)) nwrites++;
    }
    delete [] string;

    // Set up buffers to scatter strings to global buffer
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      std::vector<int*> index(nwrites);
      std::vector<int> indexbuf(nwrites);
      iptr = &indexbuf[0];
      std::vector<int> ones(nwrites);
      char *strbuf;
      if (nwrites*p_size > 0) strbuf = new char[nwrites*p_size];
      ptr = strbuf;
      int ncnt = 0;
      for (i=0; i<nBranch; i++) {
        if (ncnt >= nwrites) break;
        if (p_network->getActiveBranch(i) &&
            p_network->getBranch(i)->serialWrite(ptr,p_size,signal)) {
          index[ncnt] = iptr;
          *(index[ncnt]) = p_network->getGlobalBranchIndex(i);
          ones[ncnt] = 1;
          ncnt++;
          ptr += p_size;
          iptr++;
        }
      }

      // Scatter data to global buffer and set mask array
      if (ncnt > 0) {
        NGA_Scatter(p_stringGA,strbuf,&index[0],nwrites);
        NGA_Scatter(p_maskGA,&ones[0],&index[0],nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
    }
    GA_Pgroup_sync(p_GAgrp);

    // String data is now stored on global array. Process 0 now retrieves data
    // from each successive processor and writes it to standard out  
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      int nprocs = GA_Pgroup_nnodes(p_GAgrp);
      int lo, hi;
      for (i=0; i<nprocs; i++) {
        NGA_Distribution(p_maskGA, i, &lo, &hi);
        int ld = hi - lo + 1;
        // Figure out how many strings are coming from process i
        std::vector<int> imask(ld);
        NGA_Get(p_maskGA,&lo,&hi,&imask[0],&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char *iobuf;
          if (p_size*nwrites > 0) iobuf = new char[p_size*nwrites];
          std::vector<int*> index(nwrites);
          std::vector<int> indexbuf(nwrites);
          iptr = &indexbuf[0];
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,&index[0],nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              out << ptr;
              ptr += p_size;
              nwrites++;
            }
          }
          if (p_size*nwrites > 0) delete [] iobuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
  }

  /**
   * Write single string to standard output. This is used to write headers for a
   * data listing. It is mostly a convenience function so that users do not have
   * to identify the head node
   * @param out stream object for output
   * @param str character string containing the header
   */
  /*void header(std::ostream & out, const char *str) const
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      out << str;
    }
  }*/
public:

  void header(const char *str)
  {
    if (p_fout)
    {
      *p_fout << str;
    } 
#ifdef USE_GOSS
    else if (m_client.isConnectionValid()) 
    {
      std::cout << "Choosing GOSS for the output stream" << std::endl;
      m_client.publish(m_topic, str);
    }
#endif
    else
      std::cout << str;

  }

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
#ifdef USE_GOSS
    gridpack::goss::GOSSClient m_client;
    std::string m_topic;
#endif
};

}   // serial_io
}   // gridpack
#endif  // _serial_io_h_
