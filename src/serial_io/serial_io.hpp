/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   serial_io.hpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:49:01 d3g096
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
    p_useFile = false;
#ifdef USE_GOSS
    p_channel = false;
    p_connection = NULL;
    p_session = NULL;
    p_destination = NULL;
    p_producer = NULL;
#endif
  }

  /**
   * Simple Destructor
   */
  ~SerialBusIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
      if (!p_useFile) p_fout.reset(new std::ofstream);
      p_fout->open(filename);
      p_useFile = true;
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
    p_useFile = true;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_fout->is_open()) p_fout->close();
    }
    p_useFile = false;
  }

  /**
   * Write output from buses to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_useFile) {
      write(*p_fout, signal);
    } else {
      write(std::cout, signal);
    }
  }

#ifdef USE_GOSS
  /**
   * Open a channel for IO
   * @param topic tag used in publish-subscribe that will identify messages
   * from this program
   * @param URI string for address of broker
   * @param username account name for server recieve messages
   * @param passwd password for server
   */
  void openChannel(const char *topic, const char *URI,
      const char *username, const char *passwd)
  {
    if (!p_channel && GA_Pgroup_nodeid(p_GAgrp)==0) {
      printf("Opening Channel\n");
      std::string brokerURI = URI;
      //std::auto_ptr<ConnectionFactory> connectionFactory(
      //ConnectionFactory::createCMSConnectionFactory(brokerURI));

      std::auto_ptr<ActiveMQConnectionFactory>
        connectionFactory(new ActiveMQConnectionFactory(brokerURI)) ;
      // Create a Connection
      std::string User = username;
      std::string Pass = passwd;
      p_connection = connectionFactory->createConnection(User, Pass);
      p_connection->start();

      // Create a Session
      p_session = p_connection->createSession(Session::AUTO_ACKNOWLEDGE);

      // Create the destination (Topic or Queue)
      p_destination = p_session->createTopic(topic);

      // Create a MessageProducer from the Session to the Topic
      p_producer = p_session->createProducer(p_destination);
      p_producer->setDeliveryMode(DeliveryMode::NON_PERSISTENT);

      p_channel = true;

      gridpack::utility::CoarseTimer *timer =
        gridpack::utility::CoarseTimer::instance();
      char sbuf[128];
      sprintf(sbuf,"Simulation started %f\n",timer->currentTime());
      std::auto_ptr<TextMessage>
        message(p_session->createTextMessage(sbuf));
      p_producer->send(message.get());
    } else {
      if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
        printf("ERROR: Channel already opened\n");
      }
    }
  }

  /**
   * Close IO channel
   */
  void closeChannel()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_connection) delete p_connection;
      if (p_session) delete p_session;
      if (p_destination) delete p_destination;
      if (p_producer) delete p_producer;
      p_channel = false;
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
    if (p_useFile) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
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
      printf(buf);
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
    int **index;
    int *indexbuf;
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      index = new int*[nwrites];
      indexbuf = new int[nwrites];
      iptr = indexbuf;
      int ones[nwrites];
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
        NGA_Scatter(p_stringGA,strbuf,index,nwrites);
        NGA_Scatter(p_maskGA,ones,index,nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
      delete [] index;
      delete [] indexbuf;
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
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          indexbuf = new int[nwrites];
          iptr = indexbuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
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
          delete [] index;
          delete [] indexbuf;
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
      std::auto_ptr<TextMessage> message(p_session->createTextMessage(p_channel_buf));
      printf("Sending message of length %d\n",p_channel_buf.length());
      p_producer->send(message.get());
      p_channel_buf.clear();
    }
  }
#endif

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
    char string[p_size];
    int nwrites = 0;
    int i;
    int one = 1;
    GA_Zero(p_maskGA);

    // Count up total strings being written from this processor
    for (i=0; i<nBus; i++) {
      if (p_network->getActiveBus(i) &&
          p_network->getBus(i)->serialWrite(string,p_size,signal)) {
        nwrites++;
      }
    }

    // Set up buffers to scatter strings to global buffer
    int **index;
    int *indexbuf;
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      index = new int*[nwrites];
      indexbuf = new int[nwrites];
      iptr = indexbuf;
      int ones[nwrites];
      char *strbuf;
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
        NGA_Scatter(p_stringGA,strbuf,index,nwrites);
        NGA_Scatter(p_maskGA,ones,index,nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
      delete [] index;
      delete [] indexbuf;
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
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          indexbuf = new int[nwrites];
          iptr = indexbuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
#ifndef USE_GOSS
              out << ptr;
#else
              if (p_channel) {
                p_channel_buf.append(ptr);
              } else {
                out << ptr;
              }
#endif
              ptr += p_size;
              nwrites++;
            }
          }
          delete [] index;
          delete [] indexbuf;
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
  void header(std::ostream & out, const char *str)
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
  }

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    bool p_useFile;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
#ifdef USE_GOSS
    bool p_channel;
    Connection *p_connection;
    Session *p_session;
    Destination *p_destination;
    MessageProducer *p_producer;
    std::string p_channel_buf;
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
    p_useFile = false;
  }

  /**
   * Simple Destructor
   */
  ~SerialBranchIO(void)
  {
    NGA_Deregister_type(p_GA_type);
    GA_Destroy(p_stringGA);
    GA_Destroy(p_maskGA);
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Redirect output to a file instead of standard out
   * @param filename name of file that output goes to
   */
  void open(const char *filename)
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
      if (!p_useFile) p_fout.reset(new std::ofstream);
      p_fout->open(filename);
      p_useFile = true;
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
    p_useFile = true;
  }

  /**
   * Close file and redirect output to standard out 
   */
  void close()
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      if (p_useFile && p_fout->is_open()) p_fout->close();
    }
  }

  /**
   * Write output from branches to standard out
   * @param signal an optional character string used to control contents of
   *                output
   */
  void write(const char *signal = NULL)
  {
    if (p_useFile) {
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
  void header(const char *str)
  {
    if (p_useFile) {
      header(*p_fout, str);
    } else {
      header(std::cout, str);
    }
  }

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
          " with allocated size: data: %d allocated: %ld\n",
          sizeof(_data_type),p_size);
      printf(buf);
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
    int **index;
    int *indexbuf;
    int *iptr;
    char *ptr;

    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      index = new int*[nwrites];
      indexbuf = new int[nwrites];
      iptr = indexbuf;
      int ones[nwrites];
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
        NGA_Scatter(p_stringGA,strbuf,index,nwrites);
        NGA_Scatter(p_maskGA,ones,index,nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
      delete [] index;
      delete [] indexbuf;
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
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          indexbuf = new int[nwrites];
          iptr = indexbuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
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
          delete [] index;
          delete [] indexbuf;
        }
      }
    }
    GA_Pgroup_sync(p_GAgrp);
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
    char string[p_size];
    int nwrites = 0;
    int i;
    int one = 1;
    GA_Zero(p_maskGA);

    // Count up total strings being written from this processor
    for (i=0; i<nBranch; i++) {
      if (p_network->getActiveBranch(i) &&
          p_network->getBranch(i)->serialWrite(string,p_size,signal)) nwrites++;
    }

    // Set up buffers to scatter strings to global buffer
    int **index;
    int *indexbuf;
    int *iptr;
    char *ptr;
    GA_Zero(p_maskGA);
    if (nwrites > 0) {
      index = new int*[nwrites];
      indexbuf = new int[nwrites];
      iptr = indexbuf;
      int ones[nwrites];
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
        NGA_Scatter(p_stringGA,strbuf,index,nwrites);
        NGA_Scatter(p_maskGA,ones,index,nwrites);
      }
      if (nwrites*p_size > 0) delete [] strbuf;
      delete [] index;
      delete [] indexbuf;
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
        int imask[ld];
        NGA_Get(p_maskGA,&lo,&hi,imask,&one);
        int j;
        nwrites = 0;
        for (j=0; j<ld; j++) {
          if (imask[j] == 1) {
            nwrites++;
          }
        }
        // Create buffers to retrieve strings from process i
        if (nwrites > 0) {
          char iobuf[p_size*nwrites];
          index = new int*[nwrites];
          indexbuf = new int[nwrites];
          iptr = indexbuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              index[nwrites] = iptr;
              *(index[nwrites]) = j + lo;
              nwrites++;
              iptr++;
            }
          }
          NGA_Gather(p_stringGA,iobuf,index,nwrites);
          ptr = iobuf;
          nwrites = 0;
          for (j=0; j<ld; j++) {
            if (imask[j] == 1) {
              out << ptr;
              ptr += p_size;
              nwrites++;
            }
          }
          delete [] index;
          delete [] indexbuf;
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
  void header(std::ostream & out, const char *str) const
  {
    if (GA_Pgroup_nodeid(p_GAgrp) == 0) {
      out << str;
    }
  }

  private:
    int p_GA_type;
    boost::shared_ptr<_network> p_network;
    int p_stringGA;
    int p_maskGA;
    int p_size;
    bool p_useFile;
    boost::shared_ptr<std::ofstream> p_fout;
    int p_GAgrp;
};

}   // serial_io
}   // gridpack
#endif  // _serial_io_h_
