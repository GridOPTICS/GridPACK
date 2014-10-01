/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   serial_io_c.cpp
 * @author Bruce Palmer
 * @date   2014-09-17 11:05:08 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/serial_io/serial_io.hpp"
#include "../component/fortran_component.hpp"

typedef gridpack::network::BaseNetwork<
  gridpack::fortran_component::FortranBusComponent,
  gridpack::fortran_component::FortranBranchComponent>
  FortranNetwork;

typedef gridpack::serial_io::SerialBusIO<FortranNetwork> FortranBusSerialIO;
typedef gridpack::serial_io::SerialBranchIO<FortranNetwork> FortranBranchSerialIO;

struct networkWrapper {
  boost::shared_ptr<FortranNetwork> network;
};

/**
 * Create a bus serial IO object
 * @param writer pointer to Fortran bus serial IO object
 * @param max_len the maximum string length written by any bus
 * @param network pointer to Fortran network object
 */
extern "C" void bus_serial_io_create(FortranBusSerialIO **writer,
    int max_len, networkWrapper *wnetwork)
{
  *writer = new FortranBusSerialIO(max_len, wnetwork->network);
}

/**
 * Destroy a bus serial IO object
 * @param writer pointer to Fortran bus serial IO object
 */
extern "C" void bus_serial_io_destroy(FortranBusSerialIO **writer)
{
  delete (*writer);
}

/**
 * Open an external file and redirect output to it
 * @param writer pointer to Fortran bus serial IO object
 * @param filename name of new file
 */
extern "C" void bus_serial_io_open(FortranBusSerialIO *writer, char *filename)
{
  writer->open(filename);
}

/**
 * Close an external file and redirect output to standard out
 * @param writer pointer to Fortran bus serial IO object
 */
extern "C" void bus_serial_io_close(FortranBusSerialIO *writer)
{
  writer->close();
}

/**
 * Write output from network to output
 * @param writer pointer to Fortran bus serial IO object
 */
extern "C" void bus_serial_io_write(FortranBusSerialIO *writer, char *signal)
{
  char *ptr;
  // make sure to pass in a null string if nothing is in signal
  if (strlen(signal) == 0) {
    ptr = signal;
  } else {
    ptr = NULL;
  }
  writer->write(ptr);
}

/**
 * Write single string to standard output. This is used to write headers for a
 * data listing. It is mostly a convenience function so that users do not have
 * to identify the head node
 * @param writer pointer to Fortran bus serial IO object
 * @param str character string containing the header
 */
extern "C" void bus_serial_io_header(FortranBusSerialIO *writer, char *str)
{
  int slen = strlen(str);
  str[slen] = '\n';
  str[slen+1] = '\0';
  writer->header(str);
}

/**
 * Create a branch serial IO object
 * @param writer pointer to Fortran branch serial IO object
 * @param max_len the maximum string length written by any branch
 * @param network pointer to Fortran network object
 */
extern "C" void branch_serial_io_create(FortranBranchSerialIO **writer,
    int max_len, networkWrapper *wnetwork)
{
  *writer = new FortranBranchSerialIO(max_len, wnetwork->network);
}

/**
 * Destroy a branch serial IO object
 * @param writer pointer to Fortran branch serial IO object
 */
extern "C" void branch_serial_io_destroy(FortranBranchSerialIO **writer)
{
  delete (*writer);
}

/**
 * Open an external file and redirect output to it
 * @param writer pointer to Fortran branch serial IO object
 * @param filename name of new file
 */
extern "C" void branch_serial_io_open(FortranBranchSerialIO *writer, char *filename)
{
  writer->open(filename);
}

/**
 * Close an external file and redirect output to standard out
 * @param writer pointer to Fortran branch serial IO object
 */
extern "C" void branch_serial_io_close(FortranBranchSerialIO *writer)
{
  writer->close();
}

/**
 * Write output from network to output
 * @param writer pointer to Fortran branch serial IO object
 */
extern "C" void branch_serial_io_write(FortranBranchSerialIO *writer, char *signal)
{
  char *ptr;
  // make sure to pass in a null string if nothing is in signal
  if (strlen(signal) == 0) {
    ptr = signal;
  } else {
    ptr = NULL;
  }
  writer->write(ptr);
}

/**
 * Write single string to standard output. This is used to write headers for a
 * data listing. It is mostly a convenience function so that users do not have
 * to identify the head node
 * @param writer pointer to Fortran branch serial IO object
 * @param str character string containing the header
 */
extern "C" void branch_serial_io_header(FortranBranchSerialIO *writer, char *str)
{
  int slen = strlen(str);
  str[slen] = '\n';
  str[slen+1] = '\0';
  writer->header(str);
}
