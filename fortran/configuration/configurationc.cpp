/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// file: configurationc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 30, 2014 by William A. Perkins
// Last Change: 2014-08-22 08:36:22 d3g096
// -------------------------------------------------------------

#include "cursor_wrapper.hpp"

extern "C" CursorWrapper *
configuration_open(gridpack::parallel::Communicator *comm, 
                   const char *conf_name)
{
  // std::cout << conf_name << std::endl;
  CursorWrapper *result(new CursorWrapper);
  result->cursor.reset();

  gridpack::utility::Configuration 
    *config(gridpack::utility::Configuration::configuration());
  
  // config->enableLogging(&std::cout);
  if (config->open(conf_name, *comm)) {
    result->cursor = config->getCursor("");
  } 
  return result;
}

extern "C" CursorWrapper *
configuration_cursor(CursorWrapper *cwrap, 
                     const char *path)
{
  CursorWrapper *result(new CursorWrapper);
  result->cursor = cwrap->cursor->getCursor(path);
  return result;
}

extern "C" bool
configuration_get_bool(CursorWrapper *cwrap, 
                       char *key,
                       bool *flag)
{
  bool result;
  result = cwrap->cursor->get(key, flag);
  // std::cout << "get_bool: " << key << ": " << *flag << std::endl;
  return result;
}

extern "C" bool
configuration_get_int(CursorWrapper *cwrap, 
                       char *key,
                       int *i)
{
  bool result;
  result = cwrap->cursor->get(key, i);
  // std::cout << "get_int: " << key << ": " << *i << std::endl;
  return result;
}

extern "C" bool
configuration_get_double(CursorWrapper *cwrap, 
                       char *key,
                       double *d)
{
  bool result;
  result = cwrap->cursor->get(key, d);
  // std::cout << "get_double: " << key << ": " << *d << std::endl;
  return result;
}

extern "C" bool
configuration_get_string(CursorWrapper *cwrap, 
                       char *key,
                       char *s,
                       int *slen)
{
  bool result;
  std::string string = s;
  result = cwrap->cursor->get(key, &string);
  strcpy(s,string.c_str());
  *slen = strlen(s);
  // std::cout << "get_string: " << key << ": " << *d << std::endl;
  return result;
}

extern "C" bool
configuration_ok(CursorWrapper *cwrap)
{
  bool result(false);
  if (cwrap != NULL) {
    result = cwrap->cursor;
  }
  return result;
}


extern "C" void
configuration_destroy(CursorWrapper **cwrap)
{
  delete *cwrap;
  *cwrap = NULL;
  return;
}
