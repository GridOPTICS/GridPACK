// -------------------------------------------------------------
// file: configurationc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 30, 2014 by William A. Perkins
// Last Change: 2014-08-14 12:32:31 d3g096
// -------------------------------------------------------------

#include <gridpack/configuration/configuration.hpp>

struct cursorWrapper {
  gridpack::utility::Configuration::CursorPtr cursor;
};

extern "C" cursorWrapper *
configuration_open(gridpack::parallel::Communicator *comm, 
                   const char *conf_name)
{
  // std::cout << conf_name << std::endl;
  cursorWrapper *result(new cursorWrapper);
  result->cursor.reset();

  gridpack::utility::Configuration 
    *config(gridpack::utility::Configuration::configuration());
  
  // config->enableLogging(&std::cout);
  if (config->open(conf_name, *comm)) {
    result->cursor = config->getCursor("");
  } 
  return result;
}

extern "C" cursorWrapper *
configuration_cursor(cursorWrapper *cwrap, 
                     const char *path)
{
  cursorWrapper *result(new cursorWrapper);
  result->cursor = cwrap->cursor->getCursor(path);
  return result;
}

extern "C" bool
configuration_get_bool(cursorWrapper *cwrap, 
                       char *key,
                       bool *flag)
{
  bool result;
  result = cwrap->cursor->get(key, flag);
  // std::cout << "get_bool: " << key << ": " << *flag << std::endl;
  return result;
}

extern "C" bool
configuration_get_int(cursorWrapper *cwrap, 
                       char *key,
                       int *i)
{
  bool result;
  result = cwrap->cursor->get(key, i);
  // std::cout << "get_int: " << key << ": " << *i << std::endl;
  return result;
}

extern "C" bool
configuration_get_double(cursorWrapper *cwrap, 
                       char *key,
                       double *d)
{
  bool result;
  result = cwrap->cursor->get(key, d);
  // std::cout << "get_double: " << key << ": " << *d << std::endl;
  return result;
}

extern "C" bool
configuration_ok(cursorWrapper *cwrap)
{
  bool result(false);
  if (cwrap != NULL) {
    result = cwrap->cursor;
  }
  return result;
}


extern "C" void
configuration_destroy(cursorWrapper **cwrap)
{
  delete *cwrap;
  *cwrap = NULL;
  return;
}
