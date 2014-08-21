// -------------------------------------------------------------
// file: thingc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created July 30, 2014 by William A. Perkins
// Last Change: 2014-08-21 09:43:20 d3g096
// -------------------------------------------------------------

#include <iostream>

extern "C" void *create_thing(int i);
extern "C" void thing_message(void *p);
extern "C" void destroy_thing(void *p);

// -------------------------------------------------------------
//  class Thing
// -------------------------------------------------------------
class Thing {
protected:

  void *my_fthing0;
  void *my_fthing1;
  void *my_fthing2;

  void my_init(void)
  {
    my_fthing0 = create_thing(0);
    my_fthing1 = create_thing(1);
    my_fthing2 = create_thing(2);
  }

public:

  /// Default constructor.
  Thing(void) 
    : my_fthing1(NULL), my_fthing2(NULL)
  {
    my_init();
  }

  /// Destructor
  ~Thing(void)
  {
    destroy_thing(my_fthing0);
    destroy_thing(my_fthing1);
    destroy_thing(my_fthing2);
  }

  void message(void)
  {
    thing_message(my_fthing0);
    thing_message(my_fthing1);
    thing_message(my_fthing2);
  }
};

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  Thing thing;
  thing.message();
  return 0;
}

