#ifndef __TASCEL_COUNTERS_H__
#define __TASCEL_COUNTERS_H__

#include <mpi.h>

#if defined(__cplusplus)
extern "C" {
#endif

  void counters_create_();
  
  void counters_destroy_();

  void counters_process_(MPI_Fint *count);


#if defined(__cplusplus)
}
#endif


#endif /* __TASCEL_COUNTERS_H__ */
