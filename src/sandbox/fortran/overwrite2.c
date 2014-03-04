#include <string.h>
void c_overwrite_(int *esize, void *a_ptr, void *b_ptr)
{
  memcpy(b_ptr, a_ptr, *esize);
}
