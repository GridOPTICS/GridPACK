#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minilzo.h"
#define KB 1024
#define debug 1
static lzo_align_t __LZO_MMODEL wrkmem [ (LZO1X_1_MEM_COMPRESS + sizeof(lzo_align_t) - 1) / sizeof(lzo_align_t) ];

void compress_and_write_(char *filename, float *array, int *array_size, int *count, int *offsets) {
  int size = *array_size;
  int out_size = size + size / 16 + 64 + 3;
  void* out_buf = malloc(out_size);
  FILE* f;
  int wrote;
  if(lzo_init() != LZO_E_OK) printf("Error initializing LZO.\n");
  if(debug) fprintf(stderr, "Compressing %d bytes...", size);
  if(lzo1x_1_compress((lzo_bytep) array,
                       size,
                       (lzo_bytep) out_buf,
                       (lzo_uintp) &out_size,
                       wrkmem) == LZO_E_OK) {
    if(debug) fprintf(stderr, "done\nCompressed size is %d bytes.\n", out_size);
  }
  else {
    printf("Failure.\nCompression did not succeed.\n");
  }
  f = fopen(filename, "w");
  fwrite(out_buf, 1, out_size, f);
  fwrite (offsets, sizeof(int), *count, f);
  fwrite(count, sizeof(int), 1, f);
  fclose(f);
  return;
}
