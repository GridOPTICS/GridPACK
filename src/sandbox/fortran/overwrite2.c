void c_overwrite_(int *esize, void *a_ptr, void *b_ptr)
{
  int i;
  char *a, *b;
  a = (char*)a_ptr;
  b = (char*)b_ptr;
  printf("esize: %d\n",*esize);
  for (i=0; i<*esize; i++) {
    b[i] = a[i];
    printf("a[%d]: %c b[%d]: %c\n",i,a[i],i,b[i]);
  }

  /*memcpy(b_ptr, a_ptr, esize); */
}
