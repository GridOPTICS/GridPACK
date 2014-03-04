void c_overwrite(int *esize, void **a_ptr, void **b_ptr)
{
  memcpy(*b_ptr, *a_ptr, *esize);
}
