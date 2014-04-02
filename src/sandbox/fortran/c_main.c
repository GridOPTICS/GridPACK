int main (int argc, char **argv)
{
  int i;
  double f;
  int max = 10;
  for (i=0; i<max; i++) {
    for_print(i,&f);
    printf("From C: f = %6.1f\n",f);
  }
}
