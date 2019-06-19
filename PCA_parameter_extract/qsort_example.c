#include <stdlib.h>
#include <stdio.h>


int cmp (double *num1, double *num2)
{
  if (*num1 < *num2) return -1;
  else if (*num1 == *num2) return 0;
  return 1;
}


int main (void)
{
  double array[] = { 6.0, 12.0, 7.0, 22.0, 47.0, 16.0, -2.0, 60.0, 55.0 };
  int arraylen = sizeof(array) / sizeof(double);
  int i;

  qsort (array,arraylen,sizeof(double),
                    (int (*)(const void *, const void *)) cmp);

  // now display the sorted array

  for (i=0; i < arraylen; i++)
    printf("%lf ",array[i]);

  return 0;
}
