#include "nl.h"

/*
*/

int main()
{
  long *p, *d;
  long j, n, nd;
  long count;

  n = 25;

  p = (long*)malloc(n * sizeof(long));

    int_primes(n, p);

  printf("\nПростые числа:\n");
  for (j = 0; j < n; j++) 
    printf("%d ", p[j]);

  printf("\n");

  free(p);

  n = 12321331;
  if (int_is_prime(n))
      printf("%d --- простое\n", n);
  else 
      printf("%d --- составное\n", n);

  n = 123341;
  if (int_is_prime(n))
      printf("%d --- простое\n", n);
  else 
      printf("%d --- составное\n", n);

  n = 12312317;
  d = (long*)malloc(32 * sizeof(long));

  int_factor(n, d, &nd);

  printf("Разложение числа %d:\n", n);
  for (j = 0; j < nd; j++) 
    printf("%d ", d[j]);
  printf("\n");

  free(d);

    printf("gcd(%d, %d) = %d\n", 221, 323, int_gcd(221, 323));
    printf("gcd(%d, %d) = %d\n", 0, 0, int_gcd(0, 0));
    printf("gcd(%d, %d) = %d\n", 0, 3, int_gcd(0, 3));
    printf("gcd(%d, %d) = %d\n", 3, 0, int_gcd(3, 0));

  return 0;
}
