/*
  Пример использования функций из модуля @qsort.h@
  Быстрая сортировка
*/

#include <stdlib.h>
#include "nl.h"

int compare(size_t i, size_t j, void* data)
{
  if (((int*)data)[i] < ((int*)data)[j])
    return -1;
  return 1;
}

void swap(size_t i, size_t j, void* data)
{
  int tmp;
  tmp = ((int*)data)[i];
  ((int*)data)[i] = ((int*)data)[j];
  ((int*)data)[j] = tmp;
} 

int main()
{
  size_t j, n;
  int* v;

  n = 10;

  v = malloc(n * sizeof(int));

  for (j = 0; j < n; j++) 
    v[j] = rand() / 1000;
  
  printf("\nИсходный массив:\n");
  for (j = 0; j < n; j++) 
    printf("  %d", v[j]);

  nl_qsort(v, n, compare, swap);
  
  printf("\nОтсортированный массив:\n");
  for (j = 0; j < n; j++) 
    printf("  %d", v[j]);

  printf("\n");

  free(v);

  return 0;
}
