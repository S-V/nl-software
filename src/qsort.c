/* 
   Original version is written by Douglas C. Schmidt.
   Some changes adopting qsort to specific problems was done by Anton Tagunov.
   This new version by Nikolai Yu. Zolotykh is more universal but works slower... */

/* Copyright (C) 1991, 1992, 1996, 1997 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Written by Douglas C. Schmidt (schmidt@ics.uci.edu).

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with the GNU C Library; see the file COPYING.LIB.  If not,
   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/*
   Для чего нужна новая функция qsort?

   Если честно, то я (NZ) не совсем уверен, что она действительно нужна...
   
   Вот аргумент "за": новая функция более универсальна и удобна. Чтобы ей
   воспользоваться, пользователю нужно написать две функции: compare и
   swap. Функция compare сравнивает два элемента. Функция swap переставляет
   два элемента. Причем на вход этих функций передаются _индексы_
   элементов. Для того, чтобы воспользоваться стандартной Функцией qsort,
   пользователь должен написать только одну функцию compare. При этом на
   вход функции передаются указатели на эти элементы. Перестановку
   элементов стандартный qsort осуществляет сам, зная их размеры.
   Рассмотрим следующую задачу. Необходимо отсортировать два массива: a[0],
   a[1],..., a[n-1] и b[0], b[1],..., b[n-1] в порядке возрастания
   отношений: a[0]/b[0], a[1]/b[1],..., a[n-1]/b[n-1]. Чтобы
   воспользоваться стандартной функцией, нет другого способа, как
   переписать эти данные в один массив a[0], b[0], a[1], b[1],..., a[n-1],
   b[n-1]. При использовании нового варианта этого можно избежать.

   Вот аргумент "против": в новом варианте появляются накладные расходы:
     - из-за вызова функции swap; 
     - функциям swap и compare передаются не указатели на элементы массива, 
       а их индексы, что также приводит к замедлению программы;
   Кроме того, в новой версии qsort не удается на заключительном этапе
   использовать метод сортировки вставками, использование которого могло
   значительно ускорить реализацию. Использование swap для сортировки
   вставками дает сортировку пузырьком, что не оптимально. Как результат,
   получаем, что новый вариант qsort примерно в 10 раз медленнее
   стандартного. Заметим, что эти проблемы легко решаются в C++
*/


#include "util.h"
#include "qsort.h"

typedef struct 
{
  size_t lo;
  size_t hi;
} stack_node;

/* эмуляция стека */
#define QSORT_STACK_SIZE 32
#define QSORT_PUSH(LOW, HIGH) do {top->lo = LOW; top++->hi = HIGH;} while (0)
#define QSORT_POP(LOW, HIGH)  do {LOW = (--top)->lo; HIGH = top->hi;} while (0)
#define QSORT_STACK_NOT_EMPTY (stack < top)

#define QSORT_COMPARE(i, j) compare(i, j, data)
#define QSORT_SWAP(i, j) do {swap(i, j, data);} while(0)

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of MAX_INT is allocated on the
      stack.  Assuming a 32-bit integer, this needs only 32 *
      sizeof(stack_node) == 136 bits.  Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (n)
      stack size is needed (actually O(1) in this case)!  */

void nl_qsort(
  void* data,
  size_t n, 
  int (*compare)(size_t, size_t, void*),
  void (*swap)(size_t, size_t, void*))
{
  size_t lo = 0;
  size_t hi = n - 1;
  
  stack_node stack[QSORT_STACK_SIZE];
  stack_node *top = stack + 1;

  if (n <= 1)
    return;
  
  while (QSORT_STACK_NOT_EMPTY) 
  {
    size_t left, right;

    /* Select median value from among LO, MID, and HI. Rearrange
       LO and HI so the three values are sorted. This lowers the
       probability of picking a pathological pivot value and
       skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

    /* v[lo] is used as buffer and stored a pivot item -- NZ */

    size_t mid = (hi + lo) / 2;
    
    if (QSORT_COMPARE(lo, mid) < 0)
    {
      QSORT_SWAP(lo, mid);
    }
    if (QSORT_COMPARE(hi, lo) < 0)
    {
      QSORT_SWAP(hi, lo);

      if(QSORT_COMPARE(lo, mid) < 0)
      {
        QSORT_SWAP(lo, mid);
      }
    } 
    
    left  = lo + 1;
    right = hi - 1;
  
    /* Here's the famous ``collapse the walls'' section of quicksort.
       Gotta like those tight inner loops!  They are the main reason
       that this algorithm runs much faster than others. */

    do 
    {
      while(QSORT_COMPARE(left, lo) < 0)
        left++;
    
      while(QSORT_COMPARE(lo, right) < 0)
        right--;
    
      if (left < right) 
      {
        QSORT_SWAP(left, right);
        left++;
        right--;
      } 
      else
      {
        if (left == right) 
        {
          left++;
          right--;
          break;
        }
      }
    } while (left <= right);

    QSORT_SWAP(lo, right); /* recall the pivot item -- NZ */
    
    /* Set up pointers for next iteration.  First determine whether
       left and right partitions are below the threshold size.  If so,
       ignore one or both.  Otherwise, push the larger partition's
       bounds on the stack and continue sorting the smaller one. */

    if (right <= lo) 
    {
      if (hi <= left)
        QSORT_POP(lo, hi); 
      else
        lo = left;
    } 
    else
    { 
      if (hi <= left)
        hi = right;
      else if (right - lo > hi - left)
      {
        QSORT_PUSH(lo, right);
        lo = left;
      } 
      else 
      {
        QSORT_PUSH(left, hi);
        hi = right;
      }
    }
  }
}

#undef QSORT_STACK_SIZE
#undef QSORT_PUSH
#undef QSORT_POP
#undef QSORT_STACK_NOT_EMPTY
#undef QSORT_COMPARE
#undef QSORT_SWAP
