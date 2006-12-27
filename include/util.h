#ifndef __UTIL_H__
#define __UTIL_H__ 

#include <malloc.h>
#include <math.h>
#include <stdio.h>

/**
  \file 
  Вспомогательные процедуры. 

  Файл содержит вспомогательные процедуры для библиотеки 
  численных методов.
*/

/*
  Тип данных для представления индексов векторов и матриц.
  Оставлено для совместимости с прежними версиями.
  Сейчас напрямую используем size_t
*/
typedef size_t size_t;

/**
  \name Обработка ошибок
*/

/**
  Коды ошибок
*/

typedef enum
{
  nl_err_OK                    = 0,	  ///< Ошибок нет
  nl_err_out_of_memory         = 1,	  ///< Недостаточно памяти
  nl_err_inconsistent_size,		  ///< Неверный размер
  nl_err_incorrect_index,		  ///< Неверный индекс
  nl_err_matrix_must_be_square,		  ///< Матрица должна быть квадратной
  nl_err_matrix_is_singular,		  ///< Матрица вырождена
  nl_err_matrix_is_not_positive_definite, ///< Матрица не положительно определенная
  nl_err_IO,				  ///< Ошибка при вводе/выводе
  nl_err_diag_elem_must_be_non_zero,	  ///< Диагональные элементы должны быть ненулевыми
  nl_err_count,				  ///< Количество типов зарегистрированных ошибок
} nl_err_enum;

/**
  Указатель на функцию, которая вызывается в случае возникновения ошибки.
  По умолчанию печатает на консоли сообщение об ошибке и в случае, 
  если #critical отлично от нуля или переменная 
  #nl_exit_when_error отлична от нуля, аварийно завершает работу программы.
  В случае необходимости пользователь может написать свой обработчик ошибок.
*/
extern void (*nl_error)(nl_err_enum error, int critical);

/**
  Если значение этой переменной отлично от нуля, то в случае возникновения ошибки
  процедура #nl_error прерывает работу программы.
*/
extern int nl_exit_when_error;



/**                         
  \name Создание, удаление и копирование векторов
*/

/**
  Создать вектор длины \f$n\f$. 
  Функция размещает память под вектор и возвращает указатель на него.
  Тип элементов - double
*/
extern double* nl_dvector_create(size_t n);

/**
  Создать вектор длины \f$n\f$. 
  Функция размещает память под вектор и возвращает указатель на него.
  Тип элементов - #size_t
*/
extern size_t* nl_xvector_create(size_t n);

/**
  Создать копию вектора \f$v\f$ длины \f$n\f$. 
  Функция размещает память под вектор, копирует значения и возвращает
  указатель на него.
  Тип элементов - double
*/
extern double* nl_dvector_create_copy(const double* v, size_t n);

/**
  Создать копию вектора \f$v\f$ длины \f$n\f$. 
  Функция размещает память под вектор, копирует значения и возвращает
  указатель на него.
  Тип элементов - #size_t
*/
extern size_t* nl_xvector_create_copy(size_t* v, size_t n);

/**
  Копирует вектор \f$b\f$ в вектор \f$a\f$.
  Тип элементов - double
*/
extern void nl_dvector_copy(double* a, const double* b, size_t n);

/**
  Копирует вектор \f$b\f$ в вектор \f$a\f$.
  Тип элементов - #size_t
*/
extern void nl_xvector_copy(size_t* a, const size_t* b, size_t n);

/**
  Освободить память, занятую вектором \f$v\f$ длины \f$n\f$. 
  Тип элементов - double
*/
extern void nl_dvector_free(double* v);

/**
  Освободить память, занятую вектором \f$v\f$ длины \f$n\f$. 
  Тип элементов - #size_t
*/
extern void nl_xvector_free(size_t* v);


/**
  \name Ввод и вывод векторов
*/

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ на экран.
  Тип элементов - double
*/
extern void nl_dvector_print(double* v, size_t n, const char* format);

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ на экран.
  Тип элементов - #size_t
*/
extern void nl_xvector_print(size_t* v, size_t n, const char* format);

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ в файл file.
  Тип элементов - double
*/
extern void nl_dvector_fprint(FILE* file, const double* v, size_t n, const char* format);

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ в файл file.
  Тип элементов - size_t
*/
extern void nl_xvector_fprint(FILE* file, const size_t* v, size_t n, const char* format);

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ в файл с именем \f$filename\f$.
  Тип элементов - double
*/
extern void nl_dvector_fwrite(const char* filename, double* v, size_t n, const char* format);

/**
  Форматированный вывод вектора \f$v\f$ длины \f$n\f$ в файл с именем \f$filename\f$.
  Тип элементов - #size_t
*/
extern void nl_xvector_fwrite(const char* filename, size_t* v, size_t n, const char* format);

/**
  Чтение вектора с клавиатуры.
  Функция считывает с клавиатуры размер вектора \f$n\f$, размещает память под
  вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - double
*/
extern double* nl_dvector_scan(size_t* n);

/**
  Чтение вектора с клавиатуры.
  Функция считывает с клавиатуры размер вектора \f$n\f$, размещает память под
  вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - #size_t
*/
extern size_t* nl_xvector_scan(size_t* n);

/**
  Чтение вектора из файла.
  Функция считывает из файла \f$file\f$ размер вектора \f$n\f$, размещает память под
  вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - double
*/
extern double* nl_dvector_fscan(FILE* file, size_t* n);

/**
  Чтение вектора из файла.
  Функция считывает из файла \f$file\f$ размер вектора \f$n\f$, размещает память под
  вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - #size_t
*/
extern size_t* nl_xvector_fscan(FILE* file, size_t* n);

/**
  Чтение вектора из файла.
  Функция считывает из файла с именем \f$filename\f$ размер вектора \f$n\f$, размещает
  память под вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - double 
*/
extern double* nl_dvector_fread(const char* filename, size_t* n);

/**
  Чтение вектора из файла.
  Функция считывает из файла с именем \f$filename\f$  размер вектора \f$n\f$, размещает
  память под вектор, считывает значения и возвращает указатель на вектор.
  Тип элементов - #size_t
*/
extern size_t* nl_xvector_fread(const char* filename, size_t* n);



/**
  \name Операции с векторами
*/

/**
	Сложение векторов.
	- Вход:
	  - \f$a\f$ - вектор
	  - \f$b\f$ - вектор
	  - \f$n\f$ - размер векторов
	- Выход:
	  - \f$a = a+b\f$
*/
extern void nl_dvector_add(double* a, const double* b, size_t n);


/**
	Вычитание из вектора \f$a\f$ вектора \f$b\f$.
	- Вход:
	  - \f$a\f$ - вектор
	  - \f$b\f$ - вектор
	  - \f$n\f$ - размер векторов
	- Выход:
	  - \f$a = a-b\f$
*/
extern void nl_dvector_sub(double* a, const double* b, size_t n);

/**
  Вычисление 1-нормы (манхеттенской) разности двух векторов.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$b\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер векторов \f$a\f$ и \f$b\f$
  - Выход:
  - Результат: \f$\|a-b\|_1\f$

*/
extern double nl_dvector_diff_norm1(const double* a, const double* b, size_t n);

/**
  Вычисление 2-нормы разности (евклидовой) двух векторов.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$b\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер векторов \f$a\f$ и \f$b\f$
  - Выход:
  - Результат: \f$\|a-b\|_2\f$

*/
extern double nl_dvector_diff_norm2(const double* a, const double* b, size_t n);

/**
  Вычисление 1-нормы (манхеттенской) вектора.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер вектора \f$a\f$
  - Выход:
  - Результат: \f$\|a\|_1\f$

*/
extern double nl_dvector_norm1(const double* a, size_t n);


/**
  Вычисление 2-нормы (евклидовой) вектора.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер вектора \f$a\f$
  - Выход:
  - Результат: \f$\|a\|_2\f$

*/
extern double nl_dvector_norm2(const double* a, size_t n);


/**
  Вычисление \f$\infty\f$-нормы (чебышевой) вектора.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер вектора \f$a\f$
  - Выход:
  - Результат: \f$\|a\|_{\infty}\f$

*/
extern double nl_dvector_norm_inf(const double* a, size_t n);

/**
  Вычисление скалярного произведения векторов.
  - Вход:
	- \f$a\f$ - вектор размера \f$n\f$
	- \f$b\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер вектора \f$a\f$
  - Выход:
  - Результат: \f$(a,b)\f$

*/
extern double nl_dvector_dot(const double* a, const double* b, size_t n);

/**
  Вычисление вектора \f$z = \alpha x + y\f$.
  - Вход:
	- \f$alpha\f$ - коэффициент перед \f$x\f$
	- \f$x\f$ - вектор размера \f$n\f$
	- \f$y\f$ - вектор размера \f$n\f$
	- \f$n\f$ - размер вектора \f$a\f$
  - Выход:
	- \f$z\f$ - вектор размера \f$n\f$
*/
extern void nl_dvector_axpy(double alpha, const double *x, const double *y, size_t n, double *z);

/**
  Перестановка элементов вектора.
  - Вход:
	- \f$a\f$ - вектор длины \f$n\f$
	- \f$P\f$ - перестановка чисел \f$1,2,\dots,n\f$
	- \f$n\f$ - длина векторов \f$a\f$, \f$P\f$
  - Выход:
	- \f$b\f$ - вектор, составленный из компонент вектора \f$a\f$, переставленных согласно \f$P\f$
*/

extern void nl_dvector_permute(const double *a, const size_t *P, size_t n, double *b);



/**
  \name Создание, удаление и копирование матриц
*/

/**
  Создать матрицу с размерами \f$m\times n\f$. 
  Функция размещает память под матрицу и возвращает указатель на него.
  Тип элементов - double
*/
extern double** nl_dmatrix_create(size_t m, size_t n);

/**
  Создать копию матрицы \f$A\f$ с размерами \f$m\times n\f$. 
  Функция размещает память под матрицу, копирует значения и возвращает
  указатель на него.
  Тип элементов - double
*/
extern double** nl_dmatrix_create_copy(const double** A, size_t m, size_t n);

/**
  Копирует \f$m\times n\f$-матрцу B матрицу \f$A\f$.
  Под элементы матрицы \f$A\f$ должно быть выделено место.
  Тип элементов - double
*/
extern void nl_dmatrix_copy(double** A, const double** B, size_t m, size_t n);

/**
  Освободить память, занятую матрицей \f$A\f$ с размерами \f$m\times n\f$. 
  Тип элементов - double
*/
extern void nl_dmatrix_free(double** A, size_t m);

/**
  Конвертация представления матрицы.
  Матрица \f$A\f$ с размерами \f$m\times n\f$, представленная в виде списка своих
  элементов (по строкам), конвертируется в представление,
  принятое в библиотеке. Для результирующе матрицы выделяется память.
  Тип элементов - double
*/
extern double** nl_dmatrix_convert(double* A, size_t m, size_t n);



/**
  \name Ввод и вывод матриц
*/


/**
  Форматированный вывод \f$m\times n\f$ матрицы \f$A\f$ на экран.
  Тип элементов - double
*/
extern void nl_dmatrix_print(double **A, size_t m, size_t n, const char* format);

/**
  Форматированный вывод \f$m\times n\f$ матрицы \f$A\f$ в файл \f$file\f$.
  Тип элементов - double
*/
extern void nl_dmatrix_fprint(FILE* file, double** A, size_t m, size_t n, const char* format);

/**
  Форматированный вывод \f$m\times n\f$ матрицы \f$A\f$ в файл с именем \f$filename\f$.
  Тип элементов - double
*/
extern void nl_dmatrix_fwrite(const char* filename, double** A, size_t m, size_t n, const char* format);

/**
  Чтение матрицы с клавиатуры.
  Функция считывает с клавиатуры размеры матрицы \f$m\f$, \f$n\f$, размещает память под
  матрицу, считывает значения и возвращает указатель на матрицу.
  Тип элементов - double
*/
extern double** nl_dmatrix_scan(size_t* m, size_t* n);

/**
  Чтение матрицы из файла.
  Функция считывает из файла \f$file\f$ размеры матрицы \f$m\f$, \f$n\f$, размещает память под
  матрицу, считывает значения и возвращает указатель на матрицу.
  Тип элементов - double
*/
extern double** nl_dmatrix_fscan(FILE* file, size_t* m, size_t* n);

/**
  Чтение матрицы из файла.
  Функция считывает из файла с именем \f$filename\f$ размеры матрицы \f$m\f$, \f$n\f$,
  размещает память под матрицу, считывает значения и возвращает указатель
  на матрицу.
  Тип элементов - double 
*/
extern double** nl_dmatrix_fread(const char* filename, size_t* m, size_t* n);


/**
  \name Операции с матрицами 
*/

/**
	Умножает матрицу \f$A\f$ размера \f$m\times n\f$ на матрицу B \f$n \times k\f$
	Результат в \f$Mul\f$ размера \f$m\times k\f$. 
*/
extern void nl_dmatrix_mult(size_t m, size_t n, size_t k, double** A, double** B, double** Mul);

/**
	Сложение двух матриц.
	- Вход:
	  - \f$m\f$ - число строк матриц
	  - \f$n\f$ - число столбцов матриц
	  - \f$A\f$ - матрица
	  - \f$B\f$ - матрица
	- Выход:
	  - \f$A = A+B\f$
*/
extern void nl_dmatrix_add(size_t m, size_t n, double** A, double** B);

/**
	Вычитает из матрицы \f$A\f$ размера \f$m\times n\f$ матрицу B размера \f$m\times n\f$
	Результат в матрице \f$А\f$.
*/
extern void nl_dmatrix_sub(size_t m, size_t n, double** A, double** B);

/**
  Вычисление чебышевой нормы разности двух матриц.
  - Вход:
    - \f$A\f$ - матрица размера \f$m\times n\f$
    - \f$B\f$ - матрица размера \f$m\times n\f$
    - \f$m\f$, \f$n\f$ - размеры матриц \f$A\f$ и \f$B\f$
  - Выход:
  - Результат: \f$\|A-B\|_{\infty}\f$

*/
extern double nl_dmatrix_diff_norm_inf(double** A, double** B, size_t m, size_t n);


/**
  \name Операции с матрицами и векторами
*/


/**
  Умножение вектора-строки \f$b\f$ на матрицу \f$A\f$.
  Функция выделяет память под результирующий вектор и записывает в него результат
  умножения.
  - Вход:
       - \f$b\f$ - вектор-строка
       - \f$A\f$ - исходная матрица
       - \f$m\f$ - число строк матрицы \f$A\f$
       - \f$n\f$ - число столбцов матрицы \f$A\f$
  - Результат:
       - \f$A\cdot b\f$
*/
extern double* nl_dvector_mult_dmatrix(double* b, double** A, size_t m, size_t n);


/**
  Умножение матрицы \f$A\f$ на вектор \f$b\f$ справа.

  - Вход:
       - \f$A\f$ - исходная матрица
       - \f$m\f$ - число строк матрицы \f$A\f$
       - \f$n\f$ - число столбцов матрицы \f$A\f$
       - \f$b\f$ - вектор-столбец
  - Выход:
       - \f$res = A\cdot b\f$
*/
extern void nl_dmatrix_mult_dvector(double** A, size_t m, size_t n, const double* b, double* res);

extern const char double_format[];
extern const char integer_format[];


#endif
