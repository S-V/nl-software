#include <math.h>
#include <stdio.h>

/*
  Тест на простоту и факторизация
*/

long int_gcd(long a, long b)
{
    long tmp;

    a = abs(a);
    b = abs(b);

    while(b != 0)
    {
        tmp = b;
        b = a % b;
        a = tmp;
    }

    return a;
}

void int_primes(long n, long* primes)  // n >= 2
{

    long p, j, k, pk, sqrtp;

    primes[0] = 2;
    primes[1] = 3;
    p = 3;

    for (j = 2; j < n; j++)
    {
        for(;;)
        {
            p += 2;
            sqrtp = sqrt(p);
            for(k = 1; ; k++)
            {
                pk = primes[k];
                if (p%pk == 0)
                    break;
                if (pk >= sqrtp)
                    goto prime_found;
            }
        }
prime_found:
        primes[j] = p;

    }
}

int int_is_prime(long n) // n != 0, 1
{

    long p, k, kmax;

    if (n == 2 || n == 3 || n == 5 || n == 7)
        return 1;

    if (n%2 == 0)
        return 0;
    if (n%3 == 0)
        return 0;

    kmax = (sqrt(n) + 1)/6;
    for (k = 1; k <= kmax; k++)
    {
        p = 6*k - 1;
        if (n%p == 0)
            return 0;
        p = 6*k + 1;
        if (n%p == 0)
            return 0;

    }

    return 1;
}

void int_factor(long n, long* divisors, long* nd)
{

    long p, k, kmax;

    *nd = 0;

    while (n%2 == 0)
    {
        divisors[(*nd)++] = 2;
        n /= 2;
    }
    while (n%3 == 0)
    {
        divisors[(*nd)++] = 3;
        n /= 3;
    }
        
    for (k = 1; ; k++)
    {
        p = 6*k - 1;
        if (n%p == 0)
        {
            divisors[(*nd)++] = p;
            n /= p;
        }
        if (p*p > n)
            break;
        p = 6*k + 1;
        if (n%p == 0)
        {
            divisors[(*nd)++] = p;
            n /= p;
        }
        if (p*p > n)
            break;
    }

    divisors[(*nd)++] = n;
}
