#include <math.h>
#include "util.h"

void fft_transform(double* re, double* im, size_t n)
{
    size_t i, j, k, kh;
    size_t len, halflen;
    double theta, pi;
    double w_re, w_im, u_re, u_im, t_re, t_im;
    double tmp, tmpr, tmpi;

    pi = 3.1415926535897932384;

    j = 0;
    for (i = 0; i < n - 1; i++)
    {
        if (i < j)
        {
            tmp = re[i]; re[i] = re[j]; re[j] = tmp;
            tmp = im[i]; im[i] = im[j]; im[j] = tmp;
        }
        k = n / 2;
        while (k <= j)
        {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    for (halflen = 1; halflen < n; halflen *= 2)
    {
        len = 2*halflen;
        theta  = -pi/halflen;
        w_re   = 1.0;
        w_im   = 0.0; 
        if (halflen == 1)
        {
            u_re = -1.0;
            u_im = 0.0;
        }
        else if (halflen == 2)
        {
            u_re = 0.0;
            u_im = -1.0;
        }  
        else 
        {
            u_re = cos(theta);
            u_im = sin(theta);
        }
        for (j = 0; j < halflen; j++)
        {
            for (k = j; k < n; k += len)
            {
                kh = k + halflen;

                t_re = re[kh]*w_re - im[kh]*w_im;
                t_im = re[kh]*w_im + im[kh]*w_re;

                re[kh] = re[k]  - t_re;
                im[kh] = im[k]  - t_im;

                re[k] += t_re;
                im[k] += t_im;
            }
            tmpr = u_re*w_re - u_im*w_im;
            tmpi = u_im*w_re + u_re*w_im;
            w_re = tmpr;
            w_im = tmpi;
        }
    }
}

void fft_inverse(double* re, double* im, size_t n)
{
    size_t j;

    for (j = 0 ; j < n; j++)
        im[j] = -im[j];

    fft_transform(re, im, n);

    for (j = 0 ; j < n; j++)
    {
        re[j] /= n;
        im[j] = -im[j]/n;
    }
}
