#include <stdio.h>
#include "bn.h"
char print_buf[8192];

int gcdEx(int a, int b, int *x, int *y)
{
    if (b == 0)
    {
        *x = 1, *y = 0;
        return a;
    }
    else
    {
        int r = gcdEx(b, a % b, x, y);
        int t = *x;
        *x = *y;
        *y = t - a / b * *y;
        return r;
    }
}

bignum gcdEx2(bignum a, bignum b, bignum *x, bignum *y)
{
    if (bignum_is_zero(&b))
    {
        bignum_from_int(x, 1);
        bignum_from_int(y, 0);
        return a;
    }
    else
    {
        bignum a_mod_b;
        bignum_init(&a_mod_b);
        bignum_mod(&a, &b, &a_mod_b);
        bignum r = gcdEx2(b, a_mod_b, x, y);
        bignum t;
        bignum_assign(&t, x);
        bignum_assign(x, y);

        // *y = t - a/b * (*y)
        bignum a_div_b;
        bignum_div(&a, &b, &a_div_b);
        bignum a_div_b_mul_y;
        bignum_mul(&a_div_b, y, &a_div_b_mul_y);
        bignum tmp;
        bignum_sub(&t, &a_div_b_mul_y, &tmp);
        bignum_assign(y, &tmp);
        return r;
    }
}

int main()
{
    // int a = 3, b = 20;
    // int x = 0, y = 0;
    // int res = gcdEx(a, b, &x, &y);
    // printf("the mod rev : %d\n", x);

    bignum a, b, x, y;
    bignum_from_int(&a, 9);
    bignum_from_int(&b, 20);
    bignum_init(&x);
    bignum_init(&y);
    gcdEx2(a, b, &x, &y);
    bignum_to_string(&x, print_buf, 8192);
    printf("the x is : 0x%s\n", print_buf);
    printf("hello!\n");
}