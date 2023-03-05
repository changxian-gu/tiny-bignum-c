#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memcpy */
#include "bn.h"

/* O(log n) */
void pow_mod_faster(struct bn* a, struct bn* b, struct bn* n, struct bn* res)
{
  bignum_from_int(res, 1); /* r = 1 */

  struct bn tmpa;
  struct bn tmpb;
  struct bn tmp;
  bignum_assign(&tmpa, a);
  bignum_assign(&tmpb, b);

  while (1)
  {
    if (tmpb.array[0] & 1)     /* if (b % 2) */
    {
      bignum_mul(res, &tmpa, &tmp);  /*   r = r * a % m */
      bignum_mod(&tmp, n, res);
    }
    bignum_rshift(&tmpb, &tmp, 1); /* b /= 2 */
    bignum_assign(&tmpb, &tmp);

    if (bignum_is_zero(&tmpb))
      break;

    bignum_mul(&tmpa, &tmpa, &tmp);
    bignum_mod(&tmp, n, &tmpa);
  }
}


// 如果num是偶数直接排除
int miller_rabin(bignum *num)
{
  char print_buf[8192];
  bignum one, two, three;
  bignum_from_int(&one, 1);
  bignum_from_int(&two, 2);
  bignum_from_int(&three, 3);
  if (bignum_cmp(num, &two) < 0)
    return 0;
  if (bignum_cmp(num, &two) == EQUAL || bignum_cmp(num, &three) == EQUAL)
    return 1;
  bignum tmp;
  bignum_init(&tmp);
  bignum_mod(num, &two, &tmp);
  if (bignum_is_zero(&tmp))
    return 0;

  bignum temp, s;
  bignum_init(&tmp);
  bignum_init(&s);
  bignum_init(&temp);
  bignum_dec(num);
  bignum_assign(&s, num);
  bignum_inc(num);

  uint64_t t = 0;

  bignum_from_int(&temp, 2);
  bignum_mod(&s, &temp, &tmp);

  while (bignum_is_zero(&tmp))
  {
    bignum_init(&temp);
    bignum_rshift(&s, &temp, 1);
    bignum_assign(&s, &temp);
    bignum_assign(&tmp, &s);
    t++;
  }

  int trials = 5;
  char buf[256];
  FILE *fp = fopen("/dev/urandom", "r");
  for (int i = 0; i < trials; i++)
  {
    bignum a;
    bignum tmp, tmp2;
    fread(buf, 256, 1, fp);
    bignum_from_string(&a, buf, 256);
    // # a 是2 和 num - 1 之间的数字
    // # a = (rand() % (num - 1 - 2)) + 2
    bignum_from_int(&tmp, 3);
    bignum_sub(num, &tmp, &tmp2);
    bignum_init(&tmp);
    bignum_mod(&a, &tmp2, &tmp);
    bignum_inc(&tmp);
    bignum_inc(&tmp);
    bignum_assign(&a, &tmp);

    // # v = a ^ s % num
    bignum v;
    bignum_init(&v);
    pow_mod_faster(&a, &s, num, &v);

    // if v != 1
    if (bignum_cmp(&v, &one) != EQUAL)
    {
      uint64_t i = 0;
      bignum num_sub_1;
      bignum_sub(num, &one, &num_sub_1);
      // while (v != num - 1):
      while (bignum_cmp(&v, &num_sub_1) != EQUAL)
      {
        if (i == t - 1)
          return 0;
        else
        {
          i++;
          // v = (v ^ 2) % num
          pow_mod_faster(&v, &two, num, &tmp);
          bignum_assign(&v, &tmp);
        }
      }
    }
  }
  fclose(fp);
  return 1;
}

int is_prime(bignum *num)
{
  bignum t, t2;
  bignum_init(&t);
  bignum_init(&t2);
  bignum_from_int(&t, 2);
  if (bignum_cmp(num, &t) < 0)
    return 0;
  // 168 个 素数
  int low_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                      103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
                      229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
                      367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
                      503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
                      653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
                      821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971,
                      977, 983, 991, 997};

  for (int i = 0; i < sizeof(low_primes) / sizeof(int); i++)
  {
    bignum_from_int(&t, low_primes[i]);
    if (bignum_cmp(num, &t) == EQUAL)
      return 1;
  }

  for (int i = 0; i < sizeof(low_primes) / sizeof(int); i++)
  {
    bignum_from_int(&t, low_primes[i]);
    bignum_mod(num, &t, &t2);
    if (bignum_is_zero(&t2))
      return 0;
  }

  return miller_rabin(num);
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

void modInverse(bignum *e, bignum *n, bignum* d)
{
  bignum x, y;
  bignum_init(&x);
  bignum_init(&y);
  gcdEx2(*e, *n, &x, &y);
  bignum_assign(d, &x);
}

void rsa_encrypt(bignum* msg, bignum* e, bignum* n, bignum* cipher) {
  pow_mod_faster(msg, e, n, cipher);
}

void rsa_decrypt(bignum* cipher, bignum* d, bignum* n, bignum* msg) {
  pow_mod_faster(cipher, d, n, msg);
}


static void test_rsa1024(void)
{
  char pool[8192];
  char print_buf[8192];
  bignum p, q;
  char pstr[] = "232be540505425be64d7e6d5e0acbd71e1dc0abe0ab447a14b997b8946da44ec56031d233a2e9747bebbca6145b58b3fbb6e46a1ac9c3602b39198f53e6e0489";
  char qstr[] = "df5f34028c341f0a73c40dcc61efd9391049e30c487b0e37f13a163f712e4f9bb02e550eaad9b34ee3002feaa4dbedc5185444c77fb12d605001c43987169b49";
  bignum_from_string(&p, pstr, 128);
  bignum_from_string(&q, qstr, 128);
  /*
    生成大素数
  */
  // bignum_init(&p);
  // bignum_init(&q);
  // FILE *fp = fopen("/dev/urandom", "r");
  // fread(pool, 8192, 1, fp);
  // int index = 0;
  // while (1) {
  //   memcpy(p.array, pool + index, 64);
  //   bignum_to_string(&p, print_buf, 8192);
  //   // printf("the p is : %s\n", print_buf);
  //   index += 64;
  //   if (is_prime(&p))
  //     break;
  //   if (index == 8192) {
  //     fread(pool, 8192, 1, fp);
  //     index = 0;
  //   }
  // }
  bignum_to_string(&p, print_buf, 8192);
  printf("the p : 0x%s\n", print_buf);

  // fread(pool, 8192, 1, fp);
  // while (1) {
  //   memcpy(q.array, pool + index, 64);
  //   bignum_to_string(&q, print_buf, 8192);
  //   // printf("the q is : %s\n", print_buf);
  //   index += 64;
  //   if (is_prime(&q))
  //     break;
  //   if (index == 8192) {
  //     fread(pool, 8192, 1, fp);
  //     index = 0;
  //   }
  // }
  bignum_to_string(&q, print_buf, 8192);
  printf("the q : 0x%s\n", print_buf);
  // fclose(fp);

  // n = pq
  bignum n;
  bignum_init(&n);
  bignum_mul(&p, &q, &n);
  bignum_to_string(&n, print_buf, 8192);
  printf("the n is : 0x%s\n", print_buf);

  // fi = (p - 1) (q - 1)
  bignum fi;
  bignum_init(&fi);
  bignum_dec(&p);
  bignum_dec(&q);
  bignum_mul(&p, &q, &fi);
  bignum_to_string(&fi, print_buf, 8192);
  printf("the fi is : 0x%s\n", print_buf);

  bignum e;
  bignum_from_int(&e, 65537);
  bignum_to_string(&e, print_buf, 8192);
  printf("the e is : 0x%s\n", print_buf);

  bignum d;
  char dstr[] = "166bc90b2b1b6374f465509f5cddbe1f1e29e70a6b5dda2766d048cd8caa29dfc7c63ad781364354da2771c20641ef9318aeac13fee2eca2aa28b0f7846ff97f0a34da3c797cbdb860d631e48cfb23afbff0a2c998740392ec31d33e1a83259c9ec8f3a855e33abdf36cc0e5c19cf609b479e59fc77e3701e64660870ab45041";
  // bignum_init(&d);
  // modInverse(&e, &fi, &d);
  bignum_from_string(&d, dstr, 256);
  bignum_to_string(&d, print_buf, 8192);
  printf("the d is : 0x%s\n", print_buf);

  bignum msg;
  bignum_from_int(&msg, 0x5722);
  bignum_to_string(&msg, print_buf, 8192);
  printf("the msg is : 0x%s\n", print_buf);

  bignum cipher;
  rsa_encrypt(&msg, &e, &n, &cipher);
  bignum_to_string(&cipher, print_buf, 8192);
  printf("the cipher is : 0x%s\n", print_buf);

  rsa_decrypt(&cipher, &d, &n, &msg);
  bignum_to_string(&msg, print_buf, 8192);
  printf("the msg is : 0x%s\n", print_buf);
}

int main()
{
  test_rsa1024();
  return 0;
}
