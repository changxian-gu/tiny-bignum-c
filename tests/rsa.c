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

void print_byte(char* buf, int len) {
  for (int i = 0; i < 1; i++) {
    printf("%02x", *((unsigned char*)buf + i));
  }
  printf("\n");
}

int miller_rabin(struct bn *num) {
  struct bn one, two;
  bignum_from_int(&one, 1);
  bignum_from_int(&two, 2);
  struct bn tmp, s;
  struct bn temp;
  bignum_init(&tmp);
  bignum_init(&s);
  bignum_init(&temp);
  bignum_mod(num, &two, &tmp);
  if (bignum_is_zero(&tmp))
    return 0;
  bignum_init(&tmp);

  bignum_dec(num);
  bignum_assign(&s, num);
  bignum_inc(num);

  uint64_t t = 0;
  bignum_from_int(&temp, 2);
  bignum_mod(&s, &temp, &tmp);
  while (bignum_is_zero(&tmp)) {
    bignum_init(&temp);
    bignum_rshift(&s, &temp, 1);
    bignum_assign(&s, &temp);
    t++;
  }

  int trials = 3;
  char buf[256];
  FILE* fp = fopen("/dev/urandom", "r");
  for (int i = 0; i < trials; i++) {
    struct bn a;
    struct bn tmp, tmp2;
    fread(buf, 256, 1, fp);
    bignum_from_string(&a, buf, 256);
    // a 是2 和 num - 1 之间的数字
    // a = (rand() % (num - 1 - 2)) + 2
    bignum_from_int(&tmp, 3);
    bignum_sub(num, &tmp, &tmp2);
    bignum_init(&tmp);
    bignum_mod(&a, &tmp2, &tmp);
    bignum_inc(&tmp);
    bignum_inc(&tmp);
    bignum_assign(&a, &tmp);
    // v = a ^ s % num
    struct bn v;
    bignum_init(&v);
    pow_mod_faster(&a, &s, num, &v);
    
    // if v != 1
    if (bignum_cmp(&v, &one) != 0) {
      uint64_t i = 0;
      struct bn num_sub_1;
      bignum_sub(num, &one, &num_sub_1);
      // while (v != num - 1):
      while (bignum_cmp(&v, &num_sub_1) != EQUAL) {
        if (i == t - 1)
          return 0;
        else {
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

int is_prime(struct bn *num) {
  struct bn t, t2;
  bignum_init(&t);
  bignum_init(&t2);
  bignum_from_int(&t, 2);
  if (bignum_cmp(num, &t) < 0)
    return 0;
  int low_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 
                    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
                    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
                    367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
                    503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
                    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
                    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971,
                    977, 983, 991, 997
  };

  for (int i = 0; i < sizeof(low_primes) / 4; i++) {
    bignum_from_int(&t, low_primes[i]);
    if (bignum_cmp(num, &t) == EQUAL)
      return 1;
  }

  for (int i = 0; i < sizeof(low_primes); i++) {
    bignum_from_int(&t, low_primes[i]);
    bignum_mod(num, &t, &t2);
    if(bignum_is_zero(&t2))
      return 0;
  }

  return miller_rabin(num);
}


static void test_rsa2048(void)
{
  char tmp[128];
  struct bn p, q;
  struct bn pub;
  FILE* fp = fopen("/dev/urandom", "r");
  fread(tmp, 128, 1, fp);
  bignum_from_string(&p, tmp, 128);
  printf("the p is : ");
  print_byte(tmp, 128);

  fread(tmp, 128, 1, fp);
  bignum_from_string(&q, tmp, 128);
  printf("the q is : ");
  print_byte(tmp, 128);

  fclose(fp);
  bignum_mul(&p, &q, &pub);
  char public[]  = "a15f36fc7f8d188057fc51751962a5977118fa2ad4ced249c039ce36c8d1bd275273f1edd821892fa75680b1ae38749fff9268bf06b3c2af02bbdb52a0d05c2ae2384aa1002391c4b16b87caea8296cfd43757bb51373412e8fe5df2e56370505b692cf8d966e3f16bc62629874a0464a9710e4a0718637a68442e0eb1648ec5";
  char private[] = "3f5cc8956a6bf773e598604faf71097e265d5d55560c038c0bdb66ba222e20ac80f69fc6f93769cb795440e2037b8d67898d6e6d9b6f180169fc6348d5761ac9e81f6b8879529bc07c28dc92609eb8a4d15ac4ba3168a331403c689b1e82f62518c38601d58fd628fcb7009f139fb98e61ef7a23bee4e3d50af709638c24133d";
  char buf[8192];
  bignum_to_string(&pub, buf, 8192);
  printf("the n is : %s", buf);
  

  struct bn n; /* public  key */
  struct bn d; /* private key */
  struct bn e; /* public exponent */
  struct bn m; /* clear text message */
  struct bn c; /* cipher text */

  //int len_pub = strlen(public);
  //int len_prv = strlen(private);

  int x = 54321;

  bignum_init(&n);
  bignum_init(&d);
  bignum_init(&e);
  bignum_init(&m);
  bignum_init(&c);

  bignum_from_string(&n, public,  256);
  bignum_from_string(&d, private, 256);
  bignum_from_int(&e, 65537);
  bignum_init(&m);
  bignum_init(&c);

  bignum_from_int(&m, x);
  bignum_to_string(&m, buf, sizeof(buf));
  printf("m = %s \n", buf);

//printf("  Copied %d bytes into m\n", i);

  printf("  Encrypting number x = %d \n", x);
  pow_mod_faster(&m, &e, &n, &c);
  printf("  Done...\n\n");

  bignum_to_string(&c, buf, sizeof(buf));
  printf("  Decrypting cipher text '");
  int i = 0;
  while (buf[i] != 0)
  {
    printf("%c", buf[i]);
    i += 1;
  }
  printf("'\n");

  /* Clear m */
  bignum_init(&m); 

  pow_mod_faster(&c, &d, &n, &m);
  printf("  Done...\n\n");


  bignum_to_string(&m, buf, sizeof(buf));
  printf("m = %s \n", buf);
}


int main()
{
  printf("\n");
  printf("Testing RSA encryption implemented with bignum. \n");

  test_rsa2048();

  printf("\n");
  printf("\n");

  int a;
  while (scanf("%d", &a) != EOF) {
    struct bn num;
    bignum_from_int(&num, a);
    if (miller_rabin(&num)) {
      printf("yes, it's prime\n");
    } else {
      printf("not prime\n");
    }

  }


  return 0;
}

