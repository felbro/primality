#include "prime_lib.h"
#include <iostream>

/**
* Uses trial division to see if it's prime.
*/
bool prime_trial(ULLI p) { return (prime_trial(p, sqrt(p))); }

/**
* Uses trial division up to a given value.
*/
bool prime_trial(ULLI p, ULLI up_to) {
  for (ULLI i = 2; i <= up_to; i++) {
    if (p % i == 0)
      return false;
  }
  return true;
}

/**
* Strong primality test in base 2. Returns true if p pseudoprime
*/
bool miller_rabin(ULLI p) {
  ULLI d = p - 1;
  ULLI s = 0;
  while (d % 2 == 0) {
    d >>= 1;
    s++;
  }
  if (mod_pow(2, d, p) == 1) {
    return true;
  }
  for (ULLI r = 0; r < s; r++) {
    if (mod_pow(2, d * pow(2, r), p) == p - 1)
      return true;
  }
  return false;
}

/**
* Generates the value D for which the jacobi function = -1
**/
int generate_D(ULLI p) {
  int D = 5;
  for (int i = 0; i < 7; i++) {
    if (jacobi_symbol(D, p) == -1)
      return D;
    D = -1 * (std::abs(D) + 2);
  }
  if (perfect_square(p))
    return 0;

  while (1) {
    if (jacobi_symbol(D, p) == -1)
      return D;
    D = -1 * (std::abs(D) + 2);
  }
  return 0;
}

/**
* Returns true if p is a perfect square
**/
bool perfect_square(ULLI p) {
  ULLI res = i_sqrt(p);
  return (res * res == p);
}

/**
* Calculates the integer square root
**/
ULLI i_sqrt(ULLI p) {
  long long int shift = 2;
  ULLI p_shift = p >> shift;
  while (p_shift != 0 && p_shift != p) {
    shift += 2;
    p_shift = p >> shift;
  }
  shift -= 2;
  ULLI res = 0;

  while (shift >= 0) {
    res = res << 1;
    ULLI c_res = res + 1;
    if (c_res * c_res <= p >> shift) {
      res = c_res;
    }
    shift -= 2;
  }
  return res;
}

/**
* Calculates the jacobi symbol (a/n)
*/
int jacobi_symbol(long long int a, long long int n) {
  int mult = 1;
  if (a < 0) {
    a *= -1;
    if ((n % 4) == 3)
      mult *= -1;
  }
  while (a != 0) {
    while (a % 2 == 0) {
      a /= 2;
      char mod_8 = (n % 8);
      if (mod_8 == 3 || mod_8 == 5)
        mult *= -1;
    }
    long long int t = a;
    a = n;
    n = t;

    if ((a % 4) == 3 && (n % 4) == 3)
      mult *= -1;
    a %= n;
  }
  if (n == 1)
    return mult;

  return 0;
}

/**
* Adapted from
* Schneier, Bruce (1996). Applied Cryptography: Protocols, Algorithms, and
* Source Code in C, Second Edition (2nd ed.). Wiley. ISBN 978-0-471-11709-4.
**/
ULLI mod_pow(ULLI base, ULLI exp, ULLI modulus) {
  base %= modulus;
  ULLI result = 1;
  while (exp > 0) {
    if (exp & 1)
      result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp >>= 1;
  }
  return result;
}
