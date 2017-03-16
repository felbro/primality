#include "prime_lib.h"

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
int mod_pow(int base, ULLI exp, int modulus) {
  base %= modulus;
  int result = 1;
  while (exp > 0) {
    if (exp & 1)
      result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp >>= 1;
  }
  return result;
}
