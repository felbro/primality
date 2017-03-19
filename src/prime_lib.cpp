/**
* Primality test library for a Braille-PSW test along with
* necessary functions and tests, such as Miller-Rabin and
* Lucas primality tests. The Braille-PSW test is 100%
* accurate for values < 2^64 and has no known pseudoprimes for
* values < 2^128. For larger values, the Braille-PSW test
* is to be used as a probability test since there may
* be pseudoprimes. However, if any of the tests
* (braille_psw, miller_rabin, strong_lucas_prime) returns
* False, the value is certain to be False.
*
* @author Felix Broberg
* @version 2017-03-19
*/

#include "prime_lib.h"

/**
* Naive trial division to test primality
*
* @param p  the number to test
* @return True if p is prime. False otherwise
*/
bool prime_trial(const mpz_class p) {
  mpz_class up_to;
  i_sqrt(up_to, p);
  return (prime_trial(p, up_to));
}

/**
* Naive trial division to test primality up to a given value.
*
* @param p  the number to test
* @param up_to  Up to which value to perfrom trial division.
* @return True if p is prime. False otherwise
*/
bool prime_trial(const mpz_class p, const mpz_class up_to) {
  mpz_class i;
  for (i = 2; i <= up_to; i++) {
    if (p % i == 0)
      return false;
  }
  return true;
}

/**
* Braille-PSW implementation. Performs trial divisions
* up to 1000, a Miller-Rabin test as well as a strong
* Lucas primality test. Will produce a correct result for
* all values < 2^64, and for values < 2^128 no pseudoprimes
* have been found.
*
* @param p  the number to test
* @return True if p is a probable prime, False if it is not a prime.
*/
bool braille_psw(const mpz_class p) {
  if (!prime_trial_1000(p))
    return false;
  if (p < f_k_primes[F_K_PRIMES_SIZE - 1] * f_k_primes[F_K_PRIMES_SIZE - 1])
    return true;
  if (!miller_rabin(p))
    return false;
  if (!strong_lucas_prime(p))
    return false;
  return true;
}

/**
* Trial division using primes under 1000.
* If p < 1000, a list is searched for the value.
* if p >= 1000, trial division is used.
*
* @param p  the number to test
* @return True if p could be a prime based on the
*         primes < 1000. False if it cannot be a prime.
*         If p <= 994009, a return value "true" is a certain
*         true. Else it is only probable.
*/
bool prime_trial_1000(const mpz_class p) {
  if ((p & 1) == 1 && p < 1000) {
    return std::binary_search(f_k_primes, f_k_primes + F_K_PRIMES_SIZE, p);
  }
  for (size_t i = 0; i < F_K_PRIMES_SIZE; i++) {
    if (p % f_k_primes[i] == 0)
      return false;
  }
  return true;
}

/**
* Strong lucas primality test.
* Will, for most odd values, return a correct result.
* The first strong lucas pseudoprimes are:
* 5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519.
*
* @param n  the number to test
* @return true if n is a probable prime. False if
*         it is not a prime.
*/
bool strong_lucas_prime(const mpz_class n) {
  mpz_class D;
  if (!generate_D(D, n))
    return false;
  mpz_class delta_n(n + 1);
  mpz_class P(1);
  mpz_class Q((1 - D) / 4);
  mpz_class Qk(Q);
  // Split delta_n into delta_n = d*2^s
  size_t s = 0;
  while ((delta_n & 1) == 0) {
    delta_n >>= 1;
    s++;
  }

  size_t d = mpz_sizeinbase(delta_n.get_mpz_t(), 2);

  mpz_class U, V;
  U = 1;
  V = 1;
  // For U_d = 0 (mod n)
  for (size_t i = 2; i <= d; i++) {
    U = (U * V) % n;
    V = (V * V - 2 * Qk) % n;
    Qk = (Qk * Qk) % n;
    // Current bit == 1
    if (((delta_n >> (d - i)) & 1) == 1) {
      mpz_class u_res(P * U + V);
      mpz_class new_U(((u_res + (u_res & 1) * n) / 2) % n);

      mpz_class v_res(D * U + P * V);
      mpz_class new_V(((v_res + (v_res & 1) * n) / 2) % n);

      U = new_U;
      V = new_V;
      Qk = (Qk * Q) % n;
    }
  }
  // For U_d = 0 (mod n)
  if (U == 0 || V == 0)
    return true;
  // For V_(d*2^r) = 0 (mod n), where 0 <= r < s.
  for (size_t i = 0; i < s; i++) {
    V = (V * V - 2 * Qk) % n;
    Qk = (Qk * Qk) % n;
    if (V == 0)
      return true;
  }
  return false;
}

/**
* Miller-Rabin primality test for base 2.
* Higher rate of pseudoprimes than the strong
* lucas primality test, but they have no known
* overlap, making it a suitable combination.
*
* @param p  the number to test
* @return true if n is a probable prime. False if
*         it is not a prime.
*/
bool miller_rabin(const mpz_class p) {
  mpz_class d(p - 1);
  unsigned long long s = 0;
  while (d % 2 == 0) {
    d >>= 1;
    s++;
  }
  mpz_class res;
  mod_pow(res, 2, d, p);
  if (res == 1) {
    return true;
  }
  unsigned long long r;
  for (r = 0; r < s; r++) {
    mod_pow(res, 2, d * pow(2, r), p);
    if (res == p - 1)
      return true;
  }
  return false;
}

/**
* Generates the value D in (D/p) such that (D/p) = -1.
* (D/p) is the jacobi symbol.
*
* @param D  Reference to the generated D.
* @param p  Value to be tested
* @return true if a value D can be found. Else,
*         p will be a perfect square.
**/
bool generate_D(mpz_class &D, const mpz_class p) {
  D = 5;
  for (int i = 0;; i++) {
    if (jacobi_symbol(D, p) == -1)
      return true;
    if (D < 0)
      D = -1 * D + 2;
    else
      D = -1 * (D + 2);

    if (i == 5 && perfect_square(p))
      break;
  }

  D = 0;
  return false;
}

/**
* Calculates the jacobi symbol (a/n) defined as
* (a/p) = 0 if a = 0 (mod p)
* (a/p) = 1 if a != 0 (mod p) and a = x^2 (mod p)
*                              for some integer x
* (a/p) = -1 if a != 0 (mod p) and no such x exists
*
* As the jacobi symbol is only defined for odd primes p,
* some further calculations have to be made.
*
* @param a  the numerator
* @param n  the denominator
* @return the jacobi symbol
*/
int jacobi_symbol(mpz_class a, mpz_class n) {
  int mult = 1;
  if (a < 0) {
    a *= -1;
    if ((n % 4) == 3)
      mult *= -1;
  }
  while (a != 0) {
    while (a % 2 == 0) {
      a >>= 1;
      mpz_class mod_8(n % 8);
      if (mod_8 == 3 || mod_8 == 5)
        mult *= -1;
    }
    swap(a, n);

    if ((a % 4) == 3 && (n % 4) == 3)
      mult *= -1;
    a %= n;
  }
  if (n == 1)
    return mult;
  return 0;
}

/**
* Checks if the value p is a perfect square
*
* @param p  the value to test
* @return True if p is a perfect square.
**/
bool perfect_square(const mpz_class p) {
  mpz_class res(0);
  i_sqrt(res, p);
  return (res * res == p);
}

/**
* Calculates the integer square root of p and
* puts the result into res.
*
* @param res  reference to the "return" value
*             - the integer square root of p.
* @param p  the value to test
**/
void i_sqrt(mpz_class &res, const mpz_class p) {
  long long int shift = 2;
  mpz_class p_shift(p >> shift);
  while (p_shift != 0 && p_shift != p) {
    shift += 2;
    p_shift = p >> shift;
  }
  shift -= 2;
  res = 0;

  while (shift >= 0) {
    res = res << 1;
    mpz_class c_res(res + 1);
    if (c_res * c_res <= p >> shift) {
      res = c_res;
    }
    shift -= 2;
  }
}

/**
* Calculates modded power of a value.
*
* Adapted from
* Schneier, Bruce (1996). Applied Cryptography: Protocols, Algorithms, and
* Source Code in C, Second Edition (2nd ed.). Wiley. ISBN 978-0-471-11709-4.
*
* @param res  reference to the result - the modded power of a value.
* @param base the base value.
* @param exp  the exponent value.
* @param mod  the mod value.
**/
void mod_pow(mpz_class &res, const mpz_class base, const mpz_class exp,
             const mpz_class mod) {
  mpz_class m_base(base % mod);
  mpz_class m_exp(exp);
  res = 1;
  while (m_exp > 0) {
    if ((m_exp & 1) == 1)
      res = (res * m_base) % mod;
    m_base = (m_base * m_base) % mod;
    m_exp >>= 1;
  }
}
