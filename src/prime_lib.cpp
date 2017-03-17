#include "prime_lib.h"
#include <iostream>

/**
* Uses trial division to see if it's prime.
*/
bool prime_trial(const mpz_class p) {
  mpz_class up_to;
  i_sqrt(up_to, p);
  return (prime_trial(p, up_to));
}

/**
* Uses trial division up to a given value.
*/
bool prime_trial(const mpz_class p, const mpz_class up_to) {
  mpz_class i;
  for (i = 2; i <= up_to; i++) {
    if (p % i == 0)
      return false;
  }
  return true;
}

// TODO:
bool strong_lucas_prime(const mpz_class n, mpz_class D) {
  mpz_class delta_n(n + 1);

  return false;
}

bool weak_lucas_prime(const mpz_class n, mpz_class D) {
  mpz_class delta_n(n + 1);
  mpz_class P(1);
  mpz_class Q((1 - D) / 4);
  std::vector<mpz_class> Vals;
  Vals.push_back(1);

  size_t s = mpz_sizeinbase(delta_n.get_mpz_t(), 2) - 1;
  while (s > 0) {
    s--;
    mpz_class to_add;
    to_add = delta_n >> s;
    Vals.push_back(to_add & ~1);
    if ((to_add & 1) == 1)
      Vals.push_back(to_add);
  }

  mpz_class U, V;
  mpz_class new_U, new_V;
  U = V = new_U = new_V = 1;

  for (size_t i = 1; i < Vals.size(); i++) {
    if (Vals[i] == 2 * Vals[i - 1]) {
      new_U = (U * V) % n;

      mpz_class q_res;

      // to_power(q_res, Q, Vals[i - 1]);

      mod_pow(q_res, Q, Vals[i - 1], n);

      // SLOW
      new_V = (V * V - 2 * q_res) % n;

    } else {
      mpz_class u_res(P * U + V);
      new_U = ((u_res + (u_res & 1) * n) / 2) % n;

      mpz_class v_res(D * U + P * V);

      new_V = ((v_res + (v_res & 1) * n) / 2) % n;
    }
    U = new_U;
    V = new_V;
  }

  return U % n == 0;
}

void to_power(mpz_class &res, const mpz_class base, mpz_class exp) {
  res = 1;
  while (exp > 0) {
    res *= base;
    exp--;
  }
}

/**
* Generates the value D for which the jacobi function = -1
**/
void generate_D(mpz_class &D, const mpz_class p) {
  D = 5;
  for (int i = 0; i < 7; i++) {
    if (jacobi_symbol(D, p) == -1)
      return;
    D = -1 * (abs(D) + 2);
  }
  if (perfect_square(p)) {
    D = 0;
    return;
  }
  while (1) {
    if (jacobi_symbol(D, p) == -1)
      return;
    D = -1 * (abs(D) + 2);
  }
  D = 0;
}

/**
* Strong primality test in base 2. Returns true if p pseudoprime
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
* Calculates the jacobi symbol (a/n)
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
* Returns true if p is a perfect square
**/
bool perfect_square(const mpz_class p) {
  mpz_class res(0);
  i_sqrt(res, p);
  return (res * res == p);
}

/**
* Calculates the integer square root
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
* Adapted from
* Schneier, Bruce (1996). Applied Cryptography: Protocols, Algorithms, and
* Source Code in C, Second Edition (2nd ed.). Wiley. ISBN 978-0-471-11709-4.
**/
void mod_pow(mpz_class &res, const mpz_class base, const mpz_class exp,
             const mpz_class mod) {
  // base %= modulus;
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
