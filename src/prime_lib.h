#ifndef PRIME_LIB
#define PRIME_LIB

#include <cmath>
#include <gmpxx.h>

bool prime_trial(const mpz_class p);
bool prime_trial(const mpz_class p, const mpz_class up_to);
// works
void generate_D(mpz_class &D, const mpz_class p);
// Works
bool miller_rabin(const mpz_class p);
int jacobi_symbol(mpz_class a, mpz_class n);

// works
bool perfect_square(const mpz_class p);
void i_sqrt(mpz_class &res, const mpz_class p);

// works
void mod_pow(mpz_class &res, const mpz_class base, const mpz_class exp,
             const mpz_class mod);
#endif
