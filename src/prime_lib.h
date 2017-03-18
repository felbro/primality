#ifndef PRIME_LIB
#define PRIME_LIB

#include "prime_1000.h"
#include <cmath>
#include <gmpxx.h>

bool prime_trial(const mpz_class p);
bool prime_trial(const mpz_class p, const mpz_class up_to);
bool prime_trial_1000(const mpz_class p);

bool generate_D(mpz_class &D, const mpz_class p);

bool miller_rabin(const mpz_class p);
int jacobi_symbol(mpz_class a, mpz_class n);

bool perfect_square(const mpz_class p);
void i_sqrt(mpz_class &res, const mpz_class p);

void mod_pow(mpz_class &res, const mpz_class base, const mpz_class exp,
             const mpz_class mod);

bool strong_lucas_prime(const mpz_class n);

#endif
