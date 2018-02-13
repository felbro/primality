#pragma once
/**
* Header file for prime_lib.cpp which is a primality test library
* for a Braille-PSW test along with
* necessary functions and tests, such as Miller-Rabin and
* Lucas primality tests. The Braille-PSW test is 100%
* accurate for values < 2^64 and has no known pseudoprimes for
* values < 2^128. For larger values, the Braille-PSW test
* is to be used as a probability test since there may
* be pseudoprimes. However, if any of the tests
* (baillie_psw, miller_rabin, strong_lucas_prime) returns
* False, the value is certain to be False.
*
* @author Felix Broberg
* @version 2017-03-19
*/

#include <gmpxx.h>
#include <cmath>
#include "prime_1000.h"

bool prime_trial(const mpz_class p);
bool prime_trial(const mpz_class p, const mpz_class up_to);
bool prime_trial_1000(const mpz_class p);

bool baillie_psw(const mpz_class p);
bool strong_lucas_prime(const mpz_class n);
bool miller_rabin(const mpz_class p);

bool generate_D(mpz_class &D, const mpz_class p);
int jacobi_symbol(mpz_class a, mpz_class n);

bool perfect_square(const mpz_class p);
void i_sqrt(mpz_class &res, const mpz_class p);
void mod_pow(mpz_class &res, const mpz_class base, const mpz_class exp,
             const mpz_class mod);
