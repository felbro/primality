#ifndef PRIME_LIB
#define PRIME_LIB

#include <cmath>
#define ULLI unsigned long long int

bool prime_trial(ULLI p);
bool prime_trial(ULLI p, ULLI up_to);

int generate_D(ULLI p);

bool perfect_square(ULLI p);
ULLI i_sqrt(ULLI p);

int jacobi_symbol(long long int a, long long int n);
ULLI mod_pow(ULLI b, ULLI e, ULLI m);

bool miller_rabin(ULLI p);

#endif
