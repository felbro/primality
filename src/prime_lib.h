#ifndef PRIME_LIB
#define PRIME_LIB

#include <cmath>
#define ULLI unsigned long long int

bool prime_trial(ULLI p);
bool prime_trial(ULLI p, ULLI up_to);
int jacobi_symbol(long long int a, long long int n);
int mod_pow(int b, ULLI e, int m);

#endif
