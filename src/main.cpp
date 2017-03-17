#include "prime_lib.h"
#include <iostream>

int main(int argc, char const *argv[]) {
  mpz_class D, p;
  p = "203956878356401977405765866929034577280193993314348263094772646453283062"
      "722701277632936616063144088173312372882677123879538709400158306567338328"
      "279154499698366071906766440037074217117805690872792848149112022286332144"
      "876183376326512083574821647933992961249"
      "917319836219304274280243803104015000563790123";
  generate_D(D, p);
  // std::cout << weak_lucas_prime(p, D) << '\n';
  std::cout << weak_lucas_prime(p, D) << '\n';

  // std::cout << prime_trial(p) << '\n';
  /*for (int i = 3; i < 11000; i++) {
    generate_D(D, p);
    if (weak_lucas_prime(i, D)) {
      if (!prime_trial(i)) {
        std::cout << "Wrong at:" << i << '\n';
      }
    }
  }*/
  // std::cout << miller_rabin(p) << '\n';

  return 0;
}
