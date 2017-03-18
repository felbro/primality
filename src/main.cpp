#include "prime_lib.h"
#include <iostream>
#include <regex>
#include <string>

void parse_input(mpz_class &p, std::string input) {

  std::regex not_a_number("[^\\d]");
  std::regex not_a_sign("[^\\d\\-]");
  std::string res = "";
  regex_replace(std::back_inserter(res), input.begin(), input.end(), not_a_sign,
                "");
  input = "";
  regex_replace(std::back_inserter(input), res.begin() + 1, res.end(),
                not_a_number, "");
  input = res[0] + input;
  p = input;
}

int main(int argc, char const *argv[]) {

  mpz_class p;
  std::string input;
  std::getline(std::cin, input);
  parse_input(p, input);

  if (p <= 0)
    std::cout << "A prime has to be positive" << '\n';
  else if (p == 1 || p == 2)
    std::cout << "prime" << '\n';
  else if (!prime_trial_1000(p))
    std::cout << "prime trial" << '\n';
  else if (!miller_rabin(p))
    std::cout << "miller rabin" << '\n';
  else if (!strong_lucas_prime(p))
    std::cout << "lucas" << '\n';
  else
    std::cout << "most likely a prime" << '\n';
  return 0;
}
