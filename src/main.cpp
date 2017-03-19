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

void prime_print(std::string input, const mpz_class p, bool is_prime) {
  mpz_class max(std::to_string(ULLONG_MAX));
  if (is_prime && p <= max)
    std::cout << input << " is a prime." << '\n';
  else if (is_prime)
    std::cout << input << " is most likely a prime." << '\n';
  else
    std::cout << input << " is not a prime." << '\n';
}

int main(int argc, char const *argv[]) {

  mpz_class p;
  std::string input;
  /**********************************/
  /* Please, use only one of these: */
  std::getline(std::cin, input);
  // input = "5";
  /**********************************/

  parse_input(p, input);
  std::cout << "" << '\n';
  if (p <= 0)
    std::cout << "Only positive integers can be prime." << '\n';
  else if (p == 1 || p == 2)
    prime_print(input, p, true);
  else
    prime_print(input, p, braille_psw(p));
  return 0;
}
