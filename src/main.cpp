/**
* Main file for a primality tester. Reads input and tries to parse
* into a number. Input is either read from std::cin or from a variable
* input (see marked lines in main function).
*
* Performs checks to ensure input > 0.
* Calls on a Braille-PSW test through baillie_psw() for input > 2.
* Values 1 and 2 are defined as prime from the start, and are not
* tested.
*
* @author Felix Broberg
* @version 2017-03-19
*/

#include "prime_lib.h"
#include <iostream>
#include <regex>
#include <string>

/**
* Parses an input string by removing everything but numbers and
* a possible "-" in the beginning of the string. Stores the
* result in p.
*
* @param p  Where to store the parsed input.
* @param input Input string to be parsed.
*/
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

/**
* Prints an input string. If p < 2^64, an is_prime input
* is certain to true and is printed as such. Else, if is_prime
* true but p >= 2^64, the printed result is said to be "likely"
* instead of certain.
*
* @param input  input value to be printed again along with answer.
* @param p      value that has been tested.
* @param is_prime True if p is a prime number. False otherwise.
*/
void prime_print(std::string input, const mpz_class p, bool is_prime) {
  mpz_class max(std::to_string(ULLONG_MAX));
  if (is_prime && p <= max)
    std::cout << input << " is a prime." << '\n';
  else if (is_prime)
    std::cout << input << " is most likely a prime." << '\n';
  else
    std::cout << input << " is not a prime." << '\n';
}

/**
* Reads input and on a Braille-PSW test through baillie_psw() for input > 2.
* Values 1 and 2 are defined as prime from the start, and are not
* tested. The results are then printed.
*/
int main(int argc, char const *argv[]) {

  mpz_class p;
  std::string input = "";
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
    prime_print(input, p, baillie_psw(p));
  return 0;
}
