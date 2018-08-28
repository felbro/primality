# Primality test
## Intro
An implementation of the Baillie-PSW primality test. The Braille-PSW test is 100% accurate for values < 2^64 and has no known pseudoprimes for any values. It has been theorised that it's, infact, 100% correct for all numbers up to atleast 11 000 digits . The test is performed using a combination of trial division, Miller-Rabin primality test for base 2 and a strong Lucas primality test.

## Build
To build, simply type "make" when standing in the "/src" directory. A requirement is that the "GNU Multiple Precision Arithmetic Library" is installed. This is however done by default on many Linux-distros.

## Run
Both input modes require the value to be on a single line. Commas, spaces etc are allowed as the parser will remove them before computing.
### From standard input
Execute the generated file (usually done through ./main). The program will await an input and then compute the result. If the value to be checked is greater than 1500 digits, please refer to the input option below.

### From file
Execute the generated file with the prime piped through. This can easily be done through "cat file_name | ./main".
