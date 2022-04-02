#ifndef QUADRATIC_SIEVE
#define QUADRATIC_SIEVE

#include <vector>

#include "common.hpp"

std::vector<mpz_class> quadraticSieve(const bmp::mpz_int &n, std::size_t B, bool incrementalChecking);

#endif