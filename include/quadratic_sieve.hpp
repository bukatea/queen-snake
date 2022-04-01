#ifndef QUADRATIC_SIEVE
#define QUADRATIC_SIEVE

#include <gmpxx.h>
#include <vector>

std::vector<mpz_class> quadraticSieve(const mpz_class &n, std::size_t B, bool incrementalChecking);

#endif