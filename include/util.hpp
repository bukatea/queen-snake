#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdint>
#include <gmpxx.h>
#include <vector>

// Pingela's Algorithm variant
mpz_class powerMod(const mpz_class &base, const mpz_class &exponent, const mpz_class &modulus);

// Sieve of Eratosthenes
std::vector<std::size_t> factorBase(std::size_t B);

#endif