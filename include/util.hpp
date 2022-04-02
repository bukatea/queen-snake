#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdint>
#include <vector>

#include "common.hpp"

// Pingela's Algorithm variant
bmp::mpz_int powerMod(const bmp::mpz_int &base, const bmp::mpz_int &exponent, const bmp::mpz_int &modulus);

// Sieve of Eratosthenes
std::vector<std::size_t> factorBase(std::size_t B);

#endif