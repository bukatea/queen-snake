#ifndef QUADRATIC_SIEVE
#define QUADRATIC_SIEVE

#include <Eigen/Dense>
#include <numeric>
#include <iostream>
#include <vector>

#include "bsmooth_square_finder.hpp"
#include "common.hpp"
#include "util.hpp"

std::vector<bmp::mpz_int> quadraticSieve(const bmp::mpz_int &n, std::size_t A, std::size_t B);

#endif