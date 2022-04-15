#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdint>
#include <Eigen/Dense>
#include <vector>

#include "common.hpp"

// Sieve of Eratosthenes
std::vector<std::size_t> factorBase(std::size_t B);

struct PrimeFactorization {
    Eigen::VectorXi exponents;
    Eigen::VectorXi exponentsMod2;

    PrimeFactorization(const Eigen::VectorXi &exponents, const Eigen::VectorXi &exponentsMod2);
};

PrimeFactorization trialDivision(const std::vector<std::size_t> &factorBase, const bmp::mpz_int &n);

#endif