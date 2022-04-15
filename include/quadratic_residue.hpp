#ifndef QUADRATIC_RESIDUE_HPP
#define QUADRATIC_RESIDUE_HPP

#include <utility>

#include "common.hpp"
#include "util.hpp"

// Courtesy of https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#C
struct QuadraticResidueSquareRoot {
    std::vector<bmp::mpz_int> roots;
    bool exists;

    QuadraticResidueSquareRoot() : exists(false) {}
    QuadraticResidueSquareRoot(const std::vector<bmp::mpz_int> &roots, bool exists);
    QuadraticResidueSquareRoot(const std::vector<bmp::mpz_int> &&roots, bool exists);
};

QuadraticResidueSquareRoot trialSquareRoot(const bmp::mpz_int &n, std::size_t k);

QuadraticResidueSquareRoot tonelliShanks(const bmp::mpz_int &n, std::size_t p);

#endif