#ifndef TONELLI_SHANKS_HPP
#define TONELLI_SHANKS_HPP

#include "common.hpp"
#include "util.hpp"

// Courtesy of https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#C
struct QuadraticResidueSquareRoot {
    bmp::mpz_int root1;
    bmp::mpz_int root2;
    bool exists;

    QuadraticResidueSquareRoot(const bmp::mpz_int &root1, const bmp::mpz_int &root2, bool exists);
};
 
QuadraticResidueSquareRoot tonelliShanks(const bmp::mpz_int &n, std::size_t p);

#endif