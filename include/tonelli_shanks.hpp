#ifndef TONELLI_SHANKS_HPP
#define TONELLI_SHANKS_HPP

#include <gmpxx.h>
#include "util.hpp"

// Courtesy of https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#C
struct QuadraticResidue {
	mpz_class root1;
    mpz_class root2;
    bool exists;

	QuadraticResidue(const mpz_class &root1, const mpz_class &root2, bool exists);
};
 
QuadraticResidue tonelliShanks(const mpz_class &n, const mpz_class &p);

#endif