#ifndef TONELLI_SHANKS_HPP
#define TONELLI_SHANKS_HPP

#include "util.hpp"
#include <cstdint>

// Courtesy of https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#C
struct QuadraticResidue {
	std::uint64_t root1;
    std::uint64_t root2;
    bool exists;

	QuadraticResidue(std::uint64_t root1, std::uint64_t root2, bool exists) : root1(root1), root2(root2), exists(exists) {}
};
 
QuadraticResidue tonelliShanks(std::uint64_t n, std::uint64_t p);

#endif