#ifndef BSMOOTH_SQUARE_FINDER_HPP
#define BSMOOTH_SQUARE_FINDER_HPP

#include <cstdint>
#include <gmpxx.h>

class BSmoothSquareFinder {
	public:
		BSmoothSquareFinder(const mpz_class &n, const std::vector<std::size_t> &factorBase, std::size_t A);

		std::vector<> find()

	private:
		mpz_class n;
		std::vector<std::size_t> factorBase;
		std::size_t A;

		std::vector<std::size_t> sieve;
};

#endif