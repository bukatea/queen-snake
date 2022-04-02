#ifndef BSMOOTH_SQUARE_FINDER_HPP
#define BSMOOTH_SQUARE_FINDER_HPP

#include <cmath>
#include <cstdint>

#include "common.hpp"
#include "tonelli_shanks.hpp"

class BSmoothSquareFinder {
	public:
		BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A);

		std::vector<bmp::mpz_int> find();

	private:
		static const bmp::mpf_float ERROR = "0.000000000000000000000000000000000000000000000000005"

		bmp::mpz_int n;
		bmp::mpz_int s;
		std::vector<std::size_t> factorBase;
		std::size_t piB;
		std::size_t A;

		std::vector<bmp::mpf_float_50> sieve;
};

#endif