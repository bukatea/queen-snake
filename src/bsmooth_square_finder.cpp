#include "bsmooth_square_finder.hpp"

BSmoothSquareFinder::BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A)
	: n(n), s(static_cast<bmp::mpz_int>(std::sqrt(n)) + 1), factorBase(factorBase), piB(factorBase.size()), A(A), sieve(A + 1) {
	// initialize sieve with ln f(x)
	// sieve will look like [0, 1, ..., A]
	// add negatives?
	for (std::size_t x = 0; x < A + 1; x++) {
		bmp::mpz_int base = s + x;
		sieve[x] = std::log(static_cast<bmp::mpf_float_50>(base * base - n));
	}
}

std::vector<bmp::mpz_int> BSmoothSquareFinder::find() {
	std::vector<bmp::mpz_int> res;
	for (std::size_t p : factorBase) {
		QuadraticResidueSquareRoot res = tonelliShanks(n, p);
		if (!res.exists) {
			continue;
		}

		bmp::mpf_float_50 logp = std::log(static_cast<bmp::mpf_float_50>(p));

		bmp::mpz_int x = res.root1 - s;
		while (x < A + 1) {
			sieve[x] -= logp;
			if (-piB * ERROR <= sieve[x] && sieve[x] <= piB * ERROR) {
				// sieve[x] is approx. 0, B-smooth number found
				res.push_back(x + s);
			}
			x += p;
		}
		if (p != 2) {
			x = res.root2 - s;
			while (x < A + 1) {
				sieve[x] -= logp;
				if (-piB * ERROR <= sieve[x] && sieve[x] <= piB * ERROR) {
					// sieve[x] is approx. 0, B-smooth number found
					res.push_back(x + s);
				}
				x += p;
			}
		}
	}

	return res;
}
