#include "bsmooth_square_finder.hpp"

BSmoothSquareFinder::BSmoothSquareFinder(const mpz_class &n, const std::vector<std::size_t> &factorBase, std::size_t A) : n(n), factorBase(factorBase), A(A), sieve(A + 1) {
	// initialize sieve with ln f(x)
	// sieve will look like [0, 1, ..., A]
	const mpz_class s = sqrt(n) + 1;
}
