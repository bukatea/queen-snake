#include "quadratic_sieve.hpp"

std::vector<bmp::mpz_int> quadraticSieve(const bmp::mpz_int &n, std::size_t A, std::size_t B, bool incrementalChecking) {
    std::vector<std::size_t> fB = factorBase(B);

    std::size_t required = fB.size() + 1;
    BSmoothSquareFinder finder(n, fB, A);
    BSmoothSolution solution = finder.find(required);
    if (solution.nDivisor != 0) {
        return std::vector<bmp::mpz_int>{solution.nDivisor, n / solution.nDivisor};
    }
    if (solution.bSmoothSquares.size() < required) {
        // incrementalChecking
        if (incrementalChecking) {

        }
    } else {
        // matrix stuff
    }
}
