#ifndef BSMOOTH_SQUARE_FINDER_HPP
#define BSMOOTH_SQUARE_FINDER_HPP

#include <cmath>
#include <cstdint>

#include "common.hpp"
#include "tonelli_shanks.hpp"

struct BSmoothSquare {
    bmp::mpz_int x;
    // size is piB
    std::vector<std::size_t> primeFactors;

    BSmoothSquare(const bmp::mpz_int &x, const std::vector<std::size_t> &primeFactors);
};

struct BSmoothSolution {
    // should be 0 if p^k does not divide n (this is most cases)
    // otherwise will be the divisor of n (this is good, but will be very infrequent)
    bmp::mpz_int nDivisor;
    // size is numBSmoothSquares
    std::vector<BSmoothSquare> bSmoothSquares;

    BSmoothSolution(const bmp::mpz_int &nDivisor, const std::vector<BSmoothSquare> &bSmoothSquares);
};

// add support for multiple finds
class BSmoothSquareFinder {
    public:
        BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A);

        BSmoothSolution find(std::size_t numBSmoothSquares);

    private:
        static const bmp::mpf_float ERROR = "0.000000000000000000000000000000000000000000000000005"

        bmp::mpz_int n;
        bmp::mpz_int s;
        std::vector<std::size_t> factorBase;
        std::size_t piB;
        std::size_t iterations;
        std::size_t A;

        std::vector<bmp::mpf_float_50> sieve;
        std::vector<std::size_t> primeExponentsLeqA;
        std::vector<std::vector<std::size_t>> primeFactors;
};

#endif