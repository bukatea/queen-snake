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

// add support for multiple finds
class BSmoothSquareFinder {
    public:
        BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A);

        std::vector<BSmoothSquare> find(std::size_t numBSmoothSquares);

    private:
        static const bmp::mpf_float ERROR = "0.000000000000000000000000000000000000000000000000005"

        bmp::mpz_int n;
        bmp::mpz_int s;
        std::vector<std::size_t> factorBase;
        std::size_t piB;
        std::size_t A;

        std::vector<bmp::mpf_float_50> sieve;
        std::vector<std::size_t> primeExponentsLeqA;
        std::vector<std::vector<std::size_t>> primeFactors;
};

#endif