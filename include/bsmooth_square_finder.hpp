#ifndef BSMOOTH_SQUARE_FINDER_HPP
#define BSMOOTH_SQUARE_FINDER_HPP

#include <cmath>
#include <cstdint>
#include <Eigen/Dense>
#include <iostream>

#include "common.hpp"
#include "quadratic_residue.hpp"

struct BSmoothSolution {
    // should be 0 if p^k does not divide n (this is most cases)
    // otherwise will be the divisor of n (this is good, but will be very infrequent)
    bmp::mpz_int nDivisor;
    // size is numBSmoothSquares
    std::vector<bmp::mpz_int> bSmoothSquares;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> exponentMatrix;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> exponentMod2Matrix;

    BSmoothSolution(const bmp::mpz_int &nDivisor, const std::vector<bmp::mpz_int> &bSmoothSquares, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &exponentMatrix, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &exponentMod2Matrix);
};

// add support for multiple finds
class BSmoothSquareFinder {
    public:
        BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A);

        BSmoothSolution find(std::size_t numBSmoothSquares);

    private:
        bool findSolution(
            const QuadraticResidueSquareRoot &qrsr,
            std::size_t divisor,
            std::size_t numBSmoothSquares,
            std::size_t i,
            std::size_t k,
            std::vector<bmp::mpz_int> &res
        );

        bmp::mpz_int n;
        bmp::mpz_int s;
        std::vector<std::size_t> factorBase;
        std::size_t piB;
        std::size_t A;
        bmp::mpz_int bSmoothBound;

        std::vector<bmp::mpz_int> sieve;
        std::vector<std::size_t> primeExponentsLeqA;
        std::vector<QuadraticResidueSquareRoot> qrsrModP;
        std::vector<std::vector<std::size_t>> primePowers;

        std::size_t currentRow;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> exponentMatrix;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> exponentMod2Matrix;
};

#endif