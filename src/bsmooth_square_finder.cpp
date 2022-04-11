#include "bsmooth_square_finder.hpp"

BSmoothSquare::BSmoothSquare(const bmp::mpz_int &x, const std::vector<int> &primeFactors) : x(x), primeFactors(primeFactors) {}

BSmoothSolution::BSmoothSolution(const bmp::mpz_int &nDivisor, const std::vector<BSmoothSquare> &bSmoothSquares, const Eigen::MatrixXi &exponentMatrix)
    : nDivisor(nDivisor),
      bSmoothSquares(bSmoothSquares),
      exponentMatrix(exponentMatrix) {}

BSmoothSquareFinder::BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A)
    : n(n),
      s(static_cast<bmp::mpz_int>(bmp::sqrt(n)) + 1),
      factorBase(factorBase),
      piB(factorBase.size()),
      iterations(0),
      A(A),
      sieve(A + 1),
      primeExponentsLeqA(piB),
      primeFactors(A + 1, std::vector<int>(piB)),
      primeFactorsMod2(A + 1, Eigen::VectorXi::Zero(piB)),
      currentCol(0) {
    // initialize sieve with ln f(x)
    // sieve will look like [0, 1, ..., A]
    // also initialize primeExponentsLeqA, primeExponentsLeqA[i] is the largest n s.t. factorBase[i]^n <= A
    // add negatives?
    for (std::size_t x = 0; x < A + 1; x++) {
        bmp::mpz_int base = s + x;
        sieve[x] = bmp::log(static_cast<bmp::mpf_float_50>(base * base - n));
    }

    for (std::size_t i = 0; i < piB; i++) {
        primeExponentsLeqA[i] = static_cast<std::size_t>(std::log(A) / std::log(factorBase[i]));
        iterations += primeExponentsLeqA[i];
    }
}

BSmoothSolution BSmoothSquareFinder::find(std::size_t numBSmoothSquares) {
    const static bmp::mpf_float SIEVE_ERROR = iterations * bmp::mpf_float{"0.000000000000000000000000000000000000000000000000005"};

    exponentMatrix.resize(numBSmoothSquares, piB);

    std::vector<BSmoothSquare> res;
    for (std::size_t i = 0; i < piB; i++) {
        for (std::size_t k = primeExponentsLeqA[i]; k >= 1; k--) {
            std::size_t divisor = std::pow(factorBase[i], k);

            // check if p^k | n by chance
            if (n % divisor == 0)
                return BSmoothSolution(divisor, res, exponentMatrix);

            QuadraticResidueSquareRoot roots = tonelliShanks(n, divisor);
            if (!roots.exists) {
                continue;
            }

            bmp::mpf_float_50 logd = bmp::log(static_cast<bmp::mpf_float_50>(divisor));

            bmp::mpz_int x = roots.root1 - s;
            if (x < 0)
                x = (x % divisor + divisor) % divisor;
            std::size_t indexX = static_cast<std::size_t>(x);
            while (indexX < A + 1) {
                sieve[indexX] -= logd;
                primeFactors[indexX][i] += k;
                primeFactorsMod2[indexX][i] += k % 2;
                if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                    // sieve[x] is approx. 0, B-smooth number found
                    res.emplace_back(x + s, primeFactors[indexX]);
                    exponentMatrix.col(currentCol++) = primeFactorsMod2[indexX];
                    if (res.size() == numBSmoothSquares) {
                        return BSmoothSolution(0, res, exponentMatrix);
                    }
                }
                x += divisor;
                indexX += divisor;
            }
            if (factorBase[i] != 2) {
                // if the prime is 2 the roots should be the same, 1 (mod 2)
                x = roots.root2 - s;
                if (x < 0)
                    x = (x % divisor + divisor) % divisor;
                indexX = static_cast<std::size_t>(x);
                while (indexX < A + 1) {
                    sieve[indexX] -= logd;
                    primeFactors[indexX][i] += k;
                    primeFactorsMod2[indexX][i] += k % 2;
                    if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                        // sieve[x] is approx. 0, B-smooth number found
                        res.emplace_back(x + s, primeFactors[indexX]);
                        exponentMatrix.col(currentCol++) = primeFactorsMod2[indexX];
                        if (res.size() == numBSmoothSquares) {
                            return BSmoothSolution(0, res, exponentMatrix);
                        }
                    }
                    x += divisor;
                    indexX += divisor;
                }
            }
        }
    }

    return BSmoothSolution(0, res, exponentMatrix(Eigen::seqN(0, res.size()), Eigen::all));
}
