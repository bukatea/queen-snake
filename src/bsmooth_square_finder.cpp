#include "bsmooth_square_finder.hpp"

BSmoothSquare::BSmoothSquare(const bmp::mpz_int &x, const std::vector<std::size_t> &primeFactors) : x(x), primeFactors(primeFactors) {}

BSmoothSolution::BSmoothSolution(const bmp::mpz_int &nDivisor, const std::vector<BSmoothSquare> &bSmoothSquares) : nDivisor(nDivisor), bSmoothSquares(bSmoothSquares) {}

BSmoothSquareFinder::BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A)
    : n(n),
      s(static_cast<bmp::mpz_int>(std::sqrt(n)) + 1),
      factorBase(factorBase),
      piB(factorBase.size()),
      iterations(0),
      A(A),
      sieve(A + 1),
      primeExponentsLeqA(piB),
      primeFactors(A + 1, std::vector<std::size_t>(piB)) {
    // initialize sieve with ln f(x)
    // sieve will look like [0, 1, ..., A]
    // also initialize primeExponentsLeqA, primeExponentsLeqA[i] is the largest n s.t. factorBase[i]^n <= A
    // add negatives?
    for (std::size_t x = 0; x < A + 1; x++) {
        bmp::mpz_int base = s + x;
        sieve[x] = std::log(static_cast<bmp::mpf_float_50>(base * base - n));
    }

    for (std::size_t i = 0; i < piB; i++) {
        primeExponentsLeqA[i] = static_cast<std::size_t>(std::log(A) / std::log(factorBase[i]));
        iterations += primeExponentsLeqA[i];
    }
}

BSmoothSolution BSmoothSquareFinder::find(std::size_t numBSmoothSquares) {
    const static bmp::mpf_float SIEVE_ERROR = iterations * ERROR;

    std::vector<BSmoothSquare> res;
    for (std::size_t i = 0; i < piB; i++) {
        for (std::size_t k = primeExponentsLeqA[i]; k >= 1; k--) {
            std::size_t divisor = std::pow(factorBase[i], k);

            // check if p^k | n by chance
            if (n % divisor == 0)
                return BSmoothSolution(divisor, res);

            QuadraticResidueSquareRoot res = tonelliShanks(n, divisor);
            if (!res.exists) {
                continue;
            }

            bmp::mpf_float_50 logd = std::log(static_cast<bmp::mpf_float_50>(divisor));

            bmp::mpz_int x = res.root1 - s;
            if (x < 0)
                x = (x % divisor + divisor) % divisor;
            while (x < A + 1) {
                sieve[x] -= logd;
                primeFactors[x][i] += k % 2;
                if (-SIEVE_ERROR <= sieve[x] && sieve[x] <= SIEVE_ERROR) {
                    // sieve[x] is approx. 0, B-smooth number found
                    res.emplace_back(x + s, primeFactors[x]);
                    if (res.size() == numBSmoothSquares) {
                        return BSmoothSolution(0, res);
                    }
                }
                x += divisor;
            }
            if (factorBase[i] != 2) {
                // if the prime is 2 the roots should be the same, 1 (mod 2)
                x = res.root2 - s;
                if (x < 0)
                    x = (x % divisor + divisor) % divisor;
                while (x < A + 1) {
                    sieve[x] -= logd;
                    primeFactors[x][i] += k % 2;
                    if (-SIEVE_ERROR <= sieve[x] && sieve[x] <= SIEVE_ERROR) {
                        // sieve[x] is approx. 0, B-smooth number found
                        res.emplace_back(x + s, primeFactors[x]);
                        if (res.size() == numBSmoothSquares) {
                            return BSmoothSolution(0, res);
                        }
                    }
                    x += divisor;
                }
            }
        }
    }

    return res;
}
