#include "bsmooth_square_finder.hpp"

BSmoothSolution::BSmoothSolution(
    const bmp::mpz_int &nDivisor,
    const std::vector<bmp::mpz_int> &bSmoothSquares,
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &exponentMatrix,
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &exponentMod2Matrix
) : nDivisor(nDivisor),
    bSmoothSquares(bSmoothSquares),
    exponentMatrix(exponentMatrix),
    exponentMod2Matrix(exponentMod2Matrix) {}

BSmoothSquareFinder::BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A)
    : n(n),
      s(static_cast<bmp::mpz_int>(bmp::sqrt(n)) + 1),
      factorBase(factorBase),
      piB(factorBase.size()),
      A(A),
      bSmoothBound(nextProbablePrime(1)),
      sieve(A + 1),
      primeExponentsLeqA(piB),
      qrsrModP(piB),
      primePowers(piB),
      currentRow(0) {
    // initialize sieve with ln f(x)
    // sieve will look like [0, 1, ..., A]
    // also initialize primeExponentsLeqA, primeExponentsLeqA[i] is the largest n s.t. factorBase[i]^n <= A
    // add negatives?
    for (std::size_t x = 0; x < A + 1; x++) {
        bmp::mpz_int base = s + x;
        sieve[x] = base * base - n;
    }

    for (std::size_t i = 0; i < piB; i++) {
        primeExponentsLeqA[i] = static_cast<std::size_t>(std::log(A) / std::log(factorBase[i]));
        primePowers[i].resize(primeExponentsLeqA[i]);
    }
}

bool BSmoothSquareFinder::findSolution(
    const QuadraticResidueSquareRoot &qrsr,
    std::size_t divisor,
    std::size_t numBSmoothSquares,
    std::size_t i,
    std::size_t k,
    std::vector<bmp::mpz_int> &res
) {
    for (std::size_t l = 0; l < qrsr.roots.size(); l++) {
        bmp::mpz_int x = qrsr.roots[l] - s;
        if (x < 0)
            x = (x % divisor + divisor) % divisor;
        bmp::mpz_int startX = x;
        std::ptrdiff_t indexX = static_cast<std::ptrdiff_t>(x);
        std::ptrdiff_t startIndexX = indexX;
        for (std::size_t j = 0, inc = factorBase[i]; j < k; j++, inc *= factorBase[i]) {
            // forwards
            while (indexX < A + 1) {
                sieve[indexX] /= factorBase[i];
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= bSmoothBound) {
                    // sieve[x] is approx. 0, B-smooth number found
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.row(currentRow) = primeFactorization.exponents;
                    exponentMod2Matrix.row(currentRow++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x += inc;
                indexX += inc;
            }
            // backwards
            x = startX - inc;
            indexX = startIndexX - inc;
            while (indexX >= 0) {
                sieve[indexX] /= factorBase[i];
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= bSmoothBound) {
                    // sieve[x] is approx. 0, B-smooth number found
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.row(currentRow) = primeFactorization.exponents;
                    exponentMod2Matrix.row(currentRow++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x -= inc;
                indexX -= inc;
            }
            x = startX;
            indexX = startIndexX;
        }
    }
    return false;
}

BSmoothSolution BSmoothSquareFinder::find(std::size_t numBSmoothSquares) {
    exponentMatrix.resize(numBSmoothSquares, piB);
    exponentMod2Matrix.resize(numBSmoothSquares, piB);

    std::vector<bmp::mpz_int> res;
    for (std::size_t i = 0; i < piB; i++) {
        // check if p | n by chance
        if (n % factorBase[i] == 0)
            return BSmoothSolution(factorBase[i], res, exponentMatrix, exponentMod2Matrix);

        if (factorBase[i] == 2) {
            bmp::mpz_int x = 0;
            bmp::mpz_int startX = x;
            std::size_t indexX = static_cast<std::size_t>(x);
            std::size_t startIndexX = indexX;
            if (sieve[indexX] % 2 != 0) {
                x++;
                startX++;
                indexX++;
                startIndexX++;
            }
            while (indexX < A + 1) {
                while (sieve[indexX] % 2 == 0) {
                    sieve[indexX] /= 2;
                }
                if (sieve[indexX] <= bSmoothBound) {
                    // sieve[x] is approx. 0, B-smooth number found
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.row(currentRow) = primeFactorization.exponents;
                    exponentMod2Matrix.row(currentRow++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return BSmoothSolution(0, res, exponentMatrix, exponentMod2Matrix); 
                }
                x += 2;
                indexX += 2;
            }
            continue;
        }

        qrsrModP[i] = factorBase[i] > 80 ? tonelliShanks(n, factorBase[i]) : trialSquareRoot(n, factorBase[i]);
        if (!qrsrModP[i].exists)
            continue;

        for (std::size_t k = primeExponentsLeqA[i]; k >= 1; k--) {
            std::size_t divisor = primePowers[i][k - 1] == 0 ? std::pow(factorBase[i], k) : primePowers[i][k - 1];

            // https://mathoverflow.net/questions/52081/is-there-an-efficient-algorithm-for-finding-a-square-root-modulo-a-prime-power
            // https://math.stackexchange.com/questions/94842/maximum-number-of-square-roots-of-a-in-mathbbz-n
            // TODO: write the case for p=2
            QuadraticResidueSquareRoot roots;
            if (divisor < 80) {
                roots = trialSquareRoot(n, divisor);
            } else {
                primePowers[i][k - 2] = divisor / factorBase[i];
                bmp::mpz_int ae = bmp::powm(n, (divisor - 2 * primePowers[i][k - 2] + 1) / 2, divisor);
                bmp::mpz_int root1r = bmp::powm(qrsrModP[i].roots[0], primePowers[i][k - 2], divisor);
                bmp::mpz_int root2r = bmp::powm(qrsrModP[i].roots[1], primePowers[i][k - 2], divisor);
                roots.roots = std::vector<bmp::mpz_int>{
                    root1r * ae % divisor,
                    root2r * ae % divisor
                };
                roots.exists = true;
            }
            if (!roots.exists)
                continue;

            if (findSolution(roots, divisor, numBSmoothSquares, i, k, res)) {
                return BSmoothSolution(0, res, exponentMatrix, exponentMod2Matrix);
            }
            break;
        }
    }

    return BSmoothSolution(0, res, exponentMatrix(Eigen::seqN(0, res.size()), Eigen::all), exponentMod2Matrix(Eigen::seqN(0, res.size()), Eigen::all));
}
