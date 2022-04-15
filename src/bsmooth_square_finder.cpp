#include "bsmooth_square_finder.hpp"

BSmoothSolution::BSmoothSolution(
    const bmp::mpz_int &nDivisor,
    const std::vector<bmp::mpz_int> &bSmoothSquares,
    const Eigen::MatrixXi &exponentMatrix,
    const Eigen::MatrixXi &exponentMod2Matrix
) : nDivisor(nDivisor),
    bSmoothSquares(bSmoothSquares),
    exponentMatrix(exponentMatrix),
    exponentMod2Matrix(exponentMod2Matrix) {}

BSmoothSquareFinder::BSmoothSquareFinder(const bmp::mpz_int &n, const std::vector<std::size_t> &factorBase, std::size_t A)
    : n(n),
      s(static_cast<bmp::mpz_int>(bmp::sqrt(n)) + 1),
      factorBase(factorBase),
      piB(factorBase.size()),
      // iterations(0),
      A(A),
      sieve(A + 1),
      primeExponentsLeqA(piB),
      qrsrModP(piB),
      primePowers(piB),
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
        primePowers[i].resize(primeExponentsLeqA[i]);
        // iterations += primeExponentsLeqA[i];
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
    // const static bmp::mpf_float SIEVE_ERROR = iterations * bmp::mpf_float{"0.000000000000000000000000000000000000000000000000005"};
    const static bmp::mpf_float_50 BSMOOTH_BOUND = bmp::log(static_cast<bmp::mpf_float_50>(factorBase[piB - 1]));

    bmp::mpf_float_50 logp = bmp::log(static_cast<bmp::mpf_float_50>(factorBase[i]));

    for (std::size_t l = 0; l < qrsr.roots.size(); l++) {
        std::cout << s << std::endl;
        bmp::mpz_int x = qrsr.roots[l] - s;
        if (x < 0)
            x = (x % divisor + divisor) % divisor;
        bmp::mpz_int startX = x;
        std::cout << "x " << x << std::endl;
        std::ptrdiff_t indexX = static_cast<std::ptrdiff_t>(x);
        std::ptrdiff_t startIndexX = indexX;
        std::cout << "indexX " << indexX << std::endl;
        for (std::size_t j = 0, inc = factorBase[i]; j < k; j++, inc *= factorBase[i]) {
            std::cout << "inc: " << inc << std::endl;
            std::cout << "forwards" << std::endl;
            // forwards
            while (indexX < A + 1) {
                sieve[indexX] -= logp;
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= BSMOOTH_BOUND) {
                    // sieve[x] is approx. 0, B-smooth number found
                    std::cout << "askjdfhlqksdhgkjadshflkadshfjadsx: " << x << std::endl;
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.col(currentCol) = primeFactorization.exponents;
                    exponentMod2Matrix.col(currentCol++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x += inc;
                indexX += inc;
            }
            // backwards
            std::cout << "backwards" << std::endl;
            x = startX - inc;
            indexX = startIndexX - inc;
            while (indexX >= 0) {
                sieve[indexX] -= logp;
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= BSMOOTH_BOUND) {
                    // sieve[x] is approx. 0, B-smooth number found
                    std::cout << "askjdfhlqksdhgkjadshflkadshfjadsx: " << x << std::endl;
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.col(currentCol) = primeFactorization.exponents;
                    exponentMod2Matrix.col(currentCol++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x -= inc;
                indexX -= inc;
            }
            x = startX;
            indexX = startIndexX;
        }
        if (factorBase[i] == 2) {
            x = qrsr.roots[l + 1] - s;
            if (x < 0)
                x = (x % divisor + divisor) % divisor;
            startX = x;
            std::cout << "x for 2 " << x << std::endl;
            indexX = static_cast<std::ptrdiff_t>(x);
            startIndexX = indexX;
            std::cout << "indexX " << indexX << std::endl;
            std::cout << "inc for 2: " << divisor << std::endl;
            std::cout << "forwards" << std::endl;
            // forwards
            while (indexX < A + 1) {
                sieve[indexX] -= logp;
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= BSMOOTH_BOUND) {
                    // sieve[x] is approx. 0, B-smooth number found
                    std::cout << "askjdfhlqksdhgkjadshflkadshfjadsx: " << x << std::endl;
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.col(currentCol) = primeFactorization.exponents;
                    exponentMod2Matrix.col(currentCol++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x += divisor;
                indexX += divisor;
            }
            // backwards
            std::cout << "backwards" << std::endl;
            x = startX - divisor;
            indexX = startIndexX - divisor;
            while (indexX >= 0) {
                sieve[indexX] -= logp;
                // if (-SIEVE_ERROR <= sieve[indexX] && sieve[indexX] <= SIEVE_ERROR) {
                if (sieve[indexX] <= BSMOOTH_BOUND) {
                    // sieve[x] is approx. 0, B-smooth number found
                    std::cout << "askjdfhlqksdhgkjadshflkadshfjadsx: " << x << std::endl;
                    bmp::mpz_int y = x + s;
                    res.emplace_back(y);
                    PrimeFactorization primeFactorization = trialDivision(factorBase, y * y - n);
                    exponentMatrix.col(currentCol) = primeFactorization.exponents;
                    exponentMod2Matrix.col(currentCol++) = primeFactorization.exponentsMod2;
                    if (res.size() == numBSmoothSquares)
                        return true;
                }
                x -= divisor;
                indexX -= divisor;
            }
            break;
        }
    }
    return false;
}

BSmoothSolution BSmoothSquareFinder::find(std::size_t numBSmoothSquares) {
    exponentMatrix.resize(numBSmoothSquares, piB);
    exponentMod2Matrix.resize(numBSmoothSquares, piB);

    std::vector<bmp::mpz_int> res;
    for (std::size_t i = 0; i < piB; i++) {
        std::cout << "p: " << factorBase[i] << std::endl;
        // check if p | n by chance
        if (n % factorBase[i] == 0)
            return BSmoothSolution(factorBase[i], res, exponentMatrix, exponentMod2Matrix);

        qrsrModP[i] = factorBase[i] > 80 ? tonelliShanks(n, factorBase[i]) : trialSquareRoot(n, factorBase[i]);
        if (!qrsrModP[i].exists)
            continue;
        std::cout << "roots for mod p are " << std::endl;
        for (const bmp::mpz_int &root : qrsrModP[i].roots) {
            std::cout << root << " ";
        }
        std::cout << std::endl;

        for (std::size_t k = primeExponentsLeqA[i]; k >= 1; k--) {
            std::cout << "k: " << k << std::endl;
            std::size_t divisor = primePowers[i][k - 1] == 0 ? std::pow(factorBase[i], k) : primePowers[i][k - 1];
            std::cout << "divisor: " << divisor << std::endl;

            std::cout << "finding roots" << std::endl;

            // https://mathoverflow.net/questions/52081/is-there-an-efficient-algorithm-for-finding-a-square-root-modulo-a-prime-power
            // https://math.stackexchange.com/questions/94842/maximum-number-of-square-roots-of-a-in-mathbbz-n
            // TODO: write the case for p=2
            QuadraticResidueSquareRoot roots;
            if (divisor < 80) {
                roots = trialSquareRoot(n, divisor);
            } else {
                primePowers[i][k - 2] = divisor / factorBase[i];
                std::cout << "asdf: " << primePowers[i][k - 2] << std::endl;
                if (factorBase[i] == 2) {
                    if (n % 8 == 1) {
                        roots.roots = std::vector<bmp::mpz_int>{
                            1,
                            1 + primePowers[i][k - 2],
                            primePowers[i][k - 2] - 1,
                            divisor - 1
                        };
                        roots.exists = true;
                    } else {
                        k = 3;
                    }
                } else {
                    bmp::mpz_int ae = bmp::powm(n, (divisor - 2 * primePowers[i][k - 2] + 1) / 2, divisor);
                    bmp::mpz_int root1r = bmp::powm(qrsrModP[i].roots[0], primePowers[i][k - 2], divisor);
                    bmp::mpz_int root2r = bmp::powm(qrsrModP[i].roots[1], primePowers[i][k - 2], divisor);
                    roots.roots = std::vector<bmp::mpz_int>{
                        root1r * ae % divisor,
                        root2r * ae % divisor
                    };
                    roots.exists = true;
                }
            }
            std::cout << "exists " << roots.exists << std::endl;
            std::cout << "roots found ";
            for (const bmp::mpz_int &root : roots.roots) {
                std::cout << root << " ";
            }
            std::cout << std::endl;
            if (!roots.exists)
                continue;

            if (findSolution(roots, divisor, numBSmoothSquares, i, k, res)) {
                return BSmoothSolution(0, res, exponentMatrix, exponentMod2Matrix);
            }
            break;
        }
    }

    return BSmoothSolution(0, res, exponentMatrix(Eigen::all, Eigen::seqN(0, res.size())), exponentMod2Matrix(Eigen::all, Eigen::seqN(0, res.size())));
}
