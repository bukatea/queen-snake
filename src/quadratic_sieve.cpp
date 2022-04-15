#include "quadratic_sieve.hpp"

std::vector<bmp::mpz_int> quadraticSieve(const bmp::mpz_int &n, std::size_t A, std::size_t B) {
    std::vector<std::size_t> fB = factorBase(B);
    std::size_t piB = fB.size();

    std::size_t required = piB + 1;
    BSmoothSquareFinder finder(n, fB, A);
    BSmoothSolution solution = finder.find(required);
    if (solution.nDivisor != 0) {
        return std::vector<bmp::mpz_int>{solution.nDivisor, n / solution.nDivisor};
    }
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> linDependencies = mod2Echelonize(solution.exponentMod2Matrix);
    for (std::size_t i = 0; i < linDependencies.rows(); i++) {
        bmp::mpz_int x = 1;
        std::size_t j = 0;
        Eigen::VectorXi newExponents = Eigen::VectorXi::Zero(piB);
        for (auto val : linDependencies.row(i)) {
            if (val != 0) {
                newExponents += solution.exponentMatrix.row(j);
                x *= solution.bSmoothSquares[j];
                x %= n;
            }
            j++;
        }

        newExponents /= 2;
        bmp::mpz_int y = 1;
        j = 0;
        for (auto val : newExponents) {
            y *= static_cast<std::size_t>(std::pow(fB[j], val));
            y %= n;
            j++;
        }
        std::cout << "x: " << x << " y: " << y << std::endl;

        bmp::mpz_int factor = bmp::gcd(x + y, n);

        if (factor == n || factor == 1) {
            continue;
        }

        return std::vector<bmp::mpz_int>{factor, n / factor};
    }
    return std::vector<bmp::mpz_int>{};
}
