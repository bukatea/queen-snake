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
    if (solution.bSmoothSquares.size() < required) {
        std::cout << "fuck" << std::endl;
        return std::vector<bmp::mpz_int>{};
    } else {
        Eigen::MatrixXf ker = solution.exponentMatrix.cast<float>().fullPivLu().kernel();
        std::cout << ker << std::endl;
        Eigen::VectorXf firstSolution = ker.col(0);

        bmp::mpz_int x = 1;
        std::size_t i = 0;
        Eigen::VectorXi newExponents = Eigen::VectorXi::Zero(piB);
        for (auto val : firstSolution) {
            if (val != 0) {
                newExponents += Eigen::Map<Eigen::VectorXi>(&(solution.bSmoothSquares[i].primeFactors[0]), piB);
                x *= solution.bSmoothSquares[i].x;
            }
            i++;
        }

        newExponents /= 2;
        bmp::mpz_int y = 1;
        i = 0;
        for (auto val : newExponents) {
            y *= static_cast<std::size_t>(std::pow(fB[i], val));
            i++;
        }

        if (x % n == y % n || x % n == (-y % n + n) % n) {
            std::cout << "bitch" << std::endl;
            return std::vector<bmp::mpz_int>{};
        }

        return std::vector<bmp::mpz_int>{bmp::gcd(x + y, n), bmp::gcd(x > y ? x - y : y - x, n)};
    }
}
