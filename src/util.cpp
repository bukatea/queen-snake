#include "util.hpp"

std::vector<std::size_t> factorBase(std::size_t B) {
    std::vector<bool> isPrime(B + 1, true);
    for (int i = 2; i * i <= B; i++) {
        if (isPrime[i]) {
            for (int j = i * i; j <= B; j += i)
                isPrime[j] = false;
        }
    }
    std::vector<std::size_t> primes;
    for (int i = 2; i <= B; i++) {
        if (isPrime[i])
            primes.push_back(i);
    }
    return primes;
}

bmp::mpz_int nextProbablePrime(std::size_t n) {
    if (n == 1)
        return 2;
    for (bmp::mpz_int i = n + (n % 2 == 0 ? 1 : 2); ; i += 2) {
        if (miller_rabin_test(i, 25)) {
            return i;
        }
    }
}

PrimeFactorization::PrimeFactorization(const Eigen::VectorXi &exponents, const Eigen::VectorXi &exponentsMod2) : exponents(exponents), exponentsMod2(exponentsMod2) {}

PrimeFactorization trialDivision(const std::vector<std::size_t> &factorBase, const bmp::mpz_int &n) {
    bmp::mpz_int num(n);
    Eigen::VectorXi exponents = Eigen::VectorXi::Zero(factorBase.size());
    Eigen::VectorXi exponentsMod2 = Eigen::VectorXi::Zero(factorBase.size());
    for (std::size_t i = 0; i < factorBase.size(); i++) {
        while (num % factorBase[i] == 0) {
            exponents[i]++;
            ++exponentsMod2[i] %= 2;
            num /= factorBase[i];
        }
    }
    return PrimeFactorization(exponents, exponentsMod2);
}
