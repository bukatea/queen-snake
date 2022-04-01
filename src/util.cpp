#include "util.hpp"

mpz_class powerMod(const mpz_class &base, const mpz_class &exponent, const mpz_class &modulus) {
    mpz_class x = 1;
    mpz_class y = a;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            x = (x * y) % modulus;
        }
        y = (y * y) % modulus;
        exponent /= 2;
    }
    return x % modulus;
}

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
