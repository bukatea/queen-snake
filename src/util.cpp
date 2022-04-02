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
