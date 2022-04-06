#include <iostream>
#include <vector>

#include "common.hpp"
#include "quadratic_sieve.hpp"

int main() {
    const std::vector<bmp::mpz_int> numbers = {
        "16921456439215439701",
        "46839566299936919234246726809",
        "6172835808641975203638304919691358469663",
        "3744843080529615909019181510330554205500926021947"
    };

    for (const bmp::mpz_int &number : numbers) {
        // TODO: initialize B
        std::size_t B;
        auto startTime = std::chrono::high_resolution_clock::now();
        std::vector<bmp::mpz_int> factors = quadraticSieve(number, B, true);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "The factors of " << number << " are" << std::endl;
        for (const bmp::mpz_int &factor : factors) {
            std::cout << "\t" << factor << std::endl;
        }
        std::chrono::duration<double> duration = endTime - startTime;
        std::cout << "It took " << duration.count() << "ms" << std::endl << std::endl;
    }

    return 0;
}
