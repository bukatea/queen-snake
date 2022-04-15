#include <iostream>
#include <vector>

#include "common.hpp"
#include "gaussian_elimination.hpp"
#include "quadratic_sieve.hpp"

int main() {
    const std::vector<bmp::mpz_int> numbers{
        bmp::mpz_int{"16921456439215439701"},
        bmp::mpz_int{"46839566299936919234246726809"},
        bmp::mpz_int{"6172835808641975203638304919691358469663"},
        bmp::mpz_int{"3744843080529615909019181510330554205500926021947"}
    };

    for (const bmp::mpz_int &number : numbers) {
        // TODO: initialize A, B
        std::size_t A = static_cast<std::size_t>(bmp::ceil(bmp::pow(static_cast<bmp::mpf_float>(number), 0.25)));
        std::size_t B = static_cast<std::size_t>(bmp::ceil(bmp::exp((0.5 + 0.25) * bmp::sqrt((bmp::log(static_cast<bmp::mpf_float>(number)) * bmp::log(bmp::log(static_cast<bmp::mpf_float>(number))))))));
        // std::size_t B = static_cast<std::size_t>(bmp::ceil(bmp::exp((0.5 + 0.1) * bmp::sqrt((bmp::log(static_cast<bmp::mpf_float>(number)) * bmp::log(bmp::log(static_cast<bmp::mpf_float>(number))))))));
        // std::size_t A = 2 * std::pow(B, 3);
        std::cout << A << " " << B << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        std::vector<bmp::mpz_int> factors = quadraticSieve(number, A, B);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "The factors of " << number << " are" << std::endl;
        for (const bmp::mpz_int &factor : factors) {
            std::cout << "\t" << factor << std::endl;
        }
        std::chrono::duration<double> duration = endTime - startTime;
        std::cout << "It took " << duration.count() << "s" << std::endl << std::endl;
    }

    return 0;
}
