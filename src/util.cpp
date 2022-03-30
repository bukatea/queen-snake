#include "util.hpp"

std::uint64_t powerMod(std::uint64_t base, std::uint64_t exponent, std::uint64_t modulus) {
    std::uint64_t x = 1;
    std::uint64_t y = a;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            x = (x * y) % modulus; // multiplying with base
        }
        y = (y * y) % modulus; // squaring the base
        exponent /= 2;
    }
    return x % modulus;
}
