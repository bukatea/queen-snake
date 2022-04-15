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
    if (num > factorBase[factorBase.size() - 1]) {
        std::cout << "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Etiam et nunc rhoncus, fermentum quam eget, luctus mauris. Vestibulum at turpis facilisis, dignissim felis sed, fringilla augue. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Donec euismod volutpat pulvinar. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc imperdiet nisi turpis, vitae blandit tellus porta vitae. Nunc tincidunt dapibus sagittis. Nunc posuere faucibus ultricies. Pellentesque non lorem eget arcu consequat ultrices. Duis iaculis mi ac commodo egestas. Nulla facilisi. Mauris pretium tellus ac odio faucibus, in finibus felis consequat. Nulla facilisi. Nulla congue varius tortor vitae faucibus. Integer felis nulla, convallis a tincidunt fermentum, ornare in velit. Mauris efficitur nunc vel velit pellentesque posuere. Sed imperdiet lectus nec ante euismod, ac pulvinar libero tempus. In ut congue neque. Aliquam aliquet, sem sit amet finibus maximus, velit dolor scelerisque ex, eget dignissim massa ex sed orci. Vivamus eget lorem id sem malesuada scelerisque. Sed quis magna a ante congue rutrum eget eget tortor. Nunc semper elit pharetra tempor pulvinar. Duis vel felis bibendum, sagittis dolor in, mollis lectus. Cras eget porta libero. Fusce nec suscipit ligula. Fusce pellentesque augue nibh, a tempor enim imperdiet non. Suspendisse ac mollis ligula. Curabitur quis dolor eget magna dignissim maximus. Cras tempus ultricies venenatis. Aenean eu imperdiet enim. Morbi consequat egestas leo sit amet gravida. Duis feugiat sit amet velit in posuere. Pellentesque sed justo vel mi tincidunt consectetur quis sed lacus. Maecenas condimentum massa augue, eu sodales dolor finibus accumsan. Morbi auctor arcu id nulla tempor, vel porttitor est efficitur. Nam blandit, lectus at rhoncus posuere, odio enim porta lorem, at ultricies nibh felis ut metus. Duis elementum leo vel ante mollis, nec volutpat est ultricies. Interdum et malesuada fames ac ante ipsum primis in faucibus. In rhoncus lectus vel tempor accumsan. Fusce velit diam, scelerisque auctor semper sit amet, pulvinar eget metus. In semper nibh sed nisl dignissim, nec malesuada mauris auctor. Pellentesque quis iaculis nisl. Morbi ultricies non nulla id gravida. Donec at massa id erat pharetra rhoncus. Vivamus efficitur consequat tellus vel mattis. Quisque consequat eleifend gravida. Nunc sit amet lectus eu mi interdum aliquam et sed arcu. Fusce euismod dignissim sapien, ac venenatis erat pretium a. Proin et tellus et ligula tincidunt maximus at quis nunc." << std::endl;
    }
    return PrimeFactorization(exponents, exponentsMod2);
}
