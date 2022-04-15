#include "quadratic_residue.hpp"

QuadraticResidueSquareRoot::QuadraticResidueSquareRoot(const std::vector<bmp::mpz_int> &roots, bool exists) : roots(roots), exists(exists) {}

QuadraticResidueSquareRoot::QuadraticResidueSquareRoot(const std::vector<bmp::mpz_int> &&roots, bool exists) : roots(std::move(roots)), exists(exists) {}

QuadraticResidueSquareRoot trialSquareRoot(const bmp::mpz_int &n, std::size_t k) {
    bmp::mpz_int residue = n % k;
    std::vector<bmp::mpz_int> roots;
    for (std::size_t i = 0; i < k; i++) {
        if (i * i % k == residue) {
            roots.emplace_back(i);
        }
    }
    return QuadraticResidueSquareRoot(roots, !roots.empty());
}

QuadraticResidueSquareRoot tonelliShanks(const bmp::mpz_int &n, std::size_t p) {
    std::size_t q = p - 1;
    std::size_t ss = 0;
    bmp::mpz_int z = 2;
    bmp::mpz_int c, r, t;
    std::size_t m;
 
    if (bmp::powm(n, (p - 1) / 2, p) != 1) {
        return QuadraticResidueSquareRoot();
    }
 
    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }
 
    if (ss == 1) {
        bmp::mpz_int r1 = bmp::powm(n, (p + 1) / 4, p);
        return QuadraticResidueSquareRoot(std::vector<bmp::mpz_int>{r1, p - r1}, true);
    }
 
    while (bmp::powm(z, (p - 1) / 2, p) != p - 1) {
        z++;
    }
 
    c = bmp::powm(z, q, p);
    r = bmp::powm(n, (q + 1) / 2, p);
    t = bmp::powm(n, q, p);
    m = ss;
 
    while (true) {
        std::size_t i = 0;
        bmp::mpz_int zz = t;
        bmp::mpz_int b = c;
        if (t == 1) {
            return QuadraticResidueSquareRoot(std::vector<bmp::mpz_int>{r, p - r}, true);
        }
        while (zz != 1 && i < (m - 1)) {
            zz = zz * zz % p;
            i++;
        }
        std::size_t e = m - i - 1;
        while (e > 0) {
            b = b * b % p;
            e--;
        }
        r = r * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
    }
}
