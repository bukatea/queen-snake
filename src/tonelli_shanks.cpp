#include "tonelli_shanks.hpp"

QuadraticResidueSquareRoot::QuadraticResidueSquareRoot(const bmp::mpz_int &root1, const bmp::mpz_int &root2, bool exists) : root1(root1), root2(root2), exists(exists) {}

QuadraticResidueSquareRoot tonelliShanks(const bmp::mpz_int &n, std::size_t p) {
    std::size_t q = p - 1;
    std::size_t ss = 0;
    bmp::mpz_int z = 2;
    bmp::mpz_int c, r, t;
    std::size_t m;
 
    if (bmp::powm(n, (p - 1) / 2, p) != 1) {
        return QuadraticResidueSquareRoot(0, 0, false);
    }
 
    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }
 
    if (ss == 1) {
        bmp::mpz_int r1 = bmp::powm(n, (p + 1) / 4, p);
        return QuadraticResidueSquareRoot(r1, p - r1, true);
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
            return QuadraticResidueSquareRoot(r, p - r, true);
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
