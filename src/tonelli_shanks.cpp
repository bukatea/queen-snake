#include "tonelli_shanks.hpp"

QuadraticResidue::QuadraticResidue(const mpz_class &root1, const mpz_class &root2, bool exists) : root1(root1), root2(root2), exists(exists) {}

QuadraticResidue tonelliShanks(const mpz_class &n, const mpz_class &p) {
    mpz_class q = p - 1;
    mpz_class ss = 0;
    mpz_class z = 2;
    mpz_class c, r, t, m;
 
    if (powerMod(n, (p - 1) / 2, p) != 1) {
        return QuadraticResidue(0, 0, false);
    }
 
    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }
 
    if (ss == 1) {
        mpz_class r1 = powerMod(n, (p + 1) / 4, p);
        return QuadraticResidue(r1, p - r1, true);
    }
 
    while (powerMod(z, (p - 1) / 2, p) != p - 1) {
        z++;
    }
 
    c = powerMod(z, q, p);
    r = powerMod(n, (q + 1) / 2, p);
    t = powerMod(n, q, p);
    m = ss;
 
    while (true) {
        mpz_class i = 0, zz = t;
        mpz_class b = c, e;
        if (t == 1) {
            return QuadraticResidue(r, p - r, true);
        }
        while (zz != 1 && i < (m - 1)) {
            zz = zz * zz % p;
            i++;
        }
        e = m - i - 1;
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
