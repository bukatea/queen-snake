#include "tonelli_shanks.hpp"

QuadraticResidue tonelliShanks(std::uint64_t n, std::uint64_t p) {
    std::uint64_t q = p - 1;
    std::uint64_t ss = 0;
    std::uint64_t z = 2;
    std::uint64_t c, r, t, m;
 
    if (powerMod(n, (p - 1) / 2, p) != 1) {
        return QuadraticResidue(0, 0, false);
    }
 
    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }
 
    if (ss == 1) {
        std::uint64_t r1 = powerMod(n, (p + 1) / 4, p);
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
        std::uint64_t i = 0, zz = t;
        std::uint64_t b = c, e;
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