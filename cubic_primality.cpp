// -----------------------------------------------------------------------
// Collatz step function calculator
// -----------------------------------------------------------------------

#include "cubic_primality.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

typedef unsigned __int128 uint128_t;

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t n)
{
#if 0
	uint128_t tmp;
	tmp = a;
	tmp *= b;
	tmp %= n;
	return tmp;
#else
    uint64_t r;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#endif
}

// binary gcd algorithm
static uint64_t uint64_gcd(uint64_t x, uint64_t y)
{
    if (x == 0)
        return y;
    if (y == 0)
        return x;
    unsigned tu = __builtin_ctzll(x);
    unsigned tv = __builtin_ctzll(y);
    unsigned h = tu > tv ? tv : tu;
    uint64_t t;
    uint64_t u = x >> tu;
    uint64_t v = y >> tv;

    while (1)
    {
        if (u > v)
        {
            t = u;
            u = v;
            v = t;
        }
        v -= u;
        if (v == 0)
        {
            return u << h;
        }
        v >>= __builtin_ctzll(v);
    }
}

// detect a perfect cube after efficient sieving
static bool mpz_perfect_cube(mpz_t n)
{
    mpz_t ignore;
    mpz_init(ignore);
    uint64_t a = mpz_mod_ui(ignore, n, 63ull * 54ull * 61ull * 43ull * 37ull * 31ull * 19ull * 13ull);
    mpz_clear(ignore);
    if (0x3f7fffe7e7fffefcull & (1ull << (a % 63)))
        return false;
    if (0x1fafd7e3f5fafcull & (1ull << (a % 54)))
        return false;
    if (0xbcbfd99e66ff4f4ull & (1ull << (a % 61)))
        return false;
    if (0x176f79ef6e8ull & (1ull << (a % 43)))
        return false;
    if (0xf537fb2bcull & (1ull << (a % 37)))
        return false;
    if (0x177e7ee8 & (1 << (a % 31)))
        return false;
    if (0x3e67c & (1 << (a % 19)))
        return false;
    if (0xedc & (1 << (a % 13)))
        return false;

    mpz_t u;
    mpz_init(u);
    uint64_t e = mpz_root(u, n, 3); // u = cuberoot(n), e non zero if computation is exact.
    mpz_clear(u);
    return e ? true : false;
}

// detect a perfect cube after efficient sieving
static bool uint64_perfect_cube(uint64_t a)
{
    if (0x3f7fffe7e7fffefcull & (1ull << (a % 63)))
        return false;
    if (0x1fafd7e3f5fafcull & (1ull << (a % 54)))
        return false;
    if (0xbcbfd99e66ff4f4ull & (1ull << (a % 61)))
        return false;
    if (0x176f79ef6e8ull & (1ull << (a % 43)))
        return false;
    if (0xf537fb2bcull & (1ull << (a % 37)))
        return false;
    if (0x177e7ee8 & (1 << (a % 31)))
        return false;
    if (0x3e67c & (1 << (a % 19)))
        return false;
    if (0xedc & (1 << (a % 13)))
        return false;

    // quick approximation of cubic root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 3.0); // cubic root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    uint64_t c, m;

    // binary search (1 more bit of cube root per iteration)
    uint64_t r = (uint64_t)d;
    uint64_t l = (uint64_t)dl;
    uint64_t h = (uint64_t)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m * m;
        if (c == a)
        {
            return true; // perfect cube
        }
        if (c < a)
        {
            l = m + 1;
            r = m;
        }
        else
        {
            h = m - 1;
        }
    }
    c = r * r * r; // check perfect cube
    return (c == a);
}

// modular exponentiation a^e mod m
static uint64_t uint64_powm(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t n = e;
    uint64_t b = a;
    uint64_t result = 1;

    while (n)
    {
        if (n & 1)
        {
            result = mulmod(b, result, m);
        }
        b = mulmod(b, b, m);
        n >>= 1;
    }
    return result;
}

typedef enum sieve_e
{
    COMPOSITE_FOR_SURE,
    PRIME_FOR_SURE,
    UNDECIDED
} sieve_t;

// sieve small primes
static sieve_t uint64_composite_sieve(uint64_t a)
{
    if (a <= 152)
    {
        // return COMPOSITE_FOR_SURE for small numbers which are composite for sure , without checking further.
        sieve_t stooopid_prime_table[] = {
            PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE};
        return stooopid_prime_table[a];
    }
    if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return COMPOSITE_FOR_SURE; // divisible by 3
    if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return COMPOSITE_FOR_SURE; // divisible by 5
    if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return COMPOSITE_FOR_SURE; // divisible by 7
    if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return COMPOSITE_FOR_SURE; // divisible by 11
    if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return COMPOSITE_FOR_SURE; // divisible by 13
    if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return COMPOSITE_FOR_SURE; // divisible by 17
    if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return COMPOSITE_FOR_SURE; // divisible by 19
    if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return COMPOSITE_FOR_SURE; // divisible by 23
    if ((uint64_t)(a * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return COMPOSITE_FOR_SURE; // divisible by 29
    if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return COMPOSITE_FOR_SURE; // divisible by 31
    if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return COMPOSITE_FOR_SURE; // divisible by 37
    if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return COMPOSITE_FOR_SURE; // divisible by 41
    if ((uint64_t)(a * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return COMPOSITE_FOR_SURE; // divisible by 43
    if ((uint64_t)(a * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return COMPOSITE_FOR_SURE; // divisible by 47
    if ((uint64_t)(a * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return COMPOSITE_FOR_SURE; // divisible by 53
    if ((uint64_t)(a * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return COMPOSITE_FOR_SURE; // divisible by 59
    if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return COMPOSITE_FOR_SURE; // divisible by 61
    if ((uint64_t)(a * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return COMPOSITE_FOR_SURE; // divisible by 67
    if ((uint64_t)(a * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return COMPOSITE_FOR_SURE; // divisible by 71
    if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return COMPOSITE_FOR_SURE; // divisible by 73
    if ((uint64_t)(a * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return COMPOSITE_FOR_SURE; // divisible by 79
    if ((uint64_t)(a * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return COMPOSITE_FOR_SURE; // divisible by 83
    if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return COMPOSITE_FOR_SURE; // divisible by 89
    if ((uint64_t)(a * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return COMPOSITE_FOR_SURE; // divisible by 97
    if ((uint64_t)(a * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return COMPOSITE_FOR_SURE; // divisible by 101
    if ((uint64_t)(a * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return COMPOSITE_FOR_SURE; // divisible by 103
    if ((uint64_t)(a * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return COMPOSITE_FOR_SURE; // divisible by 107
    if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return COMPOSITE_FOR_SURE; // divisible by 109
    if ((uint64_t)(a * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return COMPOSITE_FOR_SURE; // divisible by 113
    if ((uint64_t)(a * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return COMPOSITE_FOR_SURE; // divisible by 127
    if ((uint64_t)(a * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return COMPOSITE_FOR_SURE; // divisible by 131
    if ((uint64_t)(a * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return COMPOSITE_FOR_SURE; // divisible by 137
    if ((uint64_t)(a * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return COMPOSITE_FOR_SURE; // divisible by 139
    if ((uint64_t)(a * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return COMPOSITE_FOR_SURE; // divisible by 149
    if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return COMPOSITE_FOR_SURE; // divisible by 151
    if (a < 151 * 151)
        return PRIME_FOR_SURE; // prime
    return UNDECIDED;
}

static sieve_t mpz_composite_sieve(mpz_t n)
{
    // detect small 64-bit numbers
    unsigned sz = mpz_sizeinbase(n, 2);
    if (sz <= 64)
    {
        uint64_t a = mpz_get_ui(n);
        sieve_t sv = uint64_composite_sieve(a);
        return sv;
    }
    else
    {
        // large number, do the modular reduction in 2 steps
        uint64_t a;
        // 2^60-1 is divisible by 3,5,7,11,13,31,41,61,151 ...
        mpz_t ignore;
        mpz_init(ignore);
        a = mpz_mod_ui(ignore, n, (1ull << 60) - 1);
        mpz_clear(ignore);

        if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
            return COMPOSITE_FOR_SURE; // divisible by 3
        if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
            return COMPOSITE_FOR_SURE; // divisible by 5
        if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
            return COMPOSITE_FOR_SURE; // divisible by 7
        if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
            return COMPOSITE_FOR_SURE; // divisible by 11
        if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
            return COMPOSITE_FOR_SURE; // divisible by 13
        if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
            return COMPOSITE_FOR_SURE; // divisible by 31
        if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
            return COMPOSITE_FOR_SURE; // divisible by 41
        if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
            return COMPOSITE_FOR_SURE; // divisible by 61
        if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
            return COMPOSITE_FOR_SURE; // divisible by 151

        // 2^36-1 is divisible by 3,5,7,19,37,73,109, ...
        mpz_init(ignore);
        a = mpz_mod_ui(ignore, n, (1ull << 36) - 1);
        mpz_clear(ignore);
        if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
            return COMPOSITE_FOR_SURE; // divisible by 19
        if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
            return COMPOSITE_FOR_SURE; // divisible by 37
        if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
            return COMPOSITE_FOR_SURE; // divisible by 73
        if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
            return COMPOSITE_FOR_SURE; // divisible by 109

        // 2^48-1 is divisible by 13, 17, 97 .....
        mpz_init(ignore);
        a = mpz_mod_ui(ignore, n, (1ull << 48) - 1);
        mpz_clear(ignore);
        if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
            return COMPOSITE_FOR_SURE; // divisible by 17
        if ((uint64_t)(a * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
            return COMPOSITE_FOR_SURE; // divisible by 97
    }

    return UNDECIDED; // might be prime
}

static void uint64_exponentiate(uint64_t &s, uint64_t &t, uint64_t &u, uint64_t e, uint64_t n, uint64_t a)
{
    int bit = 63 - __builtin_clzll(e);
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;
    uint64_t tmp;

    while (bit--)
    {
        // Double
        tmp = mulmod(s, s, n);
        s2 = (uint128_t)tmp * a;
        t2 = (uint128_t)t * t;
        u2 = (uint128_t)u * u;
        tmp = mulmod(s, t, n);
        st = (uint128_t)tmp * a;
        tu = (uint128_t)t * u;
        us = (uint128_t)u * s;
        st <<= 1;
        tu <<= 1;
        us <<= 1;
        if ((e >> bit) & 1)
        {
            // add
            tmp = (us + t2 + s2) % n;
	    uu = (uint128_t)tmp * a;
            ss = s2 + st + tu;
            tt = u2 + st + uu;
        }
        else
        {
            ss = s2 + us + t2;
            tt = s2 + st + tu;
            uu = u2 + st;
        }
        s = ss % n;
        t = tt % n;
        u = uu % n;
    }
}

static void mpz_exponentiate(mpz_t s, mpz_t t, mpz_t u, mpz_t e, mpz_t n, uint64_t a)
{
    int bit = mpz_sizeinbase(e, 2) - 1;
    mpz_t tmp, s2, t2, u2, st, tu, us, uu, ss, tt;
    mpz_inits(tmp, s2, t2, u2, st, tu, us, uu, ss, tt, 0);

    while (bit--)
    {
        // Double
        mpz_mul(s2, s, s);
        mpz_mul_ui(s2, s2, a);
        mpz_mul(t2, t, t);
        mpz_mul(u2, u, u);
        mpz_mul(st, s, t);
        mpz_mul_ui(st, st, a);
        mpz_mul(tu, t, u);
        mpz_mul(us, u, s);
        mpz_mul_2exp(st, st, 1);
        mpz_mul_2exp(tu, tu, 1);
        mpz_mul_2exp(us, us, 1);
        if (mpz_tstbit(e, bit))
        {
            // add
            mpz_add(uu, us, t2);
            mpz_add(uu, uu, s2);
            mpz_mul_ui(uu, uu, a);
            mpz_add(ss, s2, st);
            mpz_add(ss, ss, tu);
            mpz_add(tt, u2, st);
            mpz_add(tt, tt, uu);
        }
        else
        {
            mpz_add(ss, s2, us);
            mpz_add(ss, ss, t2);
            mpz_add(tt, s2, st);
            mpz_add(tt, tt, tu);
            mpz_add(uu, u2, st);
        }
        mpz_mod(s, ss, n);
        mpz_mod(t, tt, n);
        mpz_mod(u, uu, n);
    }
    mpz_clears(tmp, s2, t2, u2, st, tu, us, uu, ss, tt, 0);
}

static bool uint64_cubic_primality(uint64_t n)
{
    if (n >> 61)
    {
        // the cubic test might overflow for numbers > 61 bits along the additions
        // Better use another slower deterministic test for numbers <= 64 bits
        // More precise constraint is n < 2^64/6
        mpz_t v;
        mpz_init_set_ui(v, n);
        bool r = mpz_cubic_primality(v);
        mpz_clear(v);
        return r;
    }

    if ((n & 1) == 0)
    {
        return n == 2; // even
    }
    // sieve composites with small factors
    sieve_t sv = uint64_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        return false; // composite
    case PRIME_FOR_SURE:
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    if (uint64_perfect_cube(n))
    {
        return false; // composite
    }

    uint64_t k = 0, a, bs, bt, bu;

    while (1)
    {
        k += 1;
        a = 7 + k * (k - 1);
        if (!uint64_cubic_primality(a))
        {
            continue; // try another a
        }
        if (uint64_powm(n % a, (a - 1) / 3, a) == 1)
        {
            continue; // try another a
        }
        if (a == n)
        {
            return true; // small prime
        }
        if (k >= 5405)
        {
            // overflow, u > 2^64 next line
            printf("internal overflow, k>5404, a >= 29208627, u > 2^64\n");
            exit(1);
        }
        uint64_t u = (2 * k - 1) * a * (2 * a - 1);
        uint64_t g = uint64_gcd(u, n);
        if (g == n)
        {
            continue; // try another a
        }
        if (g > 1)
        {
            return false; //  composite
        }
        bs = 0;
        bt = 1;
        bu = 0;
        uint64_exponentiate(bs, bt, bu, n - 1, n, a);
        if (bs == 0 && bt == 0 && bu == 1)
        {
		// todo : composite for sure
            // printf("more loops k %lu a %lu n %lu\n", k, a, n);
            continue; // B == 1 try another a
        }
        break;
    }
    uint64_t bs2 = bs, bt2 = bt, bu2 = bu;
    uint64_exponentiate(bs2, bt2, bu2, 2, n, a);
    bs = (bs + bs2) % n;
    bt = (bt + bt2) % n;
    bu = (bu + bu2 + 1) % n;
    if (bs != n - 1 || bt != 1 || bu != a)
        return false; // composite
    return true;      // might be prime
}

bool mpz_cubic_primality(mpz_t n)
{
    if (mpz_cmp_ui(n, 1ull << 61) < 0)
    {
        // the cubic test will run in 64 bits calculations
        return uint64_cubic_primality(mpz_get_ui(n));
    }

    if (mpz_tstbit(n, 0) == 0)
    {
        return (mpz_cmp_ui(n, 2) == 0); // even
    }
    // detects small primes, small composites
    // detects smooth composites
    sieve_t sv = mpz_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        return false; // composite
    case PRIME_FOR_SURE:
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    if (mpz_perfect_cube(n))
    {
        return false; // composite
    }

    uint64_t a, k = 0;
    bool res = false;
    mpz_t bs, bt, bu, ignore, e;
    mpz_inits(bs, bt, bu, ignore, e, 0);

    while (1)
    {
        k += 1;
        a = 7 + k * (k - 1);
        if (!uint64_cubic_primality(a))
        {
            continue; // try another a
        }
        if (uint64_powm(mpz_mod_ui(ignore, n, a), (a - 1) / 3, a) == 1)
        {
            continue; // try another a
        }
        if (mpz_cmp_ui(n, a) == 0)
        {
            res = true; // small prime
            goto done2;
        }

        if (k >= 5405)
        {
            // overflow, u > 2^64 next line
            printf("internal overflow, k>5404, a >= 29208627, u > 2^64\n");
            exit(1);
        }
        uint64_t u = (2 * k - 1) * a * (2 * a - 1);
        uint64_t v = mpz_mod_ui(ignore, n, u);
        uint64_t g = uint64_gcd(u, v);
        if (g > 1)
        {
            res = false; // composite
            goto done2;
        }
        mpz_set_ui(bs, 0);
        mpz_set_ui(bt, 1);
        mpz_set_ui(bu, 0);
        mpz_sub_ui(e, n, 1);
        mpz_exponentiate(bs, bt, bu, e, n, a);
        if (mpz_sgn(bs) == 0 && mpz_sgn(bt) == 0 && mpz_cmp_ui(bu, 1) == 0)
        {
		// todo : composite for sure
            // gmp_printf("more loops k %lu a %lu n %Zu\n", k, a, n);
            continue; // B == 1 try another a
        }
        break;
    }
    mpz_t bs2, bt2, bu2;
    mpz_init_set(bs2, bs);
    mpz_init_set(bt2, bt);
    mpz_init_set(bu2, bu);
    mpz_set_ui(e, 2);

    mpz_exponentiate(bs2, bt2, bu2, e, n, a);
    mpz_add(bs, bs, bs2);
    mpz_add_ui(bs, bs, 1);
    mpz_mod(bs, bs, n);
    if (mpz_sgn(bs) == 0)
    {
        mpz_add(bt, bt, bt2);
        mpz_mod(bt, bt, n);
        if (mpz_cmp_ui(bt, 1) == 0)
        {
            mpz_add(bu, bu, bu2);
            mpz_add_ui(bu, bu, 1);
            mpz_mod(bu, bu, n);
            if (mpz_cmp_ui(bu, a) == 0)
            {
                res = true;
                goto done1;
            }
        }
    }

done1:
    mpz_clears(bs2, bt2, bu2, 0);
done2:
    mpz_clears(bs, bt, bu, ignore, e, 0);
    return res;
}

void cubic_primality_self_test(void)
{
    uint64_t a, b, c, m;

    printf("Gcd ...\n");
    a = 12;
    b = 15;
    c = uint64_gcd(a, b);
    assert(c == 3);

    a = 120;
    b = 150;
    c = uint64_gcd(a, b);
    assert(c == 30);

    printf("Modexp ...\n");
    m = 101;
    a = 200;
    b = 300;
    c = uint64_powm(a, b, m);
    assert(c == 1);

    m = 101;
    a = 201;
    b = 301;
    c = uint64_powm(a, b, m);
    assert(c == 100);

    printf("Perfect cube ...\n");
    mpz_t ma, mb;
    bool mc;
    mpz_init_set_ui(mb, 1);
    mpz_init_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 11213);
    mpz_sub_ui(ma, ma, 1);
    mpz_mul(mb, mb, ma);
    mpz_mul(mb, mb, ma);
    mpz_mul(mb, mb, ma);
    mc = mpz_perfect_cube(mb);
    assert(mc == true);
    mpz_sub_ui(mb, mb, 1);
    mc = mpz_perfect_cube(mb);
    assert(mc == false);

    printf("Sieve ...\n");
    assert(uint64_composite_sieve(101) == PRIME_FOR_SURE);
    assert(uint64_composite_sieve(1661) == COMPOSITE_FOR_SURE);
    assert(uint64_composite_sieve(281474976710677ull) == UNDECIDED);

    mpz_init_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 127);
    mpz_sub_ui(ma, ma, 1);
    assert(mpz_composite_sieve(ma) == UNDECIDED);
    mpz_add_ui(ma, ma, 2);
    assert(mpz_composite_sieve(ma) == COMPOSITE_FOR_SURE);

    mpz_clears(ma, mb, 0);

    printf("Small primes\n");
    assert(uint64_cubic_primality(16777259ull) == true);
    assert(uint64_cubic_primality(16777265ull) == false);
    assert(uint64_cubic_primality(281474976710677ull) == true);
    assert(uint64_cubic_primality(281474976710683ull) == false);
}
