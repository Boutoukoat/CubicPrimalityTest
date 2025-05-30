// -----------------------------------------------------------------------
// Collatz step function calculator
// -----------------------------------------------------------------------

#include "cubic_primality.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

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
    // detects more than 99% of non-perfect cubes.
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

    // compute the rounded-down integer cube root, check if the solution is exact.
    mpz_t u;
    mpz_init(u);
    uint64_t e = mpz_root(u, n, 3); // u = cuberoot(n), e non zero if computation is exact.
    mpz_clear(u);
    return e ? true : false;
}

// detect a perfect cube after efficient sieving
static bool uint64_perfect_cube(uint64_t a)
{
    // detects more than 99% of non-perfect cubes.
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

    // binary search (provides at least 1 more bit of accuracy per iteration)
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

    // divisibility is based on Barrett modular reductions by constants, use modular multiplications
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

    // no factor less than 151
    if (a <= 151 * 151)
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
        // step 1 :
        //   reduce by a multiple of small factors
        // step 2:
        //    divisibility is based on Barrett modular reductions by constants, use modular multiplications

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

        // 2^44-1 is divisible by 3, 5, 23, 89, .....
        mpz_init(ignore);
        a = mpz_mod_ui(ignore, n, (1ull << 44) - 1);
        mpz_clear(ignore);
        if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
            return COMPOSITE_FOR_SURE; // divisible by 23
        if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
            return COMPOSITE_FOR_SURE; // divisible by 89
    }

    // no trivial small factor detected
    return UNDECIDED; // might be prime
}

// Iterate a third order linear recurrence using "double and add" steps
// Computes (s,t,u)^e mod n
// Assume input s,t,u,a < n < 2^61
// Make output s,t,u < n < 2^61
static void uint64_exponentiate(uint64_t &s, uint64_t &t, uint64_t &u, uint64_t e, uint64_t n, uint64_t a)
{
    int bit = 63 - __builtin_clzll(e);
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;
    uint64_t tmp;

    while (bit--)
    {
        // Double
        tmp = mulmod(s, s, n);
        s2 = (uint128_t)tmp * a; // s2 < 2^122
        t2 = (uint128_t)t * t;   // t2 < 2^122
        u2 = (uint128_t)u * u;   // u2 < 2^122
        tmp = mulmod(s, t, n);
        st = (uint128_t)tmp * a; // st < 2^122
        tu = (uint128_t)t * u;   // tu < 2^122
        us = (uint128_t)u * s;   // us < 2^122
        st <<= 1;                // st < 2^123
        tu <<= 1;                // tu < 2^123
        us <<= 1;                // us < 2^123
        if ((e >> bit) & 1)
        {
            // add
            tmp = (us + t2 + s2) % n; // us + t2 + s2 < 2^125
            uu = (uint128_t)tmp * a;  // uu < 2^122
            ss = s2 + st + tu;        // ss < 2^125
            tt = u2 + st + uu;        // tt < 2^125
        }
        else
        {
            ss = s2 + us + t2; // ss < 2^125
            tt = s2 + st + tu; // tt < 2^125
            uu = u2 + st;      // uu < 2^123
        }
        s = ss % n; // s < n < 2^61
        t = tt % n; // t < n < 2^61
        u = uu % n; // u < n < 2^61
    }
}

// Iterate a third order linear recurrence using "double and add" steps
// Computes (s,t,u)^e mod n
// Make output s,t,u < n
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

static bool uint64_cubic_primality(uint64_t n, bool verbose = false)
{
    if (n >> 61)
    {
        // the cubic test might overflow for numbers > 61 bits along the additions
        // Better use another slower deterministic test for numbers <= 64 bits
        // More precise constraint is n < 2^64/6
        mpz_t v;
        mpz_init_set_ui(v, n);
        bool r = mpz_cubic_primality(v, verbose);
        mpz_clear(v);
        return r;
    }

    if ((n & 1) == 0)
    {
        return n == 2; // even
    }
    if (n == 1) return false;  // one is not prime 

    // sieve composites with small factors
    sieve_t sv = uint64_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        if (verbose)
            printf("Number has a small factor\n");
        return false; // composite
    case PRIME_FOR_SURE:
        if (verbose)
            printf("Small number has no small factor\n");
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    if (uint64_perfect_cube(n))
    {
        if (verbose)
            printf("Number is a perfect cube\n");
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
            if (verbose)
                printf("Number is a small prime\n");
            return true; // small prime
        }
        if (k >= 5405)
        {
            // This constraint seems to be unreachable
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
            if (verbose)
                printf("Number has a small factor\n");
            return false; //  composite
        }
        bs = 0;
        bt = 1;
        bu = 0;
        uint64_exponentiate(bs, bt, bu, n - 1, n, a);
        if (bs == 0 && bt == 0 && bu == 1)
        {
            // todo : composite for sure
            if (verbose)
                printf("Composite for sure ? more loops k %lu a %lu n %lu\n", k, a, n);
            continue; // B == 1 try another a
        }
        break;
    }

    // Now B = (0,1,0)^(n-1), compute B^2 as (bs2, bt2, bu2) = (bs, bt, bu)^2
    uint64_t bs2 = bs, bt2 = bt, bu2 = bu;
    uint64_exponentiate(bs2, bt2, bu2, 2, n, a);
    // check the final condition (B^2 + B + 1) == (-1, 1, a)
    bs = (bs + bs2) % n;
    bt = (bt + bt2) % n;
    bu = (bu + bu2 + 1) % n;
    if (bs != n - 1 || bt != 1 || bu != a)
    {
        if (verbose)
            printf("Number is composite\n");
        return false; // composite for sure
    }
    return true; // might be prime
}

bool mpz_cubic_primality(mpz_t n, bool verbose)
{
    if (verbose)
        gmp_printf("Testing 0x%Zx\n", n);
    if (mpz_cmp_ui(n, 1ull << 61) < 0)
    {
        // the cubic test will run in 64 bits calculations
        return uint64_cubic_primality(mpz_get_ui(n), verbose);
    }

    if (mpz_tstbit(n, 0) == 0)
    {
        if (verbose)
            printf("Number is even\n");
        return (mpz_cmp_ui(n, 2) == 0); // even
    }

    // detects small primes, small composites
    // detects smooth composites
    sieve_t sv = mpz_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        if (verbose)
            printf("Number has a small factor\n");
        return false; // composite
    case PRIME_FOR_SURE:
        if (verbose)
            printf("Small number has no small factor\n");
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    if (mpz_perfect_cube(n))
    {
        if (verbose)
            printf("Number is a perfect cube\n");
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
            if (verbose)
                printf("Number is a small prime\n");
            res = true; // small prime
            goto done2;
        }

        if (k >= 5405)
        {
            // This constraint seems to be unreachable
            // overflow, u > 2^64 next line
            printf("internal overflow, k>5404, a >= 29208627, u > 2^64\n");
            exit(1);
        }
        uint64_t u = (2 * k - 1) * a * (2 * a - 1);
        uint64_t v = mpz_mod_ui(ignore, n, u);
        uint64_t g = uint64_gcd(u, v);
        if (g > 1)
        {
            if (verbose)
                printf("Number has a small factor\n");
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
            if (verbose)
                gmp_printf("Composite ? more loops k %lu a %lu n %Zu\n", k, a, n);
            continue; // B == 1 try another a
        }
        break;
    }

    // Now B = (0,1,0)^(n-1), compute B^2 as (bs2, bt2, bu2) = (bs, bt, bu)^2
    mpz_t bs2, bt2, bu2;
    mpz_init_set(bs2, bs);
    mpz_init_set(bt2, bt);
    mpz_init_set(bu2, bu);
    mpz_set_ui(e, 2);

    mpz_exponentiate(bs2, bt2, bu2, e, n, a);
    // check the final condition (B^2 + B + 1) == (-1, 1, a)
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
                if (verbose)
                    printf("Number passed all test and is likely not a composite one\n");
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

// ------------------------------------------------------------------------------
// Simple foolguard unit tests
// ------------------------------------------------------------------------------

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

    // 2^127 - 1 (a prime)
    mpz_init_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 127);
    mpz_sub_ui(ma, ma, 1);
    assert(mpz_composite_sieve(ma) == UNDECIDED);
    // 2^127 + 1
    mpz_add_ui(ma, ma, 2);
    assert(mpz_composite_sieve(ma) == COMPOSITE_FOR_SURE);

    printf("Small primes (uint64)\n");
    assert(uint64_cubic_primality(16777259ull) == true);
    assert(uint64_cubic_primality(281474976710677ull) == true);

    printf("Small composites (uint64)\n");
    assert(uint64_cubic_primality(16777265ull) == false);
    assert(uint64_cubic_primality(281474976710683ull) == false);

    uint32_t smallq[] = {
        1,   1,   1,   3,   1,   5,  3,   3,  1,   9,  7,   5,   3,  17,  27,  3,  1,  29,  3,   21,  7,  17,  15,  9,
        43,  35,  15,  29,  3,   11, 3,   11, 15,  17, 25,  53,  31, 9,   7,   23, 15, 27,  15,  29,  7,  59,  15,  5,
        21,  69,  55,  21,  21,  5,  159, 3,  81,  9,  69,  131, 33, 15,  135, 29, 13, 131, 9,   3,   33, 29,  25,  11,
        15,  29,  37,  33,  15,  11, 7,   23, 13,  17, 9,   75,  3,  171, 27,  39, 7,  29,  133, 59,  25, 105, 129, 9,
        61,  105, 7,   255, 277, 81, 267, 81, 111, 39, 99,  39,  33, 147, 27,  51, 25, 281, 43,  71,  33, 29,  25,  9,
        451, 41,  277, 165, 67,  27, 7,   29, 51,  17, 169, 39,  67, 27,  27,  33, 85, 155, 87,  155, 37, 0};
    printf("Medium primes (mpz)\n");
    for (int j = 0; smallq[j]; j++)
    {
        // primes (2^j + q) must be catched
        mpz_init_set_ui(ma, 1);
        mpz_mul_2exp(ma, ma, j);
        mpz_add_ui(ma, ma, smallq[j]);
        assert(mpz_cubic_primality(ma) == true);
    }

    uint32_t smallp[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37,  41,  43,
                         47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 0};
    printf("Medium composites (mpz)\n");
    for (int j = 0; smallp[j]; j++)
    {
        // Composites p * (2^127-1) must be catched
        mpz_init_set_ui(ma, 1);
        mpz_mul_2exp(ma, ma, 127);
        mpz_sub_ui(ma, ma, 1);
        mpz_mul_ui(ma, ma, smallp[j]);
        assert(mpz_cubic_primality(ma) == false);
    }

    printf("Large primes (mpz)\n");

    // 11111...6442446...11111 (1001-digits) The smallest zeroless titanic palindromic prime
    // https://t5k.org/curios/page.php?number_id=3797
    char titanic[1002];
    memset(titanic, '1', 1001);
    titanic[497] = '6';
    titanic[498] = '4';
    titanic[499] = '4';
    titanic[500] = '2';
    titanic[501] = '4';
    titanic[502] = '4';
    titanic[503] = '6';
    titanic[1001] = 0;
    mpz_set_str(ma, titanic, 10);
    // gmp_printf("%Zd\n", ma);
    assert(mpz_cubic_primality(ma) == true);

    // 11111...0...3781 (1001-digits) The smallest Cyclops titanic prime.
    // https://t5k.org/curios/page.php?number_id=18967
    memset(titanic, '1', 1001);
    titanic[500] = '0';
    titanic[997] = '3';
    titanic[998] = '7';
    titanic[999] = '8';
    titanic[1000] = '1';
    titanic[1001] = 0;
    mpz_set_str(mb, titanic, 10);
    assert(mpz_cubic_primality(mb) == true);

    printf("Large composites (mpz)\n");
    // a semiprime out of the 2 previous tests, no small factors.
    mpz_mul(ma, mb, ma);
    assert(mpz_cubic_primality(ma) == false);

    // a large square
    mpz_mul(ma, mb, mb);
    assert(mpz_cubic_primality(ma) == false);

    // a large cube
    mpz_mul(ma, ma, mb);
    assert(mpz_cubic_primality(ma) == false);

    mpz_clears(ma, mb, 0);
}
