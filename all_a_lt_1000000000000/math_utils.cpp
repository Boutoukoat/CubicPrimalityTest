

#include <algorithm>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <vector>

#include <x86intrin.h>

using namespace std;

typedef unsigned __int128 uint128_t;
typedef signed __int128 int128_t;

// --------------------------------------------------------------------------------------
//
// C optimized implementations of simple utilities
//
// --------------------------------------------------------------------------------------

// randomness generation with period 2^64
uint64_t uint64_rnd(void)
{
    static uint64_t s = 0x1234567890123456ull;
    s = s * 137 + 13;
    return (s >> 13) ^ (s << 13) ^ s;
}

// (u + v) mod n
// Assume u+v < 2*n
static inline uint64_t uint64_add_mod(uint64_t u, uint64_t v, uint64_t n)
{
    uint128_t t = u + v;
    return t < n ? t : t - n;
}

// (u) mod n
// (u << 64 + v) mod n
static inline uint64_t uint64_long_mod(uint64_t u, uint64_t v, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(u), "1"(v), "r"(n) : "flags");
    return r;
#else
    uint128_t t = ((uint128_t)u << 64) + v;
    return t % n;
#endif
}

// (u) mod n
static inline uint64_t uint128_long_mod(uint128_t u, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r = (uint64_t)(u >> 64), a = (uint64_t)u;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    return u % n;
#endif
}

// u(hi, lo) mod n
static inline void uint128_divrem(uint64_t *q, uint64_t *r, uint128_t u, uint64_t n)
{
#ifdef __x86_64__
    uint64_t hi = (uint64_t)(u >> 64), lo = (uint64_t)u;
    asm("divq %4" : "=d"(*r), "=a"(*q) : "0"(hi), "1"(lo), "r"(n) : "flags");
#else
    *q = u / n;
    *r = u % n;
#endif
}

static inline uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    uint128_t tmp = (uint128_t)a * b;
    tmp %= n;
    return tmp;
#endif
}

static inline uint64_t square_mod(uint64_t a, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %2" : "=d"(r), "=a"(a) : "1"(a));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    uint128_t tmp = (uint128_t)a * a;
    tmp %= n;
    return tmp;
#endif
}

// (u * u + s) mod n
static inline uint64_t square_add_mod(uint64_t u, uint64_t s, uint64_t n)
{
#if 0
   // need to consider this in pollard-brent algorithm 
   // u = (u * u + s) / 2
   uint128_t t = (uint128_t)u * u;
    t += s;
    t += (t & 1) ? n : 0;
    t >>= 1;
    t -= (t >= n) ? n : 0;
    return t % n;
#endif

#ifdef __x86_64__
    uint64_t r, a;
    asm("mulq %2" : "=d"(r), "=a"(a) : "1"(u) : "flags");
    asm("addq %2, %1\n\tadcq $0, %0" : "+d"(r), "+a"(a) : "r"(s) : "flags");
    asm("divq %2" : "+d"(r), "+a"(a) : "r"(n) : "flags");
    return r;
#else
    uint128_t t = (uint128_t)u * u;
    t += s;
    return t % n;
#endif
}

// (u << s) mod n
static inline uint64_t shift_mod(uint64_t u, uint64_t s, uint64_t n)
{
#ifdef __x86_164__
    uint64_t r;
    asm("xor %0, %0\n shldq %b3, %1, %0\n shlxq %3, %1, %%rax\n divq %2"
        : "=&d"(r)
        : "r"(u), "r"(n), "c"(s)
        : "flags", "%rax");
    return r;
#else
    uint128_t t = (uint128_t)u;
    t <<= s;
    return t % n;
#endif
}

// count leading zeroed bits
static inline uint64_t uint64_lzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("lzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_clzll(a);
#endif
}

// count trailing zeroed bits
static inline uint64_t uint64_tzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("tzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_ctzll(a);
#endif
}

// count trailing zeroed bits
static inline uint64_t uint128_tzcnt(uint128_t a)
{
    uint64_t t = (uint64_t)a;
    if (t)
    {
        return uint64_tzcnt(t);
    }
    t = (a >> 64);
    return 64 + uint64_tzcnt(t);
}

// simple floor(log_2) function
// log(1) = 0
// log(2) = 1
// log(3) = 1
// log(4) = 2 ...
static inline uint64_t uint64_log_2(uint64_t a)
{
    return 63 - uint64_lzcnt(a);
}

static inline uint64_t uint128_log_2(uint128_t a)
{
    uint64_t t = (uint64_t)(a >> 64);
    if (t)
    {
        return 127 - uint64_lzcnt(t);
    }
    else
    {
        return 63 - uint64_lzcnt((uint64_t)a);
    }
}

// --------------------------------------------------------------------------------------
//
// Math functions
//
// --------------------------------------------------------------------------------------

static int uint64_jacobi(uint64_t x, uint64_t y)
{
    // assert((y & 1) == 1);
    if (y == 1 || x == 1)
    {
        return 1;
    }

    if (x == 2)
    {
        // char j[4] = { -1,-1,1,1};
        // return j[((y - 3) >> 1) % 4];
        return ((y + 2) & 4) ? -1 : 1;
    }
    if (x == 3)
    {
        char j[6] = {0, (char)-1, (char)-1, 0, 1, 1};
        return j[((y - 3) >> 1) % 6];
    }
    if (x == 5)
    {
        char j[5] = {(char)-1, 0, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 5];
    }
    if (x == 7)
    {
        char j[14] = {1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 0, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 14];
    }
    if (x == 11)
    {
        char j[22] = {(char)-1, 1,        1,        1, 0, (char)-1, (char)-1, (char)-1, 1, (char)-1, (char)-1, 1,
                      (char)-1, (char)-1, (char)-1, 0, 1, 1,        1,        (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 22];
    }
    if (x == 13)
    {
        char j[13] = {1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 13];
    }
    if (x == 17)
    {
        char j[17] = {(char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        1, 0, 1,
                      1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 17];
    }
    if (x == 19)
    {
        char j[19] = {1,        1, (char)-1, 1,        (char)-1, (char)-1, 1,        1,        0,       (char)-1,
                      (char)-1, 1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 38;
        return t >= 19 ? -j[t - 19] : j[t];
    }
    if (x == 23)
    {
        char j[23] = {(char)-1, (char)-1, 1,        1, 1,        1,        1,        (char)-1,
                      1,        (char)-1, 0,        1, (char)-1, 1,        (char)-1, (char)-1,
                      (char)-1, (char)-1, (char)-1, 1, 1,        (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 46;
        return t >= 23 ? -j[t - 23] : j[t];
    }
    if (x == 29)
    {
        char j[29] = {(char)-1, 1, 1,        1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1,
                      1,        1, (char)-1, 0, (char)-1, 1, 1,        (char)-1, (char)-1, (char)-1,
                      (char)-1, 1, (char)-1, 1, 1,        1, (char)-1, 1,        1};
        return j[((y - 3) >> 1) % 29];
    }
    if (x == 31)
    {
        char j[31] = {1,        1, (char)-1, 1,        1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                      1,        1, (char)-1, 0,        1, (char)-1, (char)-1, (char)-1, 1,        1,        1,
                      (char)-1, 1, (char)-1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 62;
        return t >= 31 ? -j[t - 31] : j[t];
    }

    int t = 1;
    uint64_t a = x;
    uint64_t n = y;
    unsigned v = n & 7;
    unsigned c = (v == 3) || (v == 5);
    while (a)
    {
        v = __builtin_ctzll(a);
        a >>= v;
        t = (c & (v & 1)) ? -t : t;

        if (a < n)
        {
            uint64_t tmp = a;
            a = n;
            n = tmp;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            c = (v == 3) || (v == 5);
        }

        a -= n;
    }

    return (n == 1) ? t : 0;
}

// Kronecker symbol (x, y) based on Stein's algorithm
static int int64_kronecker(int64_t x, int64_t y)
{
    unsigned x1 = x & 1;
    unsigned y1 = y & 1;
    if (y1 == 1)
    {
        // bit wizardry for K(+/- 2,y), K(+/- 3,y) when y odd

        if (x == 0)
        {
            return (y == 1) ? 1 : 0;
        }
        if (x == 1)
        {
            return 1;
        }
        if (x == -1)
        {
            return (y & 2) ? -1 : 1;
        }
        if (x == 2)
        {
            return ((y + 2) & 4) ? -1 : 1;
        }
        if (x == -2)
        {
            return (y & 4) ? -1 : 1;
        }
        if (x1 == 1)
        {
            if (x >= 0)
            {
                return uint64_jacobi(x, y);
            }
            else
            {
                int j = uint64_jacobi(-x, y);
                return (y & 2) ? -j : j;
            }
        }
    }
    else
    {
        // bit wizardry when y even
        if ((x1 == 0 && y1 == 0) || (y == 0 && !(x == 1 || x == -1)))
        {
            return 0;
        }
        if (y == 0)
        {
            return 1;
        }
    }

    // generic case
    int64_t a = x;
    int64_t n = y;
    int64_t a0 = a;
    unsigned v, cn, ca;
    int t = 1;

    // handle negative numbers
    if (a < 0)
    {
        if (n < 0)
        {
            n = -n;
            a = -a;
            t = -t;
        }
        else
        {
            a = -a;
        }
    }
    else if (n < 0)
    {
        n = -n;
    }

    // gulp trailing zeroes from n
    v = a & 7;
    ca = (v == 3) || (v == 5);
    v = uint64_tzcnt(n);
    n >>= v;
    t = (ca & (v & 1)) ? -t : t;

    if (a0 < 0 && (n & 3) == 3)
    {
        t = -t;
    }
    v = n & 7;
    cn = (v == 3) || (v == 5);

    // gulp trailing zeroes from a within Stein's algorithm.
    while (a > 0)
    {
        v = uint64_tzcnt(a);
        a >>= v;
        t = (cn & (v & 1)) ? -t : t;
        if (a == 1)
        {
            n = 1;
            break;
        }
        if (a < n)
        {
            int64_t r = a;
            a = n;
            n = r;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            cn = (v == 3) || (v == 5);
        }
        a = a - n;
    }
    return (n == 1) ? t : 0;
}

// integer square root (rounded down)
uint64_t uint64_isqrt(uint64_t x)
{
    // Avoid divide by zero
    if (x < 2)
    {
        return x;
    }
    // This code is based on the fact that
    // sqrt(x) == x^1/2 == 2^(log2(x)/2)
    // Unfortunately it's a little more tricky
    // when fast log2 is floored.
    uint64_t log2x = uint64_log_2(x);
    uint64_t log2y = log2x / 2;
    uint64_t y = 1 << log2y;
    uint64_t y_squared = 1 << (2 * log2y);
    int64_t sqr_diff = x - y_squared;
    // Perform lerp between powers of four
    y += (sqr_diff / 3) >> log2y;
    // The estimate is probably too low, refine it upward
    y_squared = y * y;
    sqr_diff = x - y_squared;
    y += sqr_diff / (2 * y);
    // The estimate may be too high. If so, refine it downward
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    y -= (-sqr_diff / (2 * y)) + 1;
    // The estimate may still be 1 too high
    y_squared = y * y;
    sqr_diff = x - y_squared;
    return sqr_diff < 0 ? y - 1 : y;
}

// return smallest factor of n < 157*157, or 1 if none is found.
uint64_t uint64_small_factor(uint64_t n)
{
    if (n <= 152)
    {
        // return smallest factor of n <= 151, or 1 if none is found.
        uint8_t stooopid_factor_table[] = {
            1, 1, 1, 1, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1,  2, 3, 2, 1, 2, 5, 2, 3, 2,  1, 2,
            1, 2, 3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 7, 2,  3, 2, 1, 2, 5, 2, 3, 2, 1,  2, 1,
            2, 3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 7, 2, 1, 2, 3,  2, 1, 2, 5, 2, 3, 2, 1, 2,  7, 2,
            3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2,  1, 2, 5, 2, 3, 2, 7, 2, 11, 2, 3,
            2, 5, 2, 1, 2, 3, 2, 1, 2, 7, 2, 3, 2, 1, 2, 1, 2, 3, 2, 11, 2, 5, 2, 3, 2, 1, 2, 1, 2};
        return stooopid_factor_table[n];
    }
    if (!(n & 1))
        return 2;
    if ((uint64_t)(n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return 3;
    if ((uint64_t)(n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return 5;
    if ((uint64_t)(n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return 7;
    if ((uint64_t)(n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return 11;
    if ((uint64_t)(n * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return 13;
    if ((uint64_t)(n * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return 17;
    if ((uint64_t)(n * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return 19;
    if ((uint64_t)(n * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return 23;
    if ((uint64_t)(n * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return 29;
    if ((uint64_t)(n * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return 31;
    if (n < 37 * 37)
        return 1; // prime
    if ((uint64_t)(n * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return 37;
    if ((uint64_t)(n * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return 41;
    if ((uint64_t)(n * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return 43;
    if ((uint64_t)(n * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return 47;
    if ((uint64_t)(n * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return 53;
    if ((uint64_t)(n * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return 59;
    if ((uint64_t)(n * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return 61;
    if ((uint64_t)(n * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return 67;
    if ((uint64_t)(n * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return 71;
    if ((uint64_t)(n * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return 73;
    if ((uint64_t)(n * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return 79;
    if ((uint64_t)(n * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return 83;
    if ((uint64_t)(n * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return 89;
    if ((uint64_t)(n * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return 97;
    if (n < 101 * 101)
        return 1; // prime
    if ((uint64_t)(n * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return 101;
    if ((uint64_t)(n * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return 103;
    if ((uint64_t)(n * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return 107;
    if ((uint64_t)(n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return 109;
    if ((uint64_t)(n * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return 113;
    if ((uint64_t)(n * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return 127;
    if ((uint64_t)(n * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return 131;
    if ((uint64_t)(n * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return 137;
    if ((uint64_t)(n * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return 139;
    if ((uint64_t)(n * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return 149;
    if ((uint64_t)(n * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return 151;
    return 1; // no small factor < 157
}

// modular exponentiation a^e mod m
uint64_t pow2_mod(uint64_t e, uint64_t m)
{
    uint64_t n = uint64_log_2(e);
    uint64_t s = (n >= 5) ? 5 : n;
    n -= s;
    uint64_t mask = e >> n;
    uint64_t result = shift_mod(1ull, mask, m);
    while (n >= 6)
    {
        n -= 6;
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        mask = (e >> n) & 0x3f;
        result = shift_mod(result, mask, m);
    }
    while (n > 0)
    {
        n -= 1;
        result = square_mod(result, m);
        if ((e >> n) & 1)
        {
            result <<= 1;
            result -= (result >= m) ? m : 0;
        }
    }
    return result;
}

// modular exponentiation a^e mod m
uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t n = e;
    uint64_t s = a;
    uint64_t result = 1;
    while (n)
    {
        if (n & 1)
            result = mul_mod(result, s, m);
        s = square_mod(s, m);
        n >>= 1;
    }
    return result;
}

struct barrett_t
{
    uint64_t m;   // modulus n bits
    uint64_t q;   // quotient 2^(3n/2) / m
    uint64_t r;   // remainder 2^(3n/2) % m
    uint64_t n;   // modulus size
    uint64_t n2;  // 1/2 modulus size
    uint64_t n32; // 3/2 modulus size
};

void barrett_precompute(struct barrett_t *p, uint64_t m)
{
    // precompute a variant of Barrett reduction
    p->m = m;
    p->n = 1 + uint64_log_2(m);
    if (p->n < 31)
    {
        // simple Barrett
        p->n2 = p->n << 1;
        p->n32 = 0;
        p->q = (1ull << p->n2) / m;
        return;
    }

    if (p->n < 42)
    {
        // modified Barrett
        p->n2 = p->n >> 1;
        p->n32 = p->n + p->n2;
        uint128_divrem(&p->q, &p->r, (uint128_t)1 << p->n32, m);
        return;
    }

    // no optimization
    p->n2 = 0;
    p->n32 = 0;
    p->q = 0;
    p->r = 0;
}

uint64_t barrett_mul_mod(uint64_t u, uint64_t v, const struct barrett_t &bn)
{
    if (bn.n < 31)
    {
        // for modulus up to 31 bits (3 multiplications)
	// assume u, v < n
	// makes r < n
        uint64_t r = u * v;
        uint64_t e = ((uint128_t)r * bn.q) >> bn.n2;
        r -= e * bn.m;
        return r >= bn.m ? bn.m : 0;
    }

    if (bn.n < 42)
    {
        // for modulus up to 42 bits  (4 multiplications)
	// assume u, v < 2^64
	// makes r < 2^64 
        uint128_t t = (uint128_t)u * v;
        uint64_t t_lo = t & ((1ull << bn.n32) - 1);    // up to 63 bits
        uint64_t t_hi = t >> bn.n32;                   // up to 21 bits
        uint64_t b = t_lo + t_hi * bn.r;               // up to 64 bits
        uint128_t e = ((uint128_t)b * bn.q) >> bn.n32; // up to 85 bits , down to 22 bits
        uint64_t s = (uint64_t)t - bn.m * (uint64_t)e; // barrett subtraction without underfllow
        return s;                                      // s < 3 * modulus
    }

    // no optimisation (2 instructions, including a long division)
    return mul_mod(u, v, bn.m);
}

// modular exponentiation a^e mod m with precomputations
// intermediate numbers are less than (3 * m)
uint64_t barrett_pow(uint64_t a, uint64_t e, const barrett_t &b)
{
    uint64_t bits = uint64_log_2(e) - 1;
    uint64_t result = a;
    while (bits--)
    {
        result = barrett_mul_mod(result, result, b);
        if ((e >> bits) & 1)
        {
            result = barrett_mul_mod(result, a, b);
        }
    }
    // final reduction, at most 3 subtractions.
    while (result >= b.m)
    {
        result -= b.m;
    }
    return result;
}

// MR strong test
static bool witness(uint64_t n, uint64_t s, uint64_t d, uint64_t a)
{
    uint64_t x, y;
    if (n == a)
        return true;

    if (a == 2)
    {
        x = pow2_mod(d, n);
    }
    else
    {
        x = pow_mod(a, d, n);
    }

    while (s)
    {
        y = square_mod(x, n);
        if (y == 1 && x != 1 && x != n - 1)
            return false;
        x = y;
        --s;
    }
    if (y != 1)
        return false;
    return true;
}

// deterministic primality test for n < 2^64.
// Assume that small factors are already processed, assume n > 2
bool uint64_is_prime(uint64_t n)
{
    uint64_t d = n / 2;
    uint64_t s = uint64_tzcnt(d);
    d >>= s++;

    if (n < 1373653)
        return witness(n, s, d, 2) && witness(n, s, d, 3);
    if (n < 9080191)
        return witness(n, s, d, 31) && witness(n, s, d, 73);
    if (n < 4759123141)
        return witness(n, s, d, 2) && witness(n, s, d, 7) && witness(n, s, d, 61);
    if (n < 1122004669633)
        return witness(n, s, d, 2) && witness(n, s, d, 13) && witness(n, s, d, 23) && witness(n, s, d, 1662803);
    if (n < 2152302898747)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11);
    if (n < 3474749660383)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13);
    if (n < 341550071728321)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17);
    if (n < 3825123056546413051)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19) &&
               witness(n, s, d, 23);
    // n < 318665857834031151167461
    return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
           witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19) &&
           witness(n, s, d, 23) && witness(n, s, d, 29) && witness(n, s, d, 31) && witness(n, s, d, 37);
}

// binary gcd
static uint64_t uint64_gcd(uint64_t u, uint64_t v)
{
    uint64_t t, k;

    if (u < v)
    {
        t = u;
        u = v;
        v = t;
    }
    if (v == 0)
        return u;

    // strip trailing zeroes
    k = uint64_tzcnt(u | v);
    v >>= k;

    // Stein algorithm, from odd number to odd number
    u >>= uint64_tzcnt(u);
    do
    {
        v >>= uint64_tzcnt(v);

        if (u > v)
        {
            t = u;
            u = v;
            v = t;
        }
        v -= u;
    } while (v);
    return u << k;
}

// binary gcd
static uint128_t uint128_gcd(uint128_t u, uint128_t v)
{
    uint128_t t, k;

    if (u < v)
    {
        t = u;
        u = v;
        v = t;
    }
    if (v == 0)
        return u;

    // strip trailing zeroes
    k = uint128_tzcnt(u | v);
    v >>= k;

    // Stein algorithm, from odd number to odd number
    u >>= uint128_tzcnt(u);
    do
    {
        v >>= uint128_tzcnt(v);

        if (u > v)
        {
            t = u;
            u = v;
            v = t;
        }
        v -= u;
    } while (v);
    return u << k;
}

// binary modular inverse 1/x mod m, with m odd and x < m, x and m coprime
static uint64_t uint64_mod_inv(uint64_t x, uint64_t m)
{
    if (x < 2)
        return x;
    if (m < 3)
        return 0;
    uint64_t a = x, b = m, u = 1, v = 0;
    while (a != 0)
    {
        unsigned za = uint64_tzcnt(a);
        a >>= za;
        while (za--)
        {
            u += (u & 1) ? m : 0;
            u >>= 1;
        }
        if (a < b)
        {
            uint64_t t = a;
            uint64_t s = u;
            a = b;
            u = v;
            b = t;
            v = s;
        }
        a -= b;
        u = (u >= v) ? u - v : u + m - v;
    }
    return b == 1 ? v : 0;
}

// binary modular inverse 1/x mod m, with m odd and x < m, x and m coprime
static uint128_t uint128_mod_inv(uint128_t x, uint128_t m)
{
    if (x < 2)
        return x;
    if (m < 3)
        return 0;
    uint64_t a = x, b = m, u = 1, v = 0;
    while (a != 0)
    {
        unsigned za = uint128_tzcnt(a);
        a >>= za;
        while (za--)
        {
            u += (u & 1) ? m : 0;
            u >>= 1;
        }
        if (a < b)
        {
            uint128_t t = a;
            uint128_t s = u;
            a = b;
            u = v;
            b = t;
            v = s;
        }
        a -= b;
        u = (u >= v) ? u - v : u + m - v;
    }
    return b == 1 ? v : 0;
}

static bool is_perfect_square(uint64_t a)
{
    if (0xffedfdfefdecull & (1ull << (a % 48)))
        return false;
    if (0xfdfdfdedfdfcfdecull & (1ull << (a % 64)))
        return false;
    if (0x7bfdb7cfedbafd6cull & (1ull << (a % 63)))
        return false;
    if (0x7dcfeb79ee35ccull & (1ull << (a % 55)))
        return false;
    if (0x8ec196bf5a60dc4ull & (1ull << (a % 61)))
        return false;
    if (0x5d49de7c1846d44ull & (1ull << (a % 59)))
        return false;
    if (0xd228fccfc512cull & (1ull << (a % 53)))
        return false;
    if (0x7bcae4d8ac20ull & (1ull << (a % 47)))
        return false;
    if (0x4a77c5c11acull & (1ull << (a % 43)))
        return false;
    if (0x4c7d4af8c8ull & (1ull << (a % 41)))
        return false;
    if (0x9a1dee164ull & (1ull << (a % 37)))
        return false;
    if (0x6de2b848ull & (1ull << (a % 31)))
        return false;
    if (0xc2edd0cull & (1ull << (a % 29)))
        return false;
    if (0x7acca0ull & (1ull << (a % 23)))
        return false;
    if (0x4f50cull & (1ull << (a % 19)))
        return false;
    if (0x5ce8ull & (1ull << (a % 17)))
        return false;
    if (0x9e4ull & (1ull << (a % 13)))
        return false;

    // approximation of square root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 2.0); // square root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    uint64_t c, m;
    // binary search (1 more bit of square root per iteration)
    uint64_t r = (uint64_t)d;
    uint64_t l = (uint64_t)dl;
    uint64_t h = (uint64_t)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m;
        if (c == a)
        {
            return true; // perfect square
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
    c = r * r; // check perfect square
    return (c == a);
}

static bool is_perfect_cube(uint64_t a)
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

    // approximation of cubic root with floating point accuracy
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

static bool is_perfect_sursolid(uint64_t a)
{
    if (0x1f7fef8fbff7cull & (1ull << (a % 50)))
        return false;
    if (0x7fcff5fe7fcull & (1ull << (a % 44)))
        return false;
    if (0xffa7efedfdf97fcull & (1ull << (a % 61)))
        return false;
    if (0xbef7ffbdf4ull & (1ull << (a % 41)))
        return false;
    if (0x39ffff9cull & (1ull << (a % 31)))
        return false;
    if (0x1249248ull & (1ull << (a % 27)))
        return false;
    if (0x40810204080ull & (1ull << (a % 49)))
        return false;

    // approximation of fifth root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 5.0); // fifth root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    uint64_t c, m;
    // binary search (1 more bit of fifth root per iteration)
    uint64_t r = (uint64_t)d;
    uint64_t l = (uint64_t)dl;
    uint64_t h = (uint64_t)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m * m * m * m;
        if (c == a)
        {
            return true; // perfect sursolid
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
    c = r * r * r * r * r; // check perfect sursolid
    return (c == a);
}

// pollard-rho factorization with brent variant
uint64_t brent_pollard_factor(uint64_t n)
{
    uint64_t i, x, ys, k;
    uint64_t m = 1000;
    uint64_t a = 2 + uint64_rnd() % (n - 4);
    uint64_t y = 1 + uint64_rnd() % (n - 2);
    uint64_t r = 1;
    uint64_t q = 1;
    uint64_t g = 1;

    do
    {
        x = y;
        for (i = 0; i < r; i++)
        {
            // y = y * y + a mod n
            y = square_add_mod(y, a, n);
        }

        k = 0;
        do
        {
            for (i = 0; i < m; i++)
            {
                ys = y;

                // y = y * y + a mod n
                y = square_add_mod(y, a, n);

                // q = q * |x-y| mod n
                q = mul_mod(q, (x > y) ? x - y : y - x, n);
            }
            g = uint64_gcd(q, n);
            k += m;
        } while (k < r && g == 1);

        r <<= 1;
    } while (g == 1);

    if (g == n)
    {
        // this can occur if one of gcd parameter is 0
        do
        {
            ys = square_add_mod(ys, a, n);
            g = uint64_gcd((x > ys) ? x - ys : ys - x, n);
        } while (g == 1);
    }

    return g;
}

struct factor_t
{
    uint64_t prime;
    uint64_t count;
};

bool factor_sort(const factor_t &f, const factor_t &g)
{
    return (f.prime < g.prime);
}

struct factor_find_t
{
    uint64_t prime;
    factor_find_t(uint64_t f) : prime(f)
    {
    }
    bool operator()(const factor_t &f) const
    {
        return f.prime == prime;
    }
};

typedef vector<factor_t> factor_v;
typedef vector<factor_t>::reverse_iterator reverse_iterator_v;

// search for large factors, assume sieving already done up to factor 151
// input n can be prime or composite
void uint64_large_factors(factor_v &primes, uint64_t n)
{
    unsigned i;
    uint64_t m;
    vector<uint64_t> factors;
    factors.push_back(n);

    do
    {
        m = factors.back();
        factors.pop_back();

        if (m == 1)
            continue;

        if (m < 157 * 157 || uint64_is_prime(m))
        {
            // avoid storing duplicate prime factors
            reverse_iterator_v v = find_if(primes.rbegin(), primes.rend(), factor_find_t(m));
            if (v == primes.rend())
            {
                // not found, add the factor
                factor_t f;
                f.prime = m;
                f.count = 1;
                primes.push_back(f);
            }
            else
            {
                v->count += 1;
            }
        }
        else
        {
            // m is not prime,
            // get more prime and composite factors from pollard-rho method
            uint64_t factor = brent_pollard_factor(m);
            factors.push_back(m / factor);
            factors.push_back(factor);
        }
    } while (factors.size());
}

uint64_t uint64_smallest_factor(uint64_t m)
{
    // first search a factor < 157
    uint64_t factor = uint64_small_factor(m);
    if (factor != 1)
    {
        // small prime factor < 157 found
        return factor;
    }
    if (m < 157 * 157)
    {
        // input number is prime
        return m;
    }
    // get factors in any order
    factor_v factors;
    uint64_large_factors(factors, m);
    // get the smallest factor from the vector
    sort(factors.begin(), factors.end(), factor_sort);
    return factors[0].prime;
}

bool uint64_all_factors(factor_v &factors, uint64_t m)
{
    // first search all prime factors < 157
    uint64_t factor = uint64_small_factor(m);
    while (factor != 1)
    {
        // avoid storing duplicate factors
        reverse_iterator_v v = find_if(factors.rbegin(), factors.rend(), factor_find_t(factor));
        if (v == factors.rend())
        {
            // not found, add the factor
            factor_t f;
            f.prime = factor;
            f.count = 1;
            factors.push_back(f);
        }
        else
        {
            v->count++;
        }
        m /= factor;
        factor = uint64_small_factor(m);
    }
    if (m < 157 * 157)
    {
	    // all factors < 157 have been removed
        // m is prime
        // avoid storing duplicate factors
        reverse_iterator_v v = find_if(factors.rbegin(), factors.rend(), factor_find_t(m));
        if (v == factors.rend())
        {
            // not found, add the factor
            factor_t f;
            f.prime = m;
            f.count = 1;
            factors.push_back(f);
        }
        else
        {
            v->count++;
        }
    }
    else
    {
        // get more factors in any order
        uint64_large_factors(factors, m);
    }

    // order the factors from the vector
    sort(factors.begin(), factors.end(), factor_sort);

    // search for perfect power
    for (uint64_t i = factors.size() - 1; i > 0; i -= 1)
    {
        if (factors[i].count != factors[0].count)
        {
            // not a perfect power for sure
            return false;
        }
    }
    // a perfect power
    return factors[0].count != 1;
}

// conversion of a large number into a basis-10 string.
// returns the output string length.
static unsigned uint128_sprint(char *ptr, uint128_t x)
{
    char linef[256];
    char *pb = linef;
    char *pt = ptr;

    if (x == 0)
    {
        *pb++ = '0';
    }
    while (x)
    {
        *pb++ = '0' + (char)(x % 10);
        x /= 10;
    }
    while (pb > linef)
    {
        *(pt++) = *(--pb);
    }
    *pt = 0;
    return pt - ptr;
}

unsigned int128_sprint(char *ptr, int128_t x)
{
    char *pt = ptr;

    if (x < 0)
    {
        x = -x;
        *pt++ = '-';
    }
    pt += uint128_sprint(pt, x);
    return pt - ptr;
}

