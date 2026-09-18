

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

typedef unsigned _BitInt(256) uint256_t;
typedef signed _BitInt(256) int256_t;

static bool uint64_is_perfect_square(uint64_t a);
static bool uint64_is_perfect_cube(uint64_t a);
static bool uint64_is_perfect_sursolid(uint64_t a);
static bool uint64_is_perfect_power(uint64_t a);

static bool uint128_is_perfect_square(uint128_t a);
static bool uint128_is_perfect_cube(uint128_t a);

// --------------------------------------------------------------------------------------
//
// C optimized implementations of simple utilities
//
// --------------------------------------------------------------------------------------

// randomness generation with period 2^64
static uint64_t uint64_rnd(void)
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

static inline uint128_t uint256_long_mod(uint256_t u, uint128_t n)
{
    return (uint128_t)(u % n);
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
    asm("mulq %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    uint128_t tmp = (uint128_t)a * b;
    tmp %= n;
    return (uint64_t)tmp;
#endif
}

// (a * b + c) % n
static inline uint64_t mul_add_mod(uint64_t a, uint64_t b, uint64_t c, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mulq %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("addq %2, %1\n\tadcq $0, %0" : "+d"(r), "+a"(a) : "r"(c) : "flags");
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    uint128_t tmp = (uint128_t)a * b;
    tmp += c;
    tmp %= n;
    return (uint64_t)tmp;
#endif
}

static inline uint128_t uint128_mul_mod(uint128_t a, uint128_t b, uint128_t n)
{
    uint256_t tmp = a;
    tmp *= b;
    tmp %= n;
    return (uint128_t)tmp;
}

static inline uint128_t uint128_mul_add_mod(uint128_t a, uint128_t b, uint128_t c, uint128_t n)
{
    uint256_t tmp = a;
    tmp *= b;
    tmp += c;
    tmp %= n;
    return (uint128_t)tmp;
}

static inline uint64_t square_mod(uint64_t a, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mulq %2" : "=d"(r), "=a"(a) : "1"(a));
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n) : "flags");
    return r;
#else
    uint128_t tmp = (uint128_t)a * a;
    tmp %= n;
    return tmp;
#endif
}

static inline uint128_t uint128_square_mod(uint128_t a, uint128_t n)
{
    uint256_t tmp = a;
    tmp *= a;
    tmp %= n;
    return (uint128_t)tmp;
}

static inline uint128_t uint128_square_add_mod(uint128_t a, uint128_t c, uint128_t n)
{
    uint256_t tmp = a;
    tmp *= a;
    tmp += c;
    tmp %= n;
    return (uint128_t)tmp;
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
    asm("xorq %0, %0\n shldq %b3, %1, %0\n shlxq %3, %1, %%rax\n divq %2"
        : "=&d"(r)
        : "r"(u), "r"(n), "c"(s)
        : "flags", "%rax");
    return r;
#else
    uint128_t t = (uint128_t)u;
    t <<= s;
    return (uint64_t)(t % n);
#endif
}

static inline uint128_t uint128_shift_mod(uint128_t u, uint64_t s, uint128_t n)
{
    if (s >= 128 || u >> (128 - s))
    {
    uint256_t t = u;
    t <<= s;
    return (uint128_t)(t % n);
    }
		    else
		    {
    uint128_t t = u;
    t <<= s;
    return (uint128_t)(t % n);
		    }
}

// count leading zeroed bits
static inline uint64_t uint64_lzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("lzcntq %1,%0" : "=r"(r) : "r"(a));
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
    asm("tzcntq %1,%0" : "=r"(r) : "r"(a));
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

// Euler's totient function
static uint64_t uint64_phi(uint64_t m)
{
    uint64_t result = m;
    uint64_t x = m;
    uint64_t p = 2;

    while (p * p <= x)
    {
        if (x % p == 0)
        {
            do
            {
                x /= p;
            } while (x % p == 0);
            result -= result / p;
        }
        p += 1;
    }

    if (x > 1)
    {
        result -= result / x;
    }
    return result;
}

// assume a >= 0
// assume odd b > 0
// assume a < 2^62 or b < 2^62
static int uint64_jacobi(uint64_t a, uint64_t b)
{
    static char cols[]
                    [64] = {
                        {1},                                                                               // a=0
                        {1, 1},                                                                            // a=1
                        {(char)-1, (char)-1},                                                              // a=2
                        {0, (char)-1, (char)-1},                                                           // a=3
                        {1, 1, 1, 1},                                                                      // a=4
                        {(char)-1, 0, (char)-1, 1, 1},                                                     // a=5
                        {0, 1, (char)-1, 0, (char)-1, (char)-1},                                           // a=6
                        {1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1},                                 // a=7
                        {(char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1},                              // a=8
                        {0, 1, 1, 0, 1, 1, 0, 1, 1},                                                       // a=9
                        {1, 0, (char)-1, 1, (char)-1, 1, 0, (char)-1, (char)-1, (char)-1},                 // a=10
                        {(char)-1, 1, 1, 1, 0, (char)-1, (char)-1, (char)-1, 1, (char)-1, (char)-1},       // a=11
                        {0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1},                  // a=12
                        {1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1}, // a=13
                        {(char)-1, 1, 0, 1, 1, 1, (char)-1, (char)-1, (char)-1,
                         0, (char)-1, 1, (char)-1, (char)-1},                                             // a=14
                        {0, 0, 1, 0, 1, (char)-1, 0, 1, (char)-1, 0, (char)-1, 0, 0, (char)-1, (char)-1}, // a=15
                        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},                                 // a=16
                        {(char)-1, (char)-1, (char)-1, 1, (char)-1, 1, 1, 0, 1, 1, (char)-1, 1, (char)-1, (char)-1,
                         (char)-1, 1, 1}, // a=17
                        {0, (char)-1, 1, 0, (char)-1, (char)-1,
                         0, 1, (char)-1, 0, 1, 1, 0, (char)-1, 1, 0, (char)-1, (char)-1}, // a=18
                        {1, 1, (char)-1, 1, (char)-1, (char)-1, 1, 1, 0, (char)-1, (char)-1,
                         1, 1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1}, // a=19
                        {(char)-1, 0, (char)-1, 1, 1, (char)-1, 0, (char)-1, 1, 1,
                         (char)-1, 0, (char)-1, 1, 1, (char)-1, 0, (char)-1, 1, 1}, // a=20
                        {0, 1, 0,        0,        (char)-1, (char)-1, 0, 1, (char)-1, 0, (char)-1,
                         1, 0, (char)-1, (char)-1, 0,        0,        1, 0, 1,        1}, // a=21
                        {1, (char)-1, 1, 1,        0,        1, (char)-1, (char)-1, (char)-1, 1, (char)-1, 1, 1,
                         1, (char)-1, 0, (char)-1, (char)-1, 1, (char)-1, (char)-1, (char)-1}, // a=22
                        {(char)-1, (char)-1, 1,        1, 1,        1,        1,        (char)-1,
                         1,        (char)-1, 0,        1, (char)-1, 1,        (char)-1, (char)-1,
                         (char)-1, (char)-1, (char)-1, 1, 1,        (char)-1, (char)-1}, // a=23
                        {0, 1, (char)-1, 0, (char)-1, (char)-1, 0, (char)-1, 1, 0, 1, 1,
                         0, 1, (char)-1, 0, (char)-1, (char)-1, 0, (char)-1, 1, 0, 1, 1},            // a=24
                        {1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1}, // a=25
                        {(char)-1, 1,        (char)-1, 1,        1,        0,        (char)-1, 1,        1,
                         1,        1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1,
                         0,        (char)-1, (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1}, // a=26
                        {0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1,
                         0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1}, // a=27
                        {1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 0, (char)-1, 1, 1, 1,
                         1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 0, (char)-1, 1, 1, 1}, // a=28
                        {(char)-1, 1, 1,        1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1,
                         1,        1, (char)-1, 0, (char)-1, 1, 1,        (char)-1, (char)-1, (char)-1,
                         (char)-1, 1, (char)-1, 1, 1,        1, (char)-1, 1,        1}, // a=29
                        {0, 0, 1, 0, (char)-1, 1,        0, 1,        1, 0, (char)-1, 0, 0, 1,        (char)-1,
                         0, 0, 1, 0, (char)-1, (char)-1, 0, (char)-1, 1, 0, (char)-1, 0, 0, (char)-1, (char)-1}, // a=30
                        {1,        1, (char)-1, 1,        1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                         1,        1, (char)-1, 0,        1, (char)-1, (char)-1, (char)-1, 1,        1,        1,
                         (char)-1, 1, (char)-1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1}, // a=31
                        {(char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1,
                         (char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1,
                         (char)-1, (char)-1, 1, 1, (char)-1, (char)-1, 1, 1}, // a=32
                        {0, (char)-1, (char)-1, 0, 0,        (char)-1, 0, 1, (char)-1, 0, (char)-1, 1,
                         0, 1,        1,        0, 1,        1,        0, 1, (char)-1, 0, (char)-1, 1,
                         0, (char)-1, 0,        0, (char)-1, (char)-1, 0, 1, 1}, // a=33
                        {1,        1,        (char)-1, 1,        1,        (char)-1, 1,       0,        (char)-1,
                         (char)-1, (char)-1, 1,        1,        1,        (char)-1, 1,       (char)-1, 1,
                         (char)-1, (char)-1, (char)-1, 1,        1,        1,        0,       (char)-1, 1,
                         (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1}, // a=34
                        {(char)-1, 0,        0,        1,        (char)-1, 1, 0,        1,        1,
                         0,        1,        0,        (char)-1, 1,        1, 1,        0,        (char)-1,
                         (char)-1, (char)-1, 1,        0,        (char)-1, 0, (char)-1, (char)-1, 0,
                         (char)-1, 1,        (char)-1, 0,        0,        1, (char)-1, (char)-1}, // a=35
                        {0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1,
                         0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1}, // a=36
                        {1,        (char)-1, 1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, 1,
                         (char)-1, 1,        1, (char)-1, (char)-1, 1,        (char)-1, 0,        (char)-1, 1,
                         (char)-1, (char)-1, 1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1,
                         1,        1,        1, (char)-1, 1,        1,        1}, // a=37
                        {(char)-1, (char)-1, (char)-1, 1, 1,        1,        1,        1,       0,        1,
                         (char)-1, 1,        (char)-1, 1, 1,        (char)-1, 1,        1,       (char)-1, (char)-1,
                         1,        (char)-1, (char)-1, 1, (char)-1, 1,        (char)-1, 0,       (char)-1, (char)-1,
                         (char)-1, (char)-1, (char)-1, 1, 1,        1,        (char)-1, (char)-1}, // a=38
                        {0,        1, 1,        0,        (char)-1, 0,        0,        (char)-1, 1,        0,
                         1,        1, 0,        (char)-1, 1,        0,        1,        (char)-1, 0,        1,
                         (char)-1, 0, (char)-1, 1,        0,        (char)-1, (char)-1, 0,        (char)-1, 1,
                         0,        0, 1,        0,        (char)-1, (char)-1, 0,        (char)-1, (char)-1}, // a=39
                        {1,        0, (char)-1, 1,        (char)-1, 1,        0, (char)-1, (char)-1, (char)-1,
                         (char)-1, 0, 1,        (char)-1, 1,        (char)-1, 0, 1,        1,        1,
                         1,        0, (char)-1, 1,        (char)-1, 1,        0, (char)-1, (char)-1, (char)-1,
                         (char)-1, 0, 1,        (char)-1, 1,        (char)-1, 0, 1,        1,        1}, // a=40
                        {(char)-1, 1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1,
                         1,        1, 1,        (char)-1, (char)-1, 1,        1,        (char)-1, 1,
                         1,        0, 1,        1,        (char)-1, 1,        1,        (char)-1, (char)-1,
                         1,        1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1,
                         (char)-1, 1, (char)-1, 1,        1}, // a=41
                        {0, (char)-1, 0,        0,        1,        1, 0,        1,        1,       0,        (char)-1,
                         1, 0,        1,        (char)-1, 0,        0, (char)-1, 0,        1,       (char)-1, 0,
                         1, 0,        0,        1,        (char)-1, 0, (char)-1, 1,        0,       (char)-1, (char)-1,
                         0, (char)-1, (char)-1, 0,        0,        1, 0,        (char)-1, (char)-1}, // a=42
                        {1, (char)-1, 1,        1,        (char)-1, 1,        (char)-1, 1,        1,
                         1, (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1,
                         1, 1,        0,        (char)-1, (char)-1, 1,        1,        1,        1,
                         1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1,
                         1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1}, // a=43
                        {(char)-1, 1,        1,        1,        0,        (char)-1, (char)-1, (char)-1, 1,
                         (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 0,        1,        1,
                         1,        (char)-1, 1,        1,        (char)-1, 1,        1,        1,        0,
                         (char)-1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1,
                         (char)-1, 0,        1,        1,        1,        (char)-1, 1,        1}, // a=44
                        {0, 0, (char)-1, 0, 1, (char)-1, 0, (char)-1, 1, 0, (char)-1, 0, 0, 1, 1,
                         0, 0, (char)-1, 0, 1, (char)-1, 0, (char)-1, 1, 0, (char)-1, 0, 0, 1, 1,
                         0, 0, (char)-1, 0, 1, (char)-1, 0, (char)-1, 1, 0, (char)-1, 0, 0, 1, 1}, // a=45
                        {1,        1,        1,        1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, 1,
                         0,        1,        1,        (char)-1, (char)-1, (char)-1, 1,        1,        (char)-1, 1,
                         (char)-1, 1,        (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1, 1,        1,
                         1,        (char)-1, (char)-1, 0,        (char)-1, 1,        1,        (char)-1, 1,        1,
                         (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, (char)-1}, // a=46
                        {(char)-1, (char)-1, (char)-1, 1,        1,        (char)-1, 1,        1,
                         1,        1,        1,        1,        (char)-1, (char)-1, 1,        (char)-1,
                         1,        1,        1,        (char)-1, 1,        (char)-1, 0,        1,
                         (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,
                         1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1,
                         (char)-1, (char)-1, 1,        1,        1,        (char)-1, (char)-1}, // a=47
                        {0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1,
                         0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1,
                         0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1, (char)-1, 0, 1, 1}, // a=48
                        {1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
                         1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}, // a=49
                        {(char)-1, 0, 1,        1,        (char)-1, (char)-1, 0, 1,        (char)-1, (char)-1,
                         1,        0, (char)-1, (char)-1, 1,        1,        0, (char)-1, 1,        1,
                         (char)-1, 0, 1,        1,        (char)-1, (char)-1, 0, 1,        (char)-1, (char)-1,
                         1,        0, (char)-1, (char)-1, 1,        1,        0, (char)-1, 1,        1,
                         (char)-1, 0, 1,        1,        (char)-1, (char)-1, 0, 1,        (char)-1, (char)-1}, // a=50
                        {0, 1,        1,        0, (char)-1, 1,        0, 0,        (char)-1, 0, (char)-1, 1,
                         0, 1,        1,        0, 1,        (char)-1, 0, 1,        (char)-1, 0, 1,        1,
                         0, (char)-1, (char)-1, 0, 1,        (char)-1, 0, 1,        (char)-1, 0, (char)-1, (char)-1,
                         0, (char)-1, 1,        0, 1,        0,        0, (char)-1, 1,        0, (char)-1, (char)-1,
                         0, (char)-1, (char)-1}, // a=51
                        {1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1,
                         1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1,
                         1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1,
                         1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1}, // a=52
                        {(char)-1, (char)-1, 1,        1,        1,        1,        1,        1,        (char)-1,
                         (char)-1, (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                         (char)-1, (char)-1, 1,        (char)-1, 1,        1,        (char)-1, 0,        (char)-1,
                         1,        1,        (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1,
                         (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        1,
                         1,        1,        1,        1,        (char)-1, (char)-1, 1,        1}, // a=53
                        {0, 1, (char)-1, 0,        (char)-1, (char)-1, 0,        (char)-1, 1,        0,        1,
                         1, 0, 1,        (char)-1, 0,        (char)-1, (char)-1, 0,        (char)-1, 1,        0,
                         1, 1, 0,        1,        (char)-1, 0,        (char)-1, (char)-1, 0,        (char)-1, 1,
                         0, 1, 1,        0,        1,        (char)-1, 0,        (char)-1, (char)-1, 0,        (char)-1,
                         1, 0, 1,        1,        0,        1,        (char)-1, 0,        (char)-1, (char)-1}, // a=54
                        {1,        0, (char)-1, 1,        0,        1,        0, 1,        1,        (char)-1,
                         1,        0, 1,        (char)-1, (char)-1, 0,        0, (char)-1, 1,        (char)-1,
                         (char)-1, 0, 1,        1,        1,        (char)-1, 0, 1,        (char)-1, (char)-1,
                         (char)-1, 0, 1,        1,        (char)-1, 1,        0, 0,        1,        1,
                         (char)-1, 0, (char)-1, 1,        (char)-1, (char)-1, 0, (char)-1, 0,        (char)-1,
                         1,        0, (char)-1, (char)-1, (char)-1}, // a=55
                        {(char)-1, 1,        0,        1,        1,        1,        (char)-1, (char)-1,
                         (char)-1, 0,        (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1,
                         0,        (char)-1, (char)-1, (char)-1, 1,        1,        1,        0,
                         1,        (char)-1, 1,        1,        (char)-1, 1,        0,        1,
                         1,        1,        (char)-1, (char)-1, (char)-1, 0,        (char)-1, 1,
                         (char)-1, (char)-1, 1,        (char)-1, 0,        (char)-1, (char)-1, (char)-1,
                         1,        1,        1,        0,        1,        (char)-1, 1,        1}, // a=56
                        {0, (char)-1, 1,        0, (char)-1, (char)-1, 0, (char)-1, 0,        0, (char)-1, 1,
                         0, 1,        (char)-1, 0, (char)-1, (char)-1, 0, 1,        1,        0, (char)-1, 1,
                         0, 1,        1,        0, 1,        1,        0, 1,        (char)-1, 0, 1,        1,
                         0, (char)-1, (char)-1, 0, (char)-1, 1,        0, 1,        (char)-1, 0, 0,        (char)-1,
                         0, (char)-1, (char)-1, 0, 1,        (char)-1, 0, 1,        1}, // a=57
                        {1,        (char)-1, 1,        1,        1,        (char)-1, (char)-1, (char)-1, 1,
                         1,        1,        1,        1,        0,        (char)-1, 1,        (char)-1, 1,
                         (char)-1, (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1,
                         1,        (char)-1, 1,        1,        1,        (char)-1, 1,        1,        (char)-1,
                         1,        1,        (char)-1, 1,        (char)-1, 1,        0,        (char)-1, (char)-1,
                         (char)-1, (char)-1, (char)-1, 1,        1,        1,        (char)-1, (char)-1, (char)-1,
                         1,        (char)-1, (char)-1, (char)-1}, // a=58
                        {(char)-1, 1,        (char)-1, 1,        1,        (char)-1, (char)-1, 1,        (char)-1,
                         1,        1,        1,        (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1,
                         1,        1,        1,        1,        1,        1,        (char)-1, 1,        1,
                         1,        0,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1,
                         (char)-1, (char)-1, (char)-1, 1,        1,        1,        (char)-1, (char)-1, 1,
                         (char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        1,        (char)-1, (char)-1,
                         1,        (char)-1, 1,        (char)-1, (char)-1}, // a=59
                        {0, 0,        1,        0, 1,        (char)-1, 0, 1,        (char)-1, 0, (char)-1, 0,
                         0, (char)-1, (char)-1, 0, 0,        (char)-1, 0, (char)-1, 1,        0, (char)-1, 1,
                         0, 1,        0,        0, 1,        1,        0, 0,        1,        0, 1,        (char)-1,
                         0, 1,        (char)-1, 0, (char)-1, 0,        0, (char)-1, (char)-1, 0, 0,        (char)-1,
                         0, (char)-1, 1,        0, (char)-1, 1,        0, 1,        0,        0, 1,        1}, // a=60
                        {1,        1,        (char)-1, 1,        (char)-1, 1,        1,        (char)-1, 1,
                         (char)-1, (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1,
                         1,        1,        (char)-1, 1,        1,        1,        (char)-1, (char)-1, (char)-1,
                         1,        (char)-1, 0,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                         1,        1,        (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1,
                         (char)-1, 1,        1,        (char)-1, (char)-1, 1,        (char)-1, 1,        1,
                         (char)-1, 1,        (char)-1, 1,        1,        1,        1}, // a=61
                        {(char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        1,        (char)-1, 1,
                         1,        1,        1,        (char)-1, 1,        0,        1,        1,        1,
                         (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        1,        1,        1,
                         (char)-1, 1,        1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1,
                         (char)-1, 1,        1,        1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1,
                         0,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, 1,        (char)-1,
                         (char)-1, 1,        (char)-1, 1,        1,        1,        (char)-1, (char)-1}, // a=62
                        {0, (char)-1, 0, 0, (char)-1, (char)-1, 0,        (char)-1, 1, 0, (char)-1, 1, 0, 1,
                         1, 0,        0, 1, 0,        (char)-1, (char)-1, 0,        1, 0, 0,        1, 1, 0,
                         1, (char)-1, 0, 1, (char)-1, 0,        (char)-1, (char)-1, 0, 0, (char)-1, 0, 1, 1,
                         0, (char)-1, 0, 0, (char)-1, (char)-1, 0,        (char)-1, 1, 0, (char)-1, 1, 0, 1,
                         1, 0,        0, 1, 0,        (char)-1, (char)-1}, // a=63
                    };
    static char rows[][64] = {
        {1, 1},                                                                            // b=1
        {(char)-1, 0, 1},                                                                  // b=3
        {(char)-1, (char)-1, 1, 0, 1},                                                     // b=5
        {1, (char)-1, 1, (char)-1, (char)-1, 0, 1},                                        // b=7
        {1, 0, 1, 1, 0, 1, 1, 0, 1},                                                       // b=9
        {(char)-1, 1, 1, 1, (char)-1, (char)-1, (char)-1, 1, (char)-1, 0, 1},              // b=11
        {(char)-1, 1, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 1, (char)-1, 1, 0, 1}, // b=13
        {1, 0, 1, 0, 0, (char)-1, 1, 0, 0, (char)-1, 0, (char)-1, (char)-1, 0, 1},         // b=15
        {1, (char)-1, 1, (char)-1, (char)-1, (char)-1, 1, 1, (char)-1, (char)-1, (char)-1, 1, (char)-1, 1, 1, 0,
         1}, // b=17
        {(char)-1, (char)-1, 1, 1, 1, 1, (char)-1, 1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 1,
         (char)-1, 0, 1}, // b=19
        {(char)-1, 0, 1, 1, 0, 0, (char)-1, 0, (char)-1, (char)-1, 0,
         (char)-1, 0, 0, 1, 1, 0, (char)-1, 1, 0,        1}, // b=21
        {1,        1,        1, (char)-1, 1, (char)-1, 1,        1,        (char)-1, (char)-1, 1, 1,
         (char)-1, (char)-1, 1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1, 0,        1}, // b=23
        {1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1},               // b=25
        {(char)-1, 0,        1, (char)-1, 0,        1, (char)-1, 0,        1, (char)-1, 0,        1, (char)-1, 0,
         1,        (char)-1, 0, 1,        (char)-1, 0, 1,        (char)-1, 0, 1,        (char)-1, 0, 1}, // b=27
        {(char)-1, (char)-1, 1,        1,        1,        1,        (char)-1, 1,        (char)-1, (char)-1,
         (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1,
         1,        1,        1,        1,        (char)-1, (char)-1, 1,        0,        1}, // b=29
        {1,        (char)-1, 1,        1,        (char)-1, 1,        1,        1, 1,        (char)-1, (char)-1,
         (char)-1, 1,        (char)-1, 1,        (char)-1, 1,        1,        1, (char)-1, (char)-1, (char)-1,
         (char)-1, 1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, 0, 1}, // b=31
        {1,        0, 1, (char)-1, 0, (char)-1, 1,        0, (char)-1, 0, 0, (char)-1, (char)-1, 0, 1, 1, 0, (char)-1,
         (char)-1, 0, 0, (char)-1, 0, 1,        (char)-1, 0, (char)-1, 1, 0, 1,        1,        0, 1}, // b=33
        {(char)-1, 1, 1, 0, (char)-1, 0,        (char)-1, 1, 0,        1,        1,        1,
         0,        0, 1, 1, (char)-1, (char)-1, 0,        0, (char)-1, (char)-1, (char)-1, 0,
         (char)-1, 1, 0, 1, 0,        (char)-1, (char)-1, 1, (char)-1, 0,        1}, // b=35
        {(char)-1, 1,        1,        (char)-1, (char)-1, 1,        (char)-1, 1,        1,        1,
         1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, 1,
         (char)-1, (char)-1, (char)-1, 1,        1,        1,        1,        (char)-1, 1,        (char)-1,
         (char)-1, 1,        1,        (char)-1, 1,        0,        1}, // b=37
        {1, 0,        1, 1,        0,        (char)-1, 1, 0,        1,        1, 0, 0, (char)-1, 0,
         1, (char)-1, 0, (char)-1, 1,        0,        1, (char)-1, 0,        1, 0, 0, (char)-1, (char)-1,
         0, (char)-1, 1, 0,        (char)-1, (char)-1, 0, (char)-1, (char)-1, 0, 1}, // b=39
        {1,        (char)-1, 1,        1,        (char)-1, (char)-1, 1,        1, 1, (char)-1, (char)-1,
         (char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        (char)-1, 1, 1, (char)-1, 1,
         (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 1, 1,        (char)-1,
         (char)-1, 1,        1,        (char)-1, 1,        1,        0,        1}, // b=41
        {(char)-1, (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1, 1,        1,        1,        (char)-1,
         1,        1,        1,        1,        1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,
         1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1,
         1,        1,        (char)-1, 1,        (char)-1, 1,        1,        (char)-1, 0,        1}, // b=43
        {(char)-1, 0, 1, 0, 0, (char)-1, (char)-1, 0, 0, 1, 0, (char)-1, 1, 0, 1,
         (char)-1, 0, 1, 0, 0, (char)-1, (char)-1, 0, 0, 1, 0, (char)-1, 1, 0, 1,
         (char)-1, 0, 1, 0, 0, (char)-1, (char)-1, 0, 0, 1, 0, (char)-1, 1, 0, 1}, // b=45
        {1,        1,        1,        (char)-1, 1,        1, 1,        1,        (char)-1, (char)-1,
         1,        (char)-1, 1,        (char)-1, 1,        1, 1,        (char)-1, (char)-1, 1,
         (char)-1, (char)-1, 1,        1,        (char)-1, 1, 1,        (char)-1, (char)-1, (char)-1,
         1,        (char)-1, 1,        (char)-1, 1,        1, (char)-1, (char)-1, (char)-1, (char)-1,
         1,        (char)-1, (char)-1, (char)-1, (char)-1, 0, 1}, // b=47
        {1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
         1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1}, // b=49
        {(char)-1, 0,        1,        1,        0, (char)-1, (char)-1, 0,        (char)-1, 1, 0,        1, 1,        0,
         1,        0,        0,        1,        1, 0,        (char)-1, 1,        0,        1, (char)-1, 0, (char)-1, 1,
         0,        (char)-1, (char)-1, 0,        0, (char)-1, 0,        (char)-1, (char)-1, 0, (char)-1, 1, 0,        1,
         1,        0,        (char)-1, (char)-1, 0, 1,        (char)-1, 0,        1}, // b=51
        {(char)-1, (char)-1, 1,        (char)-1, 1,        1,        (char)-1, 1,        1,        1,        (char)-1,
         1,        (char)-1, 1,        1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, (char)-1,
         1,        1,        (char)-1, (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1,
         (char)-1, 1,        1,        1,        (char)-1, 1,        (char)-1, 1,        1,        1,        (char)-1,
         1,        1,        (char)-1, 1,        (char)-1, (char)-1, 1,        0,        1}, // b=53
        {1,        (char)-1, 1,        0,        (char)-1, 1,        1,        1,        0,        0, (char)-1,
         1,        1,        0,        1,        1,        1,        (char)-1, 0,        (char)-1, 0, (char)-1,
         (char)-1, 0,        1,        (char)-1, 1,        (char)-1, 0,        1,        1,        0, 1,
         0,        1,        (char)-1, (char)-1, (char)-1, 0,        (char)-1, (char)-1, 1,        0, 0,
         (char)-1, (char)-1, (char)-1, 1,        0,        (char)-1, 1,        (char)-1, (char)-1, 0, 1}, // b=55
        {1,        0, 1,        (char)-1, 0, 1,        1,        0, (char)-1, (char)-1, 0, (char)-1,
         1,        0, 1,        (char)-1, 0, 0,        (char)-1, 0, (char)-1, (char)-1, 0, 1,
         (char)-1, 0, 1,        1,        0, (char)-1, 1,        0, (char)-1, (char)-1, 0, (char)-1,
         0,        0, (char)-1, 1,        0, 1,        (char)-1, 0, (char)-1, (char)-1, 0, 1,
         1,        0, (char)-1, 1,        0, 1,        1,        0, 1}, // b=57
        {(char)-1, 1,        1,        1,        (char)-1, 1,        (char)-1, 1,        (char)-1, (char)-1,
         1,        (char)-1, (char)-1, 1,        1,        1,        (char)-1, 1,        1,        1,
         1,        (char)-1, (char)-1, 1,        1,        1,        1,        1,        (char)-1, (char)-1,
         (char)-1, (char)-1, (char)-1, 1,        1,        (char)-1, (char)-1, (char)-1, (char)-1, 1,
         (char)-1, (char)-1, (char)-1, 1,        1,        (char)-1, 1,        1,        (char)-1, 1,
         (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, 0,        1}, // b=59
        {(char)-1, 1, 1,        1,        (char)-1, (char)-1, (char)-1, 1,        (char)-1, (char)-1, 1,
         1,        1, 1,        1,        (char)-1, (char)-1, 1,        1,        (char)-1, 1,        (char)-1,
         (char)-1, 1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, (char)-1, 1,
         (char)-1, 1, (char)-1, (char)-1, 1,        (char)-1, 1,        1,        (char)-1, (char)-1, 1,
         1,        1, 1,        1,        (char)-1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
         1,        1, (char)-1, 1,        0,        1}, // b=61
        {1, 0,        1, (char)-1, 0,        0, 1, 0, (char)-1, 1, 0,        (char)-1, 0, 0,
         1, (char)-1, 0, (char)-1, (char)-1, 0, 1, 1, 0,        1, (char)-1, 0,        0, 1,
         0, (char)-1, 1, 0,        (char)-1, 0, 0, 1, (char)-1, 0, (char)-1, (char)-1, 0, 1,
         1, 0,        1, (char)-1, 0,        0, 1, 0, (char)-1, 1, 0,        (char)-1, 0, 0,
         1, (char)-1, 0, (char)-1, (char)-1, 0, 1}, // b=63
    };

    if (a >= b)
    {
        a %= b;
    }

    // 0 <= a < b with b odd
    if (a < 3)
    {
        // (0/b) = (b == 1)
        if (a == 0)
        {
            return (b == 1) ? 1 : 0;
        }
        // (1/b) = 1
        if (a == 1)
        {
            return 1;
        }
        // (2/b) = (b % 8 == 3 || b % 8 == 5) ? -1 : 1;
        return (((b >> 2) ^ (b >> 1)) & 1) ? -1 : 1;
    }

    // 3 <= a < b with b odd
    if (a < 64)
    {
        if (a & 2)
        {
            uint64_t k = ((b - 3) >> 1) % (2 * a);
            return k < a ? cols[a][k] : -cols[a][k - a];
        }
        else
        {
            uint64_t k = ((b - 3) >> 1) % a;
            return cols[a][k];
        }
    }

    // 64 <= a < b with b odd and a < 2^62
    int t = 1;
    if ((a & 3) == 0)
    {
        b = b % a;
    }
    else if ((a & 3) == 2 && b < 3 * a)
    {
        if (b < 2 * a)
        {
            t = (b & 2) ? -t : t;
            b = b - a;
        }
        else
        {
            t = (b & 1) ? -t : t;
            b = b - 2 * a;
        }
    }
    else
    {
        b = b % (4 * a);
        if (b >= 2 * a)
        {
            b -= 2 * a;
            t = (a & 2) ? -t : t;
        }
    }

    // 64 <= a, b odd
    if (b < 64)
    {
        uint64_t k = (a - 2) % b;
        return (t == -1) ? -rows[b >> 1][k] : rows[b >> 1][k];
    }

    unsigned c = (b >> 2) ^ (b >> 1);
    while (a)
    {
        unsigned v = uint64_tzcnt(a);
        a >>= v;
        t = (c & v & 1) ? -t : t;

        if (a < b)
        {
            uint64_t k = a;
            a = b;
            b = k;
            t = ((a & b & 3) == 3) ? -t : t;
            c = (b >> 2) ^ (b >> 1);
        }

        a -= b;
    }
    return (b == 1) ? t : 0;
}

static int int64_kronecker(int64_t a, int64_t b)
{
    int v;
    int t = 1;
    if (b == 0)
    {
        return (a == -1 || a == 1) ? 1 : 0;
    }

    // make b positive
    // K(a,-b) = K(a,b)*sgn(a)
    if (b < 0)
    {
        b = -b;
        t = (a < 0) ? -t : t;
    }

    // make a positive
    // K(-a,b) = K(a,b)*(-1)^(b'>>1)
    if (a < 0)
    {
        a = -a;
        v = uint64_tzcnt(b);
        t = ((b >> v) & 2) ? -t : t;
    }

    // make b odd
    if ((b & 1) == 0)
    {
        if ((a & 1) == 0)
        {
            return 0;
        }
        v = uint64_tzcnt(b);
        if (v & 1)
        {
            t = (((a >> 1) ^ (a >> 2)) & 1) ? -t : t;
        }
        b >>= v;
    }

    // make a odd
    if ((a & 1) == 0)
    {
        if (a == 0)
        {
            return (b == 1) ? 1 : 0;
        }
        v = uint64_tzcnt(a);
        if (v & 1)
        {
            t = (((b >> 1) ^ (b >> 2)) & 1) ? -t : t;
        }
        a >>= v;
    }

    v = uint64_jacobi(a, b);
    return (t < 0) ? -v : v;
}

static int int128_kronecker(int128_t a, int128_t b)
{
    unsigned v;
    int t = 1;
    if (a == 0)
    {
        return (b == -1 || b == 1) ? 1 : 0;
    }

    if (b == 0)
    {
        return (a == -1 || a == 1) ? 1 : 0;
    }

    // K(a,-b) = K(a,b)*sgn(a)
    if (b < 0)
    {
        b = -b;
        t = (a < 0) ? -t : t;
    }
    // K(-a,b) = K(a,b)*(-1)^(b'>>1)
    if (a < 0)
    {
        a = -a;
        v = uint128_tzcnt(b);
        t = ((b >> v) & 2) ? -t : t;
    }

    // make b odd
    if ((b & 1) == 0)
    {
        if ((a & 1) == 0)
        {
            return 0;
        }
        v = uint128_tzcnt(b);
        if (v & 1)
        {
            t = (((a >> 1) ^ (a >> 2)) & 1) ? -t : t;
        }
        b >>= v;
    }

    unsigned c = (b >> 2) ^ (b >> 1);
    while ((a | b) >> 62 != 0)
    {
        v = uint128_tzcnt(a);
        a >>= v;
        t = (c & v & 1) ? -t : t;

        if (a < b)
        {
            uint128_t k = a;
            a = b;
            b = k;
            t = ((a & b & 3) == 3) ? -t : t;
            c = (b >> 2) ^ (b >> 1);
        }
        a -= b;

        if (a < 3)
        {
            if (a == 0)
            {
                return (b == 1) ? t : 0;
            }
            if (a == 1)
            {
                return t;
            }
            return (c & 1) ? -t : t;
        }
    }

    // make b odd
    if ((b & 1) == 0)
    {
        if ((a & 1) == 0)
        {
            return 0;
        }
        v = uint128_tzcnt(b);
        if (v & 1)
        {
            t = (((a >> 1) ^ (a >> 2)) & 1) ? -t : t;
        }
        b >>= v;
    }

    if (b == 1)
    {
        return t;
    }

    int r = uint64_jacobi((uint64_t)a, (uint64_t)b);
    return (t < 0) ? -r : r;
}

// integer square root (rounded down)
// assume x < 2^63
static uint64_t uint64_isqrt(uint64_t x)
{
    // Avoid divide by zero
    if (x < 2)
    {
        return x;
    }
    // This code is based on the fact that
    // sqrt(x) == x^1/2 == 2^(log2(x)/2)
    uint64_t log2x = uint64_log_2(x);
    uint64_t log2y = log2x / 2;
    uint64_t y = 1ul << log2y;
    uint64_t y_squared = 1ul << (2 * log2y);
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
    // The estimate may still be too high
    y -= (-sqr_diff / (2 * y)) + 1;
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    // The estimate may still be too high
    y -= (-sqr_diff / (2 * y)) + 1;
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    // The estimate may still be 1 too high
    return sqr_diff < 0 ? y - 1 : y;
}

// return smallest factor of n < 2^64 , exact if n < 157*157, return 1 if none is found.
static uint64_t uint64_small_factor(uint64_t n)
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

// modular exponentiation 2^e mod m
// assume e > 0, m > 0
// valid if e < 2^64, m < 2^64
static uint64_t uint64_pow2_mod(uint64_t e, uint64_t m)
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

// modular exponentiation 2^e mod m
// assume e > 0, m > 0
// valid if e < 2^128, m < 2^128
static uint128_t uint128_pow2_mod(uint128_t e, uint128_t m)
{
    uint64_t n = uint128_log_2(e);
    uint64_t s = (n >= 5) ? 5 : n;
    n -= s;
    uint128_t mask = e >> n;
    uint128_t result = uint128_shift_mod(1ull, mask, m);
    while (n >= 6)
    {
        n -= 6;
        result = uint128_square_mod(result, m);
        result = uint128_square_mod(result, m);
        result = uint128_square_mod(result, m);
        result = uint128_square_mod(result, m);
        result = uint128_square_mod(result, m);
        result = uint128_square_mod(result, m);
        mask = (e >> n) & 0x3f;
        result = uint128_shift_mod(result, mask, m);
    }
    while (n--)
    {
        result = uint128_square_mod(result, m);
        if ((e >> n) & 1)
        {
            result <<= 1;
            result -= (result >= m) ? m : 0;
        }
    }
    return result;
}

// exponentiation a^e
static uint64_t pow(uint64_t a, uint64_t e)
{
    if (e < 3)
    {
        switch (e)
        {
        case 0:
            return 1;
        case 1:
            return a;
        case 2:
            return a * a;
        }
    }

    uint64_t n = uint64_log_2(e);
    uint64_t result = a;
    while (n--)
    {
        result *= result;
        if ((e >> n) & 1)
            result *= a;
    }
    return result;
}

// modular exponentiation a^e mod m
// assume e > 0, m > 0
static uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t n = uint64_log_2(e);
    uint64_t result = a;
    while (n--)
    {
        result = square_mod(result, m);
        if ((e >> n) & 1)
            result = mul_mod(result, a, m);
    }
    return result;
}

static uint128_t uint128_pow_mod(uint128_t a, uint128_t e, uint128_t m)
{
    uint64_t n = uint128_log_2(e);
    uint128_t result = a;
    while (n--)
    {
        result = uint128_square_mod(result, m);
        if ((e >> n) & 1)
            result = uint128_mul_mod(result, a, m);
    }
    return result;
}

struct barrett_t
{
    uint64_t m;    // modulus n bits
    uint64_t q;    // quotient 2^(3n/2) / m
    uint64_t r;    // remainder 2^(3n/2) % m
    uint64_t n;    // modulus size
    uint64_t n2;   // 1/2 modulus size
    uint64_t n32;  // 3/2 modulus size
    uint64_t n321; // n32+1
};

static void barrett_precompute(struct barrett_t *p, uint64_t m)
{
    // precompute a variant of Barrett reduction
    p->m = m;
    p->n = 1 + uint64_log_2(m);
    if (p->n <= 29)
    {
        // Barrett precomputations
        p->n2 = (p->n + 1) << 1;
        p->n32 = p->n + (p->n >> 1);
        p->n321 = 0;
        p->r = ((uint128_t)1 << p->n2) % m;
        p->q = (1ull << p->n2) / m;
        return;
    }

    if (p->n < 42)
    {
        // modified Barrett precomputations
        p->n2 = p->n << 1;
        p->n32 = p->n + (p->n >> 1);
        p->n321 = p->n32 + 1;
        p->r = ((uint128_t)1 << p->n32) % m;
        p->q = ((uint128_t)1 << p->n321) / m;
        return;
    }

    // no optimization
    p->n2 = 0;
    p->n32 = 0;
    p->n321 = 0;
    p->q = 0;
    p->r = 0;
}

static uint64_t barrett_mul_mod(uint64_t u, uint64_t v, const struct barrett_t &bt)
{
    if (bt.n < 30)
    {
        // for modulus m up to 30 bits (3 multiplications)
        // assume u, v <= 2 * m
        // makes r <= 2 * m
        uint64_t r = u * v;
        uint64_t e = ((uint128_t)r * bt.q) >> bt.n2;
        r -= e * bt.m; // barrett subtraction without underflow
        return r;
    }

    if (bt.n < 42)
    {
        // for modulus up to 42 bits  (4 multiplications)
        // assume u, v <= 2 * m
        // makes r <= 2 * m
        uint128_t t = (uint128_t)u * v;                // up to 86 bits
        uint64_t t_lo = t & ((1ull << bt.n32) - 1);    // up to 63 bits
        uint64_t t_hi = t >> bt.n32;                   // up to 23 bits
        uint64_t b = t_lo + t_hi * bt.r;               // up to 64 bits
        uint64_t e = ((uint128_t)bt.q * b) >> bt.n321; // up to 85 bits and down to 22 bits
        uint64_t r = b - e * bt.m;                     // barrett subtraction without underflow r < 3*m
        return r - ((r >= bt.m) ? bt.m : 0);           // r <= 2 * m
    }

    // fall-back
    // no optimisation (2 asm instructions, including a slow long division)
    // assume u, v <= 2 * m
    // makes r <= 2 * m
    return mul_mod(u, v, bt.m); // assume u*v < (2^64-1)*m , i.e. worst case m < 63 bits
}

// modular reduction
static uint64_t barrett_long_mod(uint128_t u, const struct barrett_t &bt)
{
    if (bt.n < 30)
    {
        // makes u <= (2m)^2 approx
        while (u >> bt.n2)
        {
            uint64_t u_lo = u & ((1ull << bt.n2) - 1);
            uint64_t u_hi = u >> bt.n2;
            u = u_hi * bt.r + u_lo;
        }
        // makes r <= 2 * m
        uint64_t r = u;
        uint64_t e = ((uint128_t)r * bt.q) >> bt.n2;
        r -= e * bt.m; // barrett subtraction without underflow
        return r;
    }

    if (bt.n < 42)
    {
        // makes u <= (2*m) ^ 3/2 approx
        while (u >> bt.n32)
        {
            uint64_t u_lo = u & ((1ull << bt.n32) - 1);
            uint128_t u_hi = u >> bt.n32;
            u = u_hi * bt.r + u_lo;
        }
        // makes r <= 2 * m
        uint64_t r = (uint64_t)u;
        uint64_t e = ((uint128_t)bt.q * r) >> bt.n321;
        r -= e * bt.m; // barrett subtraction without underflow
        return r - ((r >= bt.m) ? bt.m : 0);
    }

    // fall-back
    // no optimisation (2 asm instructions, including a slow long division)
    // assume u, v <= 2 * m
    // makes r <= 2 * m
    return uint128_long_mod(u, bt.m); // assume u*v < (2^64-1)*m , i.e. worst case m < 63 bits
}

// modular exponentiation a^e mod m with precomputations
// intermediate numbers are less than (3 * m)
static uint64_t barrett_pow_mod(uint64_t a, uint64_t e, const barrett_t &bt)
{
    if (e == 0)
        return 1;
    uint64_t bits = uint64_log_2(e);
    uint64_t result = a;
    while (bits--)
    {
        result = barrett_mul_mod(result, result, bt);
        if ((e >> bits) & 1)
        {
            result = barrett_mul_mod(result, a, bt);
        }
    }
    // final reduction, at most 3 subtractions.
    while (result >= bt.m)
    {
        result -= bt.m;
    }
    return result;
}

// MR strong test
static bool uint64_witness(uint64_t n, uint64_t s, uint64_t d, uint64_t a)
{
    uint64_t x, y;
    if (n == a)
        return true;

    if (a == 2)
    {
        x = uint64_pow2_mod(d, n);
    }
    else
    {
        x = pow_mod(a, d, n);
    }

    while (s)
    {
        y = square_mod(x, n);
        if (y == 1 && x != 1 && x != n - 1)
        {
            return false;
        }
        x = y;
        --s;
    }
    if (x != 1)
    {
        return false;
    }
    return true;
}

static bool uint128_witness(uint128_t n, uint64_t s, uint128_t d, uint128_t a)
{
    uint128_t x, y;
    if (n == a)
        return true;

    if (a == 2)
    {
        x = uint128_pow2_mod(d, n);
    }
    else
    {
        x = uint128_pow_mod(a, d, n);
    }

    while (s)
    {
        y = uint128_square_mod(x, n);
        if (y == 1 && x != 1 && x != n - 1)
        {
            return false;
        }
        x = y;
        --s;
    }
    if (x != 1)
    {
        return false;
    }
    return true;
}

// deterministic primality test for n < 2^64.
// Assume that small factors are already processed, assume n > 2
bool uint64_is_prime_mr(uint64_t n)
{
    uint64_t d = n / 2;
    uint64_t s = uint64_tzcnt(d);
    d >>= s++;

    if (n < 1373653)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3);
    if (n < 9080191)
        return uint64_witness(n, s, d, 31) && uint64_witness(n, s, d, 73);
    if (n < 4759123141)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 61);
    if (n < 1122004669633)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 13) && uint64_witness(n, s, d, 23) &&
               uint64_witness(n, s, d, 1662803);
    if (n < 2152302898747)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3) && uint64_witness(n, s, d, 5) &&
               uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 11);
    if (n < 3474749660383)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3) && uint64_witness(n, s, d, 5) &&
               uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 11) && uint64_witness(n, s, d, 13);
    if (n < 341550071728321)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3) && uint64_witness(n, s, d, 5) &&
               uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 11) && uint64_witness(n, s, d, 13) &&
               uint64_witness(n, s, d, 17);
    if (n < 3825123056546413051)
        return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3) && uint64_witness(n, s, d, 5) &&
               uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 11) && uint64_witness(n, s, d, 13) &&
               uint64_witness(n, s, d, 17) && uint64_witness(n, s, d, 19) && uint64_witness(n, s, d, 23);
    // n < 318665857834031151167461
    return uint64_witness(n, s, d, 2) && uint64_witness(n, s, d, 3) && uint64_witness(n, s, d, 5) &&
           uint64_witness(n, s, d, 7) && uint64_witness(n, s, d, 11) && uint64_witness(n, s, d, 13) &&
           uint64_witness(n, s, d, 17) && uint64_witness(n, s, d, 19) && uint64_witness(n, s, d, 23) &&
           uint64_witness(n, s, d, 29) && uint64_witness(n, s, d, 31) && uint64_witness(n, s, d, 37);
}

// primality test for n < 2^64.
// Assume that small factors and small primes are already processed, assume n >= 5
// Assume n is not a perfect square

static bool uint64_lucas_nist(uint64_t n)
{
    int64_t d = 5;
    int j;
    int sgn = 1;

    // process the sequence 5, -7, 9, -11, 13, -15 .....
    while (1)
    {
        if (d != n)
        {
            j = int64_kronecker(sgn * d, n);
            if (j == 0)
            {
                return false; // composite
            }
            if (j == -1)
            {
                break; // quadratic non-residue
            }
        }
        d += 2;
        sgn = -sgn;
    }

    uint64_t e = n + 1;
    uint64_t bits = uint64_log_2(e);
    uint64_t D = sgn < 0 ? n - d % n : d % n;
    uint64_t U = 1;
    uint64_t V = 1;
    uint128_t Ut, Vt;

    while (bits--)
    {
        /* Double */
        Vt = mul_mod(D, U, n);
        Vt *= U;
        Vt += square_mod(V, n);
        Vt += (Vt & 1) ? n : 0;
        U = mul_mod(U, V, n);
        V = uint128_long_mod(Vt >> 1, n);

        if ((e >> bits) & 1)
        {
            /* Add */
            Ut = U;
            Ut += V;
            Ut += (Ut & 1) ? n : 0;
            Vt = D;
            Vt *= U;
            Vt += V;
            Vt += (Vt & 1) ? n : 0;
            U = uint128_long_mod(Ut >> 1, n);
            V = uint128_long_mod(Vt >> 1, n);
        }
    }

    return (U == 0);
}

// primality test for n < 2^64.
// Assume that small factors and small primes are already processed, assume n >= 5
static bool uint64_lucas_bpsw(uint64_t n)
{
    int j;
    uint64_t P = 3;
    uint64_t D = 5;

    // process the sequence 5, 12, 21, 32, .... (require n > 3)
    while (1)
    {
        D = square_add_mod(P, n - 4, n);
        if (D)
        {
            j = uint64_jacobi(D, n);
            if (j == 0)
            {
                return false; // composite
            }
            if (j == -1)
            {
                break; // quadratic non-residue
            }
        }
        P += 1;
    }

    uint64_t e = n + 1;
    uint64_t d = e >> 1, s = 1;
    while ((d & 1) == 0)
    {
        d >>= 1;
        s += 1;
    }

    uint64_t bits = 1 + uint64_log_2(d);
    uint64_t U = 0;
    uint64_t V = 2;
    uint128_t Ut, Vt;

    while (bits--)
    {
        /* Double */
        U = mul_mod(U, V, n);
        V = square_add_mod(V, n - 2, n);
        if ((d >> bits) & 1)
        {
            /* Add */
            Ut = P;
            Ut *= U;
            Ut += V;
            Vt = D;
            Vt *= U;
            Vt += mul_mod(P, V, n);
            Ut += (Ut & 1) ? n : 0;
            Vt += (Vt & 1) ? n : 0;
            U = uint128_long_mod(Ut >> 1, n);
            V = uint128_long_mod(Vt >> 1, n);
        }
    }
    bool b = (U == 0 && (V == 2 || V == n - 2));
    if (b)
        return true;

    while (s--)
    {
        if (V == 0)
            return true;
        V = square_add_mod(V, n - 2, n);
    }

    return false;
}

static bool uint128_lucas_nist(uint128_t n)
{
    int128_t d = 5;
    int j;
    int sgn = 1;

    // process the sequence 5, -7, 9, -11, 13, -15 .....
    while (1)
    {
        if (d != n)
        {
            j = int128_kronecker(d * sgn, n);
            if (j == 0)
            {
                return false; // composite
            }
            if (j == -1)
            {
                break; // quadratic non-residue
            }
        }
        d += 2;
        sgn = -sgn;
    }

    uint128_t D = sgn < 0 ? n - d % n : d % n;
    uint128_t e = n + 1;
    uint64_t bits = uint128_log_2(e);
    uint256_t Ut, Vt;
    uint128_t U = 1;
    uint128_t V = 1;

    while (bits--)
    {
        /* Double */
        Vt = uint128_mul_mod(D, U, n);
        Vt *= U;
        Vt += uint128_square_mod(V, n);
        Vt += (Vt & 1) ? n : 0;
        U = uint128_mul_mod(U, V, n);
        V = uint256_long_mod(Vt >> 1, n);

        if ((e >> bits) & 1)
        {
            /* Add */
            Ut = U;
            Ut += V;
            Ut += (Ut & 1) ? n : 0;
            Vt = D;
            Vt *= U;
            Vt += V;
            Vt += (Vt & 1) ? n : 0;
            U = uint256_long_mod(Ut >> 1, n);
            V = uint256_long_mod(Vt >> 1, n);
        }
    }

    return (U == 0);
}

static bool uint64_is_prime_nist(uint64_t n)
{
    if (n < 5)
    {
        return n == 2 || n == 3;
    }
    if ((n & 1) == 0)
    {
        return false;
    }

    uint64_t d = n >> 1;
    uint64_t s = uint64_tzcnt(d);
    d >>= s++;

    bool b = uint64_witness(n, s, d, 2);
    if (b != true)
    {
        return false; // composite
    }
    b = uint64_is_perfect_square(n);
    if (b == true)
    {
        return false; // composite
    }
    b = uint64_lucas_nist(n);
    if (b != true)
    {
        return false; // composite
    }
    // really prime, proven to 2^64
    return true;
}

static bool uint64_is_prime_bpsw(uint64_t n)
{
    if (n < 5)
    {
        return n == 2 || n == 3;
    }
    if ((n & 1) == 0)
    {
        return false;
    }

    uint64_t d = n >> 1;
    uint64_t s = uint64_tzcnt(d);
    d >>= s++;

    bool b = uint64_witness(n, s, d, 2);
    if (b != true)
    {
        return false; // composite
    }
    b = uint64_is_perfect_square(n);
    if (b == true)
    {
        return false; // composite
    }
    b = uint64_lucas_bpsw(n);
    if (b != true)
    {
        return false; // composite
    }
    // really prime, proven to 2^64
    return true;
}

static bool uint128_is_prime_nist(uint128_t n)
{
    if (n < 5)
    {
        return n == 2 || n == 3;
    }
    if ((n & 1) == 0)
    {
        return false;
    }

    uint128_t d = n >> 1;
    uint64_t s = uint128_tzcnt(d);
    d >>= s++;

    bool b = uint128_witness(n, s, d, 2);
    if (b != true)
    {
        return false; // composite
    }
    b = uint128_is_perfect_square(n);
    if (b == true)
    {
        return false; // composite
    }
    return true;
    b = uint128_lucas_nist(n);
    if (b != true)
    {
        return false; // composite
    }
    // really prime, proven to 2^64, very likely to 2^128
    return true;
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

// Solve:
//        x = a (mod m)
//        x = b (mod n)
//
// assuming gcd(m, n) = 1
//
// simple implementation without overflow handling
// (m * n) < 64 bits
// (m*k+a) / (m * n) < 2^64
static uint64_t uint64_crt(uint64_t a, uint64_t m, uint64_t b, uint64_t n)
{
    uint64_t k, x;
    if (b > a)
    {
        k = (b - a) % n;
    }
    else
    {
        k = n - ((a - b) % n);
    }
    // k = ((a-b)/m ) % n
    k = mul_mod(k, uint64_mod_inv(m, n), n);
    // x = (m * k + a) % (m * n)
    x = mul_add_mod(m, k, a, m * n);

    return x;
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

static bool uint64_is_perfect_square(uint64_t a)
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

static bool uint128_is_perfect_square(uint128_t a)
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
    uint64_t m;
    uint128_t c;
    // binary search (1 more bit of square root per iteration)
    uint64_t r = (uint64_t)d;
    uint64_t l = (uint64_t)dl;
    uint64_t h = (uint64_t)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m;
        c *= m;
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
    c = r;
    c *= r; // check perfect square
    return (c == a);
}

static bool uint64_is_perfect_cube(uint64_t a)
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

static bool uint64_is_perfect_sursolid(uint64_t a)
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

// detect a perfect power
// slow algorithm, assume n < 2^64/3
static bool uint64_is_perfect_power(uint64_t n)
{
    if (n < 4)
    {
        return n == 1;
    }
    uint64_t l2 = uint64_log_2(n);
    if ((1ull << l2) == n)
    {
        // perfect power of 2;
        return true;
    }
    uint64_t exponent = l2;
    while (--exponent > 1)
    {
        uint64_t lo = 1 << (l2 / exponent); // crude underestimate
        uint64_t hi = 2 * lo + 1;           // crude overestimate
                                            // binary search of the solution, lo <= root < hi
        while (hi - lo > 1)
        {
            uint64_t mid = (lo + hi) >> 1;
            uint64_t pmid = pow(mid, exponent);
            if (pmid == n)
            {
                // luckily found
                return true;
            }
            if (pmid < n)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        // verify the root
        if (pow(lo, exponent) == n)
        {
            return true;
        }
    }
    // no root found
    return false;
}

// https://en.wikipedia.org/wiki/Shanks%27s_square_forms_factorization
static uint64_t uint64_sqfof_factor(uint64_t n)
{
    static uint64_t ks[] = {1,     3,      5,      7,         11,         3 * 5,      3 * 7,      3 * 11,
                            5 * 7, 5 * 11, 7 * 11, 3 * 5 * 7, 3 * 5 * 11, 3 * 7 * 11, 5 * 7 * 11, 3 * 5 * 7 * 11,
                            0};
    static uint64_t max_n[] = {18446744073709551615ul,
                               6148914691236517205ul,
                               3689348814741910323ul,
                               2635249153387078802ul,
                               1676976733973595601ul,
                               1229782938247303441ul,
                               878416384462359600ul,
                               558992244657865200ul,
                               527049830677415760ul,
                               335395346794719120ul,
                               239568104853370800ul,
                               175683276892471920ul,
                               111798448931573040ul,
                               79856034951123600ul,
                               47913620970674160ul,
                               15971206990224720ul,
                               0ul};
    uint64_t Pi, P0, P1, Q0, Q1, Q2, b, q, b0, B, i, k, s, iks, g;
    if (n % 2 == 0)
        return 2;
    s = uint64_isqrt(n);
    if (s * s == n)
        return s; // perfect square, factor is found

    B = 3 * 2 * uint64_isqrt(2 * s);
    iks = 0;
    while (n < max_n[iks])
    {
        k = ks[iks++];
        Pi = uint64_isqrt(k * n);
        P0 = Pi;
        Q0 = 1;
        Q1 = k * n - P0 * P0;
        for (i = 2; i < B; i++)
        {
            b = (Pi + P0) / Q1;
            P1 = b * Q1 - P0;
            Q2 = Q0 + b * (P0 - P1);
            if (i % 2 == 0 && uint64_is_perfect_square(Q2))
                break;
            P0 = P1;
            Q0 = Q1;
            Q1 = Q2;
        }
        if (i == B)
            continue;

        q = uint64_isqrt(Q2);
        b0 = (Pi - P1) / q;
        P0 = b0 * q + P1;
        Q0 = q;
        Q1 = (k * n - P0 * P0) / Q0;
        while (1)
        {
            b = (Pi + P0) / Q1;
            P1 = b * Q1 - P0;
            if (P0 == P1)
                break;
            Q2 = Q0 + b * (P0 - P1);
            P0 = P1;
            Q0 = Q1;
            Q1 = Q2;
        }

        g = uint64_gcd(n, P1);
        if (g > 1 && g < n)
        {
            return g; // a factor is found
        }
    }
    return 1; // n is prime, or cannot be factored
}

// pollard-rho factorization with brent variant
static uint64_t uint64_brent_pollard_factor(uint64_t n)
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

static bool factor_sort(const factor_t &f, const factor_t &g)
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

static void uint64_add_prime_factor(factor_v &primes, uint64_t p)
{
    // avoid storing duplicate prime factors
    reverse_iterator_v v = find_if(primes.rbegin(), primes.rend(), factor_find_t(p));
    if (v == primes.rend())
    {
        // not found, add the factor
        factor_t f;
        f.prime = p;
        f.count = 1;
        primes.push_back(f);
    }
    else
    {
        v->count += 1;
    }
}

// search for large factors, assume sieving already done up to factor 151
// input n can be prime or composite , but has no factor less than 157
static void uint64_large_factors(factor_v &primes, uint64_t n)
{
    uint64_t m;
    vector<uint64_t> factors;
    factors.push_back(n);

    do
    {
        m = factors.back();
        factors.pop_back();

        if (m == 1)
            continue;

        if (m < 157 * 157 || uint64_is_prime_bpsw(m))
        {
            uint64_add_prime_factor(primes, m);
        }
        else
        {
            uint64_t factor;
            // m is not prime,
            // get more prime and composite factors from sqfof method O(n^1/4) (might fail and return 1)
            factor = uint64_sqfof_factor(m);
            if (factor == 1)
            {
                // get more prime and composite factors from pollard-rho method O(smallest factor^1/2)   <= O(n^1/4)
                // which returns only when a factor is found. Unfortunately, random parameters makes it hasardeous,
                // and it could take a long, long time to run to completion.
                factor = uint64_brent_pollard_factor(m);
            }
            factors.push_back(m / factor);
            factors.push_back(factor);
        }
    } while (factors.size());
}

static uint64_t uint64_smallest_factor(uint64_t m)
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

static void uint64_all_factors(factor_v &factors, uint64_t m)
{
    // ---------------------------------------------------------------
    // first search all prime factors < 157
    // ---------------------------------------------------------------
    uint64_t factor = uint64_small_factor(m);
    while (factor != 1)
    {
        uint64_add_prime_factor(factors, factor);
        m /= factor;
        factor = uint64_small_factor(m);
    }
    if (m < 157 * 157)
    {
        // all factors < 157 have been removed and m is prime
        uint64_add_prime_factor(factors, m);
    }
    else
    {
        // get more factors in any order
        uint64_large_factors(factors, m);
    }

    // order the factors from the vector (more or less already ordered)
    sort(factors.begin(), factors.end(), factor_sort);
}

static bool is_perfect_power(const factor_v &factors)
{
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

static bool inline is_squarefree(const factor_v &factors)
{
    for (uint64_t i = 0; i < factors.size(); i++)
    {
        if (factors[i].count > 1)
            return false;
    }
    // no square detected in factors
    return true;
}

static bool inline is_perfect_prime_power(const factor_v &factors)
{
    // a perfect power with 1 prime unique factor
    return factors.size() == 1 && factors[0].count > 1;
}

static bool inline is_semiprime(const factor_v &factors)
{
    // a semiprime (exactly 2 different proper factors)
    return factors.size() == 2 && factors[0].count == 1 && factors[1].count == 1;
}

static bool inline is_sphenic(const factor_v &factors)
{
    // a sphenic (exactly 3 different proper factors)
    return factors.size() == 3 && factors[0].count == 1 && factors[1].count == 1 && factors[2].count == 1;
}

static bool inline is_prime(const factor_v &factors)
{
    // 1 unique factor
    return factors.size() == 1 && factors[0].count == 1;
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

static unsigned int128_sprint(char *ptr, int128_t x)
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
