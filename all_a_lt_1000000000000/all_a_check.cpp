// -----------------------------------------------------------------------
// verify details of the cubic test
//
// {gettime();forstep(n=25,1000000000000,2,
// if(n%10000==1,print([n,gettime()]));
// if(n%3!=0&&!ispseudoprime(n),p=factor(n)[1,1];Q=n/p;g=gcd(p^3-1,Q^3-1);R=(n-1)%g;
//  if(R>2,A=(n^2+n+1)%g;
//   if(A>3,
//    for(k=1,(n+1)/2,
//     if(k%3!=2,a=(7+k*(k-1))%n;
//      if(Mod(a,n)^R==1,B=Mod(Mod(x,n),x^3-a*x-a)^R;
//       if(B^2+B+1==-x^2+x+a.print([n,p,Q,g,R,A]);break(7)))))))));}
//
// -----------------------------------------------------------------------

#include <algorithm>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <vector>

#include <x86intrin.h>

using namespace std;

typedef unsigned __int128 uint128_t;
typedef signed __int128 int128_t;

// randomness generation with period 2^64
uint64_t rnd(void)
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
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(u), "1"(v), "r"(n));
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
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    return u % n;
#endif
}

static inline uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
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
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
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

// pollard-rho factorization with brent variant
uint64_t brent_pollard_factor(uint64_t n)
{
    uint64_t i, x, ys, k;
    uint64_t m = 1000;
    uint64_t a = 2 + rnd() % (n - 4);
    uint64_t y = 1 + rnd() % (n - 2);
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

void uint64_large_factors(vector<uint64_t> &primes, uint64_t n)
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
            i = primes.size();
            bool found = false;
            while (i--)
            {
                if (primes[i] == m)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // push the prime factor
                primes.push_back(m);
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
    vector<uint64_t> factors;
    uint64_large_factors(factors, m);
    // get the smallest factor from the vector
    sort(factors.begin(), factors.end());
    return factors[0];
}

// compute Mod(Mod(x, n), x^3 -a*x -a)^e
// assume n is less than 60 or 61 bits
// assume a < n
// if e odd, assume s == 0 and u == 0 and t == 1
static void cubic_exponentiate(uint64_t &s, uint64_t &t, uint64_t &u, uint64_t e, uint64_t n, uint64_t a)
{
    int bit = uint64_log_2(e);
    uint64_t tmp;
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;
    while (bit--)
    {
        // start Square
        tmp = square_mod(s, n);
        s2 = (uint128_t)tmp * a;
        t2 = (uint128_t)t * t;
        u2 = (uint128_t)u * u;
        tmp = mul_mod(s, t, n);
        st = (uint128_t)tmp * a;
        tu = (uint128_t)t * u;
        us = (uint128_t)u * s;
        st <<= 1;
        tu <<= 1;
        us <<= 1;
        if (e & (1ull << bit))
        {
            // finish Square and multiply by x
            us = s2 + us + t2;
            tmp = uint128_long_mod(us, n);
            uu = (uint128_t)tmp * a;
            ss = s2 + st + tu;
            tt = uu + u2 + st;
        }
        else
        {
            // finish Square
            ss = s2 + us + t2;
            tt = s2 + st + tu;
            uu = st + u2;
        }
        s = uint128_long_mod(ss, n);
        t = uint128_long_mod(tt, n);
        u = uint128_long_mod(uu, n);
    }
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

// read file, cautious about corruped or truncated files
bool read_file(char *fn, uint64_t *p, uint64_t *Q, uint64_t *n, uint64_t *s)
{
    char line[100];
    bool p_found = false;
    bool q_found = false;
    bool n_found = false;
    bool s_found = false;
    bool valid = true;
    uint64_t p_temp = -1, q_temp = -1, n_temp = -1, s_temp = -1;
    FILE *f = fopen(fn, "rt");
    if (f)
    {
        memset(line, 0, 100);
        while (fgets(line, 100, f))
        {
            // trim the white-spaces at end of line (noisy stuff from hand-written files)
            int l = strlen(line);
            while (l > 0 && isspace(line[l - 1]))
            {
                l--;
            }
            line[l] = 0;

            if (*line == 0 || *line == '#')
            {
                // empty line, or comment, ignore
                continue;
            }
            if (!memcmp(line, "p=", 2))
            {
                p_temp = strtoull(&line[2], 0, 0);
                p_found = true;
                continue;
            }
            if (!memcmp(line, "Q=", 2))
            {
                q_temp = strtoull(&line[2], 0, 0);
                q_found = true;
                continue;
            }
            if (!memcmp(line, "n=", 2))
            {
                n_temp = strtoull(&line[2], 0, 0);
                n_found = true;
                continue;
            }
            if (!memcmp(line, "s=", 2))
            {
                s_temp = strtoull(&line[2], 0, 0);
                s_found = true;
                continue;
            }

            // not our file, give up
            valid = false;
            break;
        }
        fclose(f);
    }
    else
    {
        valid = false;
    }

    if (valid && p_found == true && q_found == true && n_found == true && s_found == true)
    {
        // file is complete, without missing field
        *p = p_temp;
        *Q = q_temp;
        *n = n_temp;
        *s = s_temp;
        return true;
    }
    return false;
}

// write file, cautious about worst cases in the life of a computer
bool write_file(char *fn, char *bak, uint64_t p, uint64_t Q, uint64_t n, uint64_t s)
{
    uint64_t p_temp = -1, q_temp = -1, n_temp = -1, s_temp = -1;
    // first verify the current file is a valid one
    if (read_file(fn, &p_temp, &q_temp, &n_temp, &s_temp))
    {
        // Create a backup file from the file just validated
        // Unix : atomic file rename.
        // In case of crash, operation is done, or not done at all, and worst case is to have duplicate files.
        rename(fn, bak);
    }

    // create or overwrite the file with new values,
    // and in case of crash and new file is incomplete or corrupted, there is a backup just done, pfew.
    FILE *f = fopen(fn, "wt");
    if (f)
    {
        // write current time
        char buf[26];
        struct tm time_m;
        time_t time_p;
        time(&time_p);
        localtime_r(&time_p, &time_m);
        asctime_r(&time_m, buf);
        fprintf(f, "# %s\n", buf);

        // write restart parameters
        fprintf(f, "p=%ld\n", p);
        fprintf(f, "Q=%ld\n", Q);
        fprintf(f, "n=%ld\n", n);
        fprintf(f, "s=%ld\n", s);
        fclose(f);

        // never paranoid enough against full disks, need to read back the file.
        // This also flushes the file in the case of remote file (read-after-write security on NFS ....)
        p_temp = -1;
        q_temp = -1;
        n_temp = -1;
        s_temp = -1;
        if (read_file(fn, &p_temp, &q_temp, &n_temp, &s_temp))
        {
            // return true if the file is correctly written
            return (p_temp == p && q_temp == Q && n_temp == n && s_temp == s);
        }
    }
    return false;
}

int main(int argc, char **argv)
{
    uint64_t p_start = 5;
    uint64_t Q_start = 2;
    uint64_t n_max = 1000000000000ul;
    char temp_filename1[100] = "cubic_all_a.txt\0";
    char temp_filename2[100] = "cubic_all_a.bak\0";
    bool use_temp_files = true;
    uint64_t interval_seconds = 600; // 10 minutes

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-p"))
        {
            p_start = strtoull(argv[++i], 0, 0);
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-Q"))
        {
            Q_start = strtoull(argv[++i], 0, 0);
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-n"))
        {
            n_max = strtoull(argv[++i], 0, 0);
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-s"))
        {
            interval_seconds = strtoull(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-f"))
        {
            if (read_file(argv[++i], &p_start, &Q_start, &n_max, &interval_seconds))
            {
                printf("Start from values p=%ld Q=%ld, to n=%ld\n", p_start, Q_start, n_max);
            }
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-t"))
        {
            sprintf(temp_filename1, "%s.txt", argv[++i]);
            sprintf(temp_filename2, "%s.txt", argv[i]);
            printf("Temp file names are %s %s\n", temp_filename1, temp_filename2);
            continue;
        }
        printf("-p ddd : starting prime p, default is %ld\n", p_start);
        printf("-Q ddd : starting cofactor Q, default is %ld\n", Q_start);
        printf("-n ddd : maximum value for n = p*Q, default is %ld\n", n_max);
        printf("-s ddd : interval in seconds between restart file updates, default is %ld\n", interval_seconds);
        printf("-f aaa : alternative filename to restart at some (p,Q,n_max)\n");
        printf("\n");
        printf("-t aaa : temp filename prefix %s %s\n", temp_filename1, temp_filename2);
        printf("\n");
        printf("n = p * Q\n");
        printf("verify cubic test for all k < (n+1)/2\n");
        printf("NO composite number can pass the test and be found to be prime ! Never, never, never .....\n");
        printf("\n");
        exit(1);
    }

    if (use_temp_files)
    {
        if (read_file(temp_filename1, &p_start, &Q_start, &n_max, &interval_seconds))
        {
            printf("Start from values p=%ld Q=%ld, to n=%ld\n", p_start, Q_start, n_max);
        }
        else
        {
            if (read_file(temp_filename2, &p_start, &Q_start, &n_max, &interval_seconds))
            {
                printf("Start from values p=%ld Q=%ld, to n=%ld\n", p_start, Q_start, n_max);
            }
        }
    }

    // round p to next prime
    p_start = (p_start < 5) ? 5 : p_start | 1;
    while (!uint64_is_prime(p_start))
    {
        p_start += 2;
    }

    // round Q > 1 to next odd, non multiple of 3
    Q_start = (Q_start < 5) ? 5 : Q_start | 1;
    while (Q_start % 3 == 0)
    {
        Q_start += 2;
    }

    printf("Start from values p=%ld Q=%ld, to n=%ld\n", p_start, Q_start, n_max);

    time_t t0 = time(NULL);

    uint64_t dp = p_start % 6 == 1 ? 4 : 2;
    uint64_t dq = Q_start % 6 == 1 ? 4 : 2;
    uint64_t p = p_start;
    uint64_t p_last = uint64_isqrt(n_max) + 1;
    uint64_t Q = Q_start;
    uint64_t display = 0;

    while (p <= p_last)
    {
        uint64_t Q_last = (n_max + p - 1) / p;
        while (Q < Q_last)
        {
            // p prime >= 5, Q > 1 is not a multiple of 3 nor a multiple of 2
            uint64_t s = uint64_smallest_factor(Q);
            if (s < p)
            {
                // the smallest factor of n = p*Q is s, a factor of Q.
                // this is already done earlier with smaller values of p
            }
            else
            {
                // single thread progress
                if (++display > 10000)
                {
                    printf("p = %ld, Q = %ld\n", p, Q);
                    display = 0;
                }
                uint64_t n = p * Q; // n = p * Q is a proven composite
                uint128_t p3 = (uint128_t)p * p * p;
                uint128_t Q3 = (uint128_t)Q * Q * Q;
                uint128_t g = uint128_gcd(p3 - 1, Q3 - 1);
                uint64_t R = g > (n - 1) ? (n - 1) : (n - 1) % g;
                if (R > 2)
                {
                    uint128_t A = (uint128_t)n * (n + 1) + 1;
                    A = g > A ? A : A % g;
                    if (A > 3)
                    {
                        //
                        // TODO : start a child thread from here
                        //
                        uint64_t a = 7;
                        //
                        // TODO : precompute barrett(n) to accelerate mod operations mod n
                        //
                        for (uint64_t d = 0; d < n; d += 2)
                        {
                            // a = 7 + k * (k-1) unrolled as  a = a + 2k
                            a += d;
                            a = (a >= n) ? a - n : a;

                            // verify fermat test Mod(a,n)^R == 1
                            // uint64_t ft = (a == 2) ? pow2_mod(ddR, n) : pow_mod(a, R, n);
                            if (pow_mod(a, R, n) == 1)
                            {
                                // run cubic test
                                // B = Mod(x, n)
                                uint64_t bs = 0;
                                uint64_t bt = 1;
                                uint64_t bu = 0;
                                // B = Mod(B, x^3 -ax -a)^R
                                cubic_exponentiate(bs, bt, bu, R, n, a);
                                // B2 = B
                                uint64_t bs2 = bs;
                                uint64_t bt2 = bt;
                                uint64_t bu2 = bu;
                                // B2 = Mod(B2, x^3 -ax -a)^2
                                cubic_exponentiate(bs2, bt2, bu2, 2, n, a);
                                // check B2+B+1 is NOT -x^2 + x + 1
                                bs2 = uint64_add_mod(bs2, bs, n);
                                if (bs2 == n - 1)
                                {
                                    bt2 = uint64_add_mod(bt2, bt, n);
                                    if (bt2 == 1)
                                    {
                                        bu2 = uint64_add_mod(bu2, bu + 1, n);
                                        if (bu2 == a)
                                        {
                                            // Aoutch, n = p*Q could be prime or pseudoprime
                                            char buff[256];
                                            char *ptr = buff;
                                            *ptr++ = '[';
                                            ptr += uint128_sprint(ptr, n);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, p);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, Q);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, g);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, R);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, A);
                                            *ptr++ = ']';
                                            *ptr = 0;
                                            printf("%s\n", buff);
                                            fflush(stdout);
                                            // got 1 counter-example, it is time to stop
                                            exit(1);
                                        }
                                        else
                                        {
                                            // n is composite for sure
                                        }
                                    }
                                    else
                                    {
                                        // n is composite for sure
                                    }
                                }
                                else
                                {
                                    // n is composite for sure
                                }
                            }
                        }
                    }
                }
            }

            // next Q, not a multiple of 3, not a multiple of 2
            Q += dq;
            dq = 6 - dq;

            // after 30 minutes, update the temp files.
            // In case of crash, just restarting the program without parameter would restart from last saved values
            time_t t1 = time(NULL);
            if (t1 - t0 >= interval_seconds)
            {
                if (!write_file(temp_filename1, temp_filename2, p, Q, n_max, interval_seconds))
                {
                    printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
                }
                else
                {
                    t0 = t1;
                }
            }
        }

        // next p, not a multiple of 3, not a multiple of 2, and prime
	do
	{
        p += dp;
        dp = 6 - dp;
	}
	while (uint64_small_factor(p) != 1 || !uint64_is_prime(p));

        // first Q >= p, not a multiple of 3, not a multiple of 2
        Q = p;
        dq = dp;
    }

    return (0);
}
