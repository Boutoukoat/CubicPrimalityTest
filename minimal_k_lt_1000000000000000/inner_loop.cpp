
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <x86intrin.h>

#include "inner_loop.h"
#include "lcg.h"
#include "tlv.h"

template <class T> unsigned log_2(const T &x)
{
    if (x)
    {
        return 63 - __builtin_clzll(x);
    }
    else
    {
        return 0;
    }
}

template <> unsigned log_2<uint128_t>(const uint128_t &x)
{
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    if (hi)
    {
        return 127 - __builtin_clzll(hi);
    }
    else if (lo)
    {
        return 63 - __builtin_clzll(lo);
    }
    else
    {
        return 0;
    }
}

template <> unsigned log_2<uint64_t>(const uint64_t &x)
{
    if (x)
    {
        return 63 - __builtin_clzll(x);
    }
    else
    {
        return 0;
    }
}

template <class T, class TT> static T square_mod(const T &u, const T &q)
{
    TT r = u;
    r *= u;
    r %= q;
    return (T)r;
}

template <> uint64_t square_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &q)
{
    uint64_t r, a;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(u), "r"(u));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(q));
    return r;
}

template <class T, class TT> static T mul_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r *= v;
    r %= q;
    return (T)r;
}

template <> uint64_t mul_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &v, const uint64_t &q)
{
    uint64_t r, a;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(u), "r"(v));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(q));
    return r;
}

template <class T> static T sub_mod(const T &u, const T &v, const T &q)
{
    T s = u - v;
    if (s <= u)
    {
        s %= q;
    }
    else
    {
        s = -s;
        s %= q;
        s = q - s;
    }
    return s;
}

template <class T, class TT> static T add_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r += v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T mod(const TT &u, const T &q)
{
    TT r = u;
    r %= q;
    return (T)r;
}

template <class T, class TT> static TT square(const T &u)
{
    TT r = u;
    r *= u;
    return (T)r;
}

template <class T, class TT> static TT mul(const T &u, const T &v)
{
    TT r = u;
    r *= v;
    return (T)r;
}

template <class T, class TT> static TT mul2(const T &u, const T &v)
{
    TT r = u;
    r *= v;
    r <<= 1;
    return (T)r;
}

template <class T> static T gcd(const T &x, const T &y)
{
    if (x == 0)
        return y;
    if (y == 0)
        return x;
    unsigned tu = __builtin_ctzll(x);
    unsigned tv = __builtin_ctzll(y);
    unsigned h = tu > tv ? tv : tu;
    T u = x >> tu;
    T v = y >> tv;
    while (1)
    {
        if (u > v)
        {
            T t = u;
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

template <class T> static int jacobi(const T &x, const T &y)
{

    int t = 1;
    T a = x;
    T n = y;
    T v = n % 8;
    bool c = (v == 3) || (v == 5);
    while (a)
    {
        while ((a & 1) == 0)
        {
            a = a / 2;
            if (c)
                t = -t;
        }

        if (a < n)
        {
            v = a;
            a = n;
            n = v;
            if ((a % 4) == 3 && (n % 4) == 3)
                t = -t;
            v = n % 8;
            c = (v == 3) || (v == 5);
        }

        a = (a - n) / 2;
        if (c)
            t = -t;
    }

    return (n == 1) ? t : 0;
}

template <class T, class TT> static T mod_inv(const T &x, const T &m)
{
    if (m < 3)
        return 0;
    if (x < 2)
        return x;
    T a = x, b = m, u = 1, v = 0;
    while (a != 0)
    {
        unsigned ta = __builtin_ctzll(a);
        a >>= ta;
        while (ta--)
        {
            if (u & 1)
            {
                TT tu = u;
                tu += m;
                tu >>= 1;
                u = (T)tu;
            }
            else
            {
                u >>= 1;
            }
        }
        if (a < b)
        {
            T t = a;
            T s = u;
            a = b;
            u = v;
            b = t;
            v = s;
        }
        a -= b;
        u = sub_mod<T>(u, v, m);
    }
    return v;
}

template <class T, class TT> static bool is_perfect_square(const T &a)
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
    T c, m;
    // binary search (1 more bit of square root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
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

template <class T, class TT> static bool is_perfect_cube(const T &a)
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
    T c, m;
    // binary search (1 more bit of cube root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
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

template <class T, class TT> static bool is_perfect_sursolid(const T &a)
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
    T c, m;
    // binary search (1 more bit of fifth root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
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

// a^e mod m
template <class T, class TT> static T power(const T &a, const T &e, const T &m)
{
    T n = e;
    T s = a;
    T result = 1;
    while (n)
    {
        if (n & 1)
            result = mul_mod<T, TT>(result, s, m);
        s = square_mod<T, TT>(s, m);
        n >>= 1;
    }
    return result;
}

// MR strong test
template <class T, class TT> static bool witness(const T &n, int s, const T &d, const T &a)
{
    T x, y;
    if (n == a)
        return true;
    x = power<T, TT>(a, d, n);
    while (s)
    {
        y = square_mod<T, TT>(x, n);
        if (y == 1 && x != 1 && x != n - 1)
            return false;
        x = y;
        --s;
    }
    if (y != 1)
        return false;
    return true;
}

// sieve small factors <= 151
template <class T> static bool sieve(const T &n)
{
    if (n <= 152)
    {
        // return false for small numbers which are composite for sure , without checking further.
        bool stooopid_prime_table[] = {
            true,  true,  true,  true,  false, true,  false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, false, false, true,  false, false, false, true,  false, false, false, false, false, true,
            false, false, false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false, false, false, true,  false, false, false, false, false, false,
            false, false, false, false, false, false, false, true,  false, false, false, true,  false, false, false,
            false, false, true,  false, true,  false, false, false, false, false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,  false,
            false, false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, false, false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false};
        return stooopid_prime_table[n];
    }
    if ((uint64_t)(n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return false; // divisible by 3
    if ((uint64_t)(n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return false; // divisible by 5
    if ((uint64_t)(n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return false; // divisible by 7
    if ((uint64_t)(n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return false; // divisible by 11
    if ((uint64_t)(n * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return false; // divisible by 13
    if ((uint64_t)(n * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return false; // divisible by 17
    if ((uint64_t)(n * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return false; // divisible by 19
    if ((uint64_t)(n * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return false; // divisible by 23
    if ((uint64_t)(n * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return false; // divisible by 29
    if ((uint64_t)(n * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return false; // divisible by 31
    if (n < 37 * 37)
        return true; // prime
    if ((uint64_t)(n * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return false; // divisible by 37
    if ((uint64_t)(n * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return false; // divisible by 41
    if ((uint64_t)(n * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return false; // divisible by 43
    if ((uint64_t)(n * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return false; // divisible by 47
    if ((uint64_t)(n * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return false; // divisible by 53
    if ((uint64_t)(n * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return false; // divisible by 59
    if ((uint64_t)(n * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return false; // divisible by 61
    if ((uint64_t)(n * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return false; // divisible by 67
    if ((uint64_t)(n * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return false; // divisible by 71
    if ((uint64_t)(n * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return false; // divisible by 73
    if ((uint64_t)(n * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return false; // divisible by 79
    if ((uint64_t)(n * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return false; // divisible by 83
    if ((uint64_t)(n * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return false; // divisible by 89
    if ((uint64_t)(n * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return false; // divisible by 97
    if (n < 101 * 101)
        return true; // prime
    if ((uint64_t)(n * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return false; // divisible by 101
    if ((uint64_t)(n * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return false; // divisible by 103
    if ((uint64_t)(n * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return false; // divisible by 107
    if ((uint64_t)(n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return false; // divisible by 109
    if ((uint64_t)(n * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return false; // divisible by 113
    if ((uint64_t)(n * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return false; // divisible by 127
    if ((uint64_t)(n * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return false; // divisible by 131
    if ((uint64_t)(n * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return false; // divisible by 137
    if ((uint64_t)(n * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return false; // divisible by 139
    if ((uint64_t)(n * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return false; // divisible by 149
    if ((uint64_t)(n * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return false; // divisible by 151
    return true;      // might be prime
}

template <class T, class TT> static bool isprime(const T &n)
{
    if (!sieve<T>(n))
        return false; // composite
    if (n < 157 * 157)
        return true; // prime
    T d = n / 2;
    int s = 1;
    while (!(d & 1))
    {
        d /= 2;
        ++s;
    }

    if (n < 1373653)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3);
    if (n < 9080191)
        return witness<T, TT>(n, s, d, 31) && witness<T, TT>(n, s, d, 73);
    if (n < 4759123141)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 61);
    if (n < 1122004669633)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 13) && witness<T, TT>(n, s, d, 23) &&
               witness<T, TT>(n, s, d, 1662803);
    if (n < 2152302898747)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11);
    if (n < 3474749660383)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13);
    if (n < 341550071728321)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
               witness<T, TT>(n, s, d, 17);
    if (n < 3825123056546413051)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
               witness<T, TT>(n, s, d, 17) && witness<T, TT>(n, s, d, 19) && witness<T, TT>(n, s, d, 23);
    // n < 318665857834031151167461
    return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
           witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
           witness<T, TT>(n, s, d, 17) && witness<T, TT>(n, s, d, 19) && witness<T, TT>(n, s, d, 23) &&
           witness<T, TT>(n, s, d, 29) && witness<T, TT>(n, s, d, 31) && witness<T, TT>(n, s, d, 37);
}

template <class T, class TT> void ok_exponentiate3(T &s, T &t, T &u, const T e, const T n, const T a)
{
    int bit = log_2<T>(e);
    T tmp;
    TT s2, t2, u2, st, tu, us, uu, ss, tt;
    while (bit--)
    {
        // Double
        tmp = square_mod<T, TT>(s, n);
        s2 = (TT)tmp * a;
        t2 = (TT)t * t;
        u2 = (TT)u * u;
        tmp = mul_mod<T, TT>(s, t, n);
        st = (TT)tmp * a;
        tu = (TT)t * u;
        us = (TT)u * s;
        st <<= 1;
        tu <<= 1;
        us <<= 1;
        if (e & ((T)1 << bit))
        {
            // add
            us = s2 + us + t2;
            tmp = mod<T, TT>(us, n);
            uu = (TT)tmp * a;
            ss = s2 + st + tu;
            tt = uu + u2 + st;
        }
        else
        {
            ss = s2 + us + t2;
            tt = s2 + st + tu;
            uu = st + u2;
        }
        s = mod<T, TT>(ss, n);
        t = mod<T, TT>(tt, n);
        u = mod<T, TT>(uu, n);
    }
}

//  Mod(Mod(x+t,n),x^2-a)^e
template <class T, class TT> void exponentiate2(T &s, T &t, const T e, const T n, const T a)
{
    T t0 = t;
    unsigned bit = log_2<T>(e);
    while (bit--)
    {
        T t2 = square_mod<T, TT>(t, n);
        T s2 = square_mod<T, TT>(s, n);
        TT ss = mul_mod<T, TT>(s, t, n);
        ss += ss;
        TT tt = (TT)s2 * a + t2;

        if (e & ((T)1 << bit))
        {
            TT tmp = ss * a;
            ss = t0 * ss + tt;
            tt = t0 * tt + tmp;
        }

        s = mod<T, TT>(ss, n);
        t = mod<T, TT>(tt, n);
    }
}

//  Mod(Mod(x+u,n),x^3-a*x-a)^e
template <class T, class TT> void exponentiate3(T &s, T &t, T &u, const T e, const T n, const T a)
{
    unsigned bit = log_2<T>(e);
    T u0 = u;

    while (bit--)
    {
        TT rt0 = (TT)a * square_mod<T, TT>(s, n);
        TT rt1 = (TT)a * mul_mod<T, TT>(s, t, n);
        rt1 += rt1;
        TT q0 = mul_mod<T, TT>(s, u, n);
        q0 += q0;
        q0 += square_mod<T, TT>(t, n);
        q0 += rt0;
        TT q1 = mul_mod<T, TT>(t, u, n);
        q1 += q1;
        q1 += rt0;
        q1 += rt1;
        TT q2 = square_mod<T, TT>(u, n);
        q2 += rt1;

        if (e & ((T)1 << bit))
        {
            TT rt2 = q0 * a;
            q0 = q0 * u0 + q1;
            q1 = q1 * u0 + q2 + rt2;
            q2 = q2 * u0 + rt2;
        }

        s = mod<T, TT>(q0, n);
        t = mod<T, TT>(q1, n);
        u = mod<T, TT>(q2, n);
    }
}

template <class T, class TT> static bool islnrc2prime(const T &n, int s = 0, int cid = 0)
{
    if (n < 23)
    {
        // prime for sure
        return (n == 1 || n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 || n == 19);
    }

    if (is_perfect_square<T, TT>(n))
    {
        return false; // composite
    }

    T k, a, b, c, bs, bt, g, j;
    for (k = 1;; k++)
    {
        a = (1 << k) + 1;
        b = a - 1;
        c = 2 * a - 1;
        g = power<T, TT>(a, 6, n) - 1;
        g = gcd<T>(g, n);
        if (g == n)
            continue; // try another a
        if (g > 1)
            return false; // composite
        g = 3 * b * b + a;
        g = gcd<T>(g, n);
        if (g == n)
            continue; // try another a
        if (g > 1)
            return false; // composite
        j = jacobi<T>(a, n);
        if (j == 0)
            return false; // composite
        if (j == 1)
            continue; // try another a

        if (power<T, TT>(2, n, n) != 2)
            return false; // composite;
        if (power<T, TT>(a, n, n) != a)
            return false; // composite;
        if (power<T, TT>(b, n, n) != b)
            return false; // composite;
        if (power<T, TT>(c, n, n) != c)
            return false; // composite;

        bs = 1;
        bt = b;
        exponentiate2<T, TT>(bs, bt, n + 1, n, a);
        return (bs == 0 && bt == (b * b - a) % n);
    }
}

template <class T, class TT> static bool islnrc3prime(const T &n, int s = 0, int cid = 0)
{
    if (n < 23)
    {
        // prime for sure
        return (n == 1 || n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 || n == 19);
    }

    if (is_perfect_cube<T, TT>(n))
    {
        return false; // composite
    }

    T k, a, bs, bt, bu;
    for (k = 1;; k++)
    {
        a = 7 + k * (k - 1);
        if (!isprime<T, TT>(a))
        {
            continue; // try another a
        }
        if (power<T, TT>(mod<T, TT>(n, a), (a - 1) / 3, a) == 1)
        {
            continue; // try another a
        }
        if (a == n)
        {
            return true; // small prime
        }
        T g = gcd((2 * k - 1) * a * (2 * a - 1), n);
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
        exponentiate3<T, TT>(bs, bt, bu, n - 1, n, a);
        if (bs == 0 && bt == 0 && bu == 1)
        {
            // log this cornercase for verification purposes
            if (s | cid)
            {
                tlv_write(s, cid, TLV_B1, n);
            }
            continue; // try another a
        }
        break;
    }
    T bs2 = bs, bt2 = bt, bu2 = bu;
    exponentiate3<T, TT>(bs2, bt2, bu2, 2, n, a);
    bs = add_mod<T, TT>(bs, bs2, n);
    bt = add_mod<T, TT>(bt, bt2, n);
    bu = add_mod<T, TT>(bu, bu2, n);
    bu = add_mod<T, TT>(bu, 1, n);
    if (bs != n - 1 || bt != 1 || bu != a)
        return false; // composite
    return true;      // might be prime
}

// Mod(Mod(x, n), x^5 - a*x - 4*a)^e
template <class T, class TT> void exponentiate5(T &vm, T &wm, T &xm, T &ym, T &zm, const T e, const T n, const T a)
{
    int bit = log_2<T>(e);

    const T p = 0;
    const T q = 0;
    const T r = 0;
    const T s = a;
    const T t = 4 * a;

    const T bp = p * p + q;
    const T bq = p * q + r;
    const T br = p * r + s;
    const T bs = p * s + t;
    const T bt = p * t;

    const T cp = bp * p + bq;
    const T cq = bp * q + br;
    const T cr = bp * r + bs;
    const T cs = bp * s + bt;
    const T ct = bp * t;

    const T dp = cp * p + cq;
    const T dq = cp * q + cr;
    const T dr = cp * r + cs;
    const T ds = cp * s + ct;
    const T dt = cp * t;

    TT tmp;
    T t0, t1, t2, t3, t4;
    TT qv, qw, qx, qy, qz;
    while (bit--)
    {
        tmp = square<T, TT>(vm);
        t0 = mod<T, TT>(tmp, n);
        tmp = mul2<T, TT>(vm, wm);
        t1 = mod<T, TT>(tmp, n);
        tmp = square<T, TT>(wm) + mul2<T, TT>(vm, xm);
        t2 = mod<T, TT>(tmp, n);
        tmp = mul2<T, TT>(wm, xm) + mul2<T, TT>(vm, ym);
        t3 = mod<T, TT>(tmp, n);
        qv = mul<T, TT>(dp, t0) + mul<T, TT>(cp, t1) + mul<T, TT>(bp, t2) + mul<T, TT>(p, t3) + square<T, TT>(xm) +
             mul2<T, TT>(wm, ym) + mul2<T, TT>(vm, zm);
        qw = mul<T, TT>(dq, t0) + mul<T, TT>(cq, t1) + mul<T, TT>(bq, t2) + mul<T, TT>(q, t3) + mul2<T, TT>(xm, ym) +
             mul2<T, TT>(wm, zm);
        qx = mul<T, TT>(dr, t0) + mul<T, TT>(cr, t1) + mul<T, TT>(br, t2) + mul<T, TT>(r, t3) + square<T, TT>(ym) +
             mul2<T, TT>(xm, zm);
        qy = mul<T, TT>(ds, t0) + mul<T, TT>(cs, t1) + mul<T, TT>(bs, t2) + mul<T, TT>(s, t3) + mul2<T, TT>(ym, zm);
        qz = mul<T, TT>(dt, t0) + mul<T, TT>(ct, t1) + mul<T, TT>(bt, t2) + mul<T, TT>(t, t3) + square<T, TT>(zm);

        if (e & ((T)1 << bit))
        {
            t4 = mod<T, TT>(qv, n);
            qv = mul<T, TT>(p, t4) + qw;
            qw = mul<T, TT>(q, t4) + qx;
            qx = mul<T, TT>(r, t4) + qy;
            qy = mul<T, TT>(s, t4) + qz;
            qz = mul<T, TT>(t, t4);
        }

        vm = mod<T, TT>(qv, n);
        wm = mod<T, TT>(qw, n);
        xm = mod<T, TT>(qx, n);
        ym = mod<T, TT>(qy, n);
        zm = mod<T, TT>(qz, n);
    }
}

template <class T, class TT> static bool islnrc5prime(const T &n, int s = 0, int cid = 0)
{
    if (n < 23)
    {
        // prime for sure
        return (n == 1 || n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 || n == 19);
    }

    if (is_perfect_sursolid<T, TT>(n))
    {
        return false; // composite
    }

    T k, a, bv, bw, bx, by, bz;
    for (k = 1;; k++)
    {
        a = 19 + k * (k - 1);
        if (a % 5 != 1)
        {
            continue; // try another a
        }
        if (!isprime<T, TT>(a))
        {
            continue; // try another a
        }
        if (power<T, TT>(mod<T, TT>(n, a), (a - 1) / 5, a) == 1)
        {
            continue; // try another a
        }
        if (a == n)
        {
            return true; // n is a small prime
        }
        T g = gcd((2 * k - 1) * a * (2 * a - 1), n);
        if (g == n)
        {
            continue; // try another a
        }
        if (g > 1)
        {
            return false; //  composite
        }
        bv = 0;
        bw = 0;
        bx = 0;
        by = 1;
        bz = 0;
        exponentiate5<T, TT>(bv, bw, bx, by, bz, n - 1, n, a);
        if (bv == 0 && bw == 0 && bx == 0 && by == 0 && bz == 1)
        {
            // log this cornercase for verification purposes
            if (s | cid)
            {
                tlv_write(s, cid, TLV_B1, n);
            }
            continue; // try another a
        }
        break;
    }
    T bv2 = bv, bw2 = bw, bx2 = bx, by2 = by, bz2 = bz;
    exponentiate5<T, TT>(bv2, bw2, bx2, by2, bz2, 2, n, a); // B^2
    bv = add_mod<T, TT>(bv, bv2, n);
    bw = add_mod<T, TT>(bw, bw2, n);
    bx = add_mod<T, TT>(bx, bx2, n);
    by = add_mod<T, TT>(by, by2, n);
    bz = add_mod<T, TT>(bz, bz2, n);

    exponentiate5<T, TT>(bv2, bw2, bx2, by2, bz2, 2, n, a); // B^4
    bv = add_mod<T, TT>(bv, bv2, n);
    bw = add_mod<T, TT>(bw, bw2, n);
    bx = add_mod<T, TT>(bx, bx2, n);
    by = add_mod<T, TT>(by, by2, n);
    bz = add_mod<T, TT>(bz, bz2, n);

    bv2 = bv, bw2 = bw, bx2 = bx, by2 = by, bz2 = bz;
    exponentiate5<T, TT>(bv2, bw2, bx2, by2, bz2, 3, n, a); // B^3
    bv = add_mod<T, TT>(bv, bv2, n);
    bw = add_mod<T, TT>(bw, bw2, n);
    bx = add_mod<T, TT>(bx, bx2, n);
    by = add_mod<T, TT>(by, by2, n);
    bz = add_mod<T, TT>(bz, bz2, n);

    bz = add_mod<T, TT>(bz, 1, n);

    if (bv != n - 1 || bw != 1 || bx != n - 1 || by != 1 || bz != a)
        return false; // composite
    return true;      // might be prime
}

int inner_loop(int s, uint16_t cid, uint128_t seed, uint64_t count)
{
    bool r, rl;
    uint128_t v = 0;
    char buff[60];
    fflush(stdout);
    Lcg u;
    u.set_seed(seed);
    v = convert_seed_to_number(seed);
    sprintf(buff, "Client %4u Started ....", (unsigned)cid);
    print128(buff, v);
    fflush(stdout);
    while (count--)
    {
        v = u.get_seed(1);
        v = convert_seed_to_number(v);
        if (v >> 61)
        {
            assert(0);
        }
        else
        {
            r = isprime<uint64_t, uint128_t>(v);
            rl = islnrc2prime<uint64_t, uint128_t>(v, s, cid);
        }
        if (r)
        {
            if (!rl)
            {
                if (s)
                    tlv_write(s, cid, TLV_PSEUDOCOMPOSITE, v);
            }
        }
        else
        {
            if (rl)
            {
                if (s)
                    tlv_write(s, cid, TLV_PSEUDOPRIME, v);
            }
        }
    }

    sprintf(buff, "Client %4u Completed ..", (unsigned)cid);
    print128(buff, v);
    fflush(stdout);
    return 0;
}

static int inner_self_test_64(void)
{
    uint64_t s, r, t;
    bool b;

    printf("Perfect square ...\n");
    b = is_perfect_square<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(64);
    if (!b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(27);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002001);
    if (!b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002000);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002002);
    if (b)
        return -1;

    printf("Perfect cube ...\n");
    b = is_perfect_cube<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(64);
    if (!b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(81);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003001);
    if (!b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003000);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003002);
    if (b)
        return -1;

    printf("Perfect sursolid ...\n");
    b = is_perfect_sursolid<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(64 * 16);
    if (!b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(81);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005001ull);
    if (!b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005000ull);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005002ull);
    if (b)
        return -1;

    printf("Gcd ...\n");
    s = 12;
    t = 15;
    r = gcd<uint64_t>(s, t);
    if (r != 3)
        return -1;
    s = 12;
    t = 30;
    r = gcd<uint64_t>(s, t);
    if (r != 6)
        return -1;

    printf("Modinv ...\n");
    s = 11;
    t = 15;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 11)
        return -1;
    s = 12;
    t = 31;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 13)
        return -1;
    s = 1234567;
    t = 87654321;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 75327931)
        return -1;

    printf("Power ...\n");
    r = power<uint64_t, uint128_t>(3, 0xaa55, 197);
    if (r != 0xa7)
        return -1;

    printf("Known primes ...\n");
    b = isprime<uint64_t, uint128_t>(200003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(2000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(20000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(2000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(20000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(200000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }

    printf("Linear recurrence second order ...\n");
    if (1)
    {
        uint64_t s, t;
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 2, 101, 13);
        if (s != 4 || t != 17)
        {
            printf("linear recurrence 2 failed _2_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 4, 101, 13);
        if (s != 35 || t != 93)
        {
            printf("linear recurrence 2 failed _4_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 3, 101, 13);
        if (s != 25 || t != 86)
        {
            printf("linear recurrence 2 failed _3_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 5, 101, 13);
        if (s != 62 || t != 35)
        {
            printf("linear recurrence 2 failed _5_\n");
            return (-1);
        }
    }
    printf("Linear recurrence third order ...\n");
    if (1)
    {
        uint64_t s, t, u;
        s = 1;
        t = 2;
        u = 3;
        exponentiate3<uint64_t, uint128_t>(s, t, u, 2, 101, 13);
        if (s != 23 || t != 77 || u != 61)
        {
            printf("linear recurrence 3 failed _2_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        u = 3;
        exponentiate3<uint64_t, uint128_t>(s, t, u, 4, 101, 13);
        if (s != 58 || t != 0 || u != 75)
        {
            printf("linear recurrence 3 failed _4_\n");
            return (-1);
        }
        s = 0;
        t = 1;
        u = 2;
        exponentiate3<uint64_t, uint128_t>(s, t, u, 3, 101, 13);
        // printf("%lu %lu %lu\n", s, t, u);
        if (s != 6 || t != 25 || u != 21)
        {
            printf("linear recurrence 3 failed _3_\n");
            return (-1);
        }
        s = 0;
        t = 1;
        u = 2;
        exponentiate3<uint64_t, uint128_t>(s, t, u, 5, 101, 13);
        // printf("%lu %lu %lu\n", s, t, u);
        if (s != 21 || t != 91 || u != 14)
        {
            printf("linear recurrence 3 failed _5_\n");
            return (-1);
        }
    }
    printf("Linear recurrence fifth order ...\n");
    if (1)
    {
        uint64_t s, t, u, v, w;
        s = 1;
        t = 2;
        u = 3;
        v = 4;
        w = 5;
        // TODO
    }

    printf("Isprime ...\n");
    t = 1;
    t <<= 3;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }
    b = islnrc3prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc3 M(3) failed\n");
        return -1;
    }
    b = islnrc5prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc5 M(3) failed\n");
        return -1;
    }

    t = 101;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }
    b = islnrc3prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc3 101 failed\n");
        return -1;
    }
    b = islnrc5prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc5 101 failed\n");
        return -1;
    }

    t = 4493;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }
    b = islnrc3prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc3 4493 failed\n");
        return -1;
    }
    b = islnrc5prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc5 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }
    b = islnrc3prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc3 M(31) failed\n");
        return -1;
    }
    b = islnrc5prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc5 M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }
    b = islnrc3prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc3 M(61) failed\n");
        return -1;
    }

#define COUNT 5000
    volatile uint64_t t0, t1;
    t0 = __rdtsc();
    inner_loop(0, 2222, 0x4000000000, COUNT);
    t1 = __rdtsc();
    double d = (double)(t1 - t0);
    d /= COUNT;
    printf("Average: %8.1f ticks/iteration\n", d);
    printf("Self-test completed\n");
    // pass
    return 0;
}

int inner_self_test(void)
{
    int rc = 0;
    rc = inner_self_test_64();
    if (rc)
    {
        printf("64/128 bit failed\n");
        return rc;
    }
    return 0;
}
