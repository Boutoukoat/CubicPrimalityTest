
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#include "inner_loop.h"
#include "tlv.h"
#include "lcg.h"

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

template <class T, class TT> static T mul_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r *= v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T add_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r += v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T sub_mod(const T &u, const T &v, const T &q)
{
    if (u >= v)
    {
        return u - v;
    }
    T s = v;
    if (s >= q)
        s %= q;
    TT r = u;
    r += q - s;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T mod(const TT &u, const T &q)
{
    TT r = u;
    r %= q;
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

    /*
if ( 0x3c & (1 << (a % 7 ))) return false;
if ( 0xfc & (1 << (a % 9 ))) return false;
if ( 0xedc & (1 << (a % 13 ))) return false;
if ( 0x3e67c & (1 << (a % 19 ))) return false;
if ( 0x177e7ee8 & (1 << (a % 31 ))) return false;
if ( 0xf537fb2bcull & (1ull << (a % 37 ))) return false;
if ( 0xbcbfd99e66ff4f4ull & (1ull << (a % 61 ))) return false;
*/
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

template <class T, class TT> static T power(const T &a, const T &e, const T &m)
{
    T n = e;
    T power = a;
    T result = 1;

    while (n)
    {
        if (n & 1)
            result = mul_mod<T, TT>(result, power, m);
        power = square_mod<T, TT>(power, m);
        n >>= 1;
    }
    return result;
}

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

    return true; // might be prime
}

template <class T, class TT> static bool isprime(const T &n)
{
    if (!sieve<T>(n))
        return false;
    if (n <= 151 * 151)
        return true;

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
#if 0
	if (n < 341531)
		return witness < T, TT > (n, s, d, 9345883071009581737ull); 
	if (n < 1050535501ull)
		return witness < T, TT > (n, s, d, 336781006125ull) && witness < T, TT > (n, s, d, 9639812373923155ull);
	if (n < 350269456337ull)
		return witness < T, TT > (n, s, d, 4230279247111683200ull) && witness < T, TT > (n, s, d, 14694767155120705706ull)
			&& witness < T, TT > (n, s, d, 16641139526367750375ull); 
	if (n < 55245642489451ull)
		return witness < T, TT > (n, s, d, 2) && witness < T, TT > (n, s, d, 141889084524735ull)
			&& witness < T, TT > (n, s, d, 1199124725622454117ull) && witness < T, TT > (n, s, d, 11096072698276303650ull);
	if (n < 7999252175582851ull)
		return witness < T, TT > (n, s, d, 2) && witness < T, TT > (n, s, d, 4130806001517ull)
			&& witness < T, TT > (n, s, d, 149795463772692060ull) && witness < T, TT > (n, s, d, 186635894390467037ull)
			&& witness < T, TT > (n, s, d, 3967304179347715805ull); if (n < 585226005592931977ull)
		return witness < T, TT > (n, s, d, 2) && witness < T, TT > (n, s, d, 123635709730000ull)
			&& witness < T, TT > (n, s, d, 9233062284813009ull) && witness < T, TT > (n, s, d, 43835965440333360ull)
			&& witness < T, TT > (n, s, d, 761179012939631437ull) && witness < T, TT > (n, s, d, 1263739024124850375ull);
	return witness < T, TT > (n, s, d, 2) && witness < T, TT > (n, s, d, 325)
		&& witness < T, TT > (n, s, d, 9375) && witness < T, TT > (n, s, d, 28178)
		&& witness < T, TT > (n, s, d, 450885) && witness < T, TT > (n, s, d, 9780504)
		&& witness < T, TT > (n, s, d, 1795265022);
#endif
}

template <class T, class TT> void exponentiate(T &s, T &t, T &u, const T e, const T n, const T a)
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
        if (e & (1ull << bit))
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

template <class T, class TT> static bool islnrcprime(const T &n)
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
        exponentiate<T, TT>(bs, bt, bu, n - 1, n, a);
        if (bs == 1 && bt == 0 && bu == 0)
        {
            continue; // try another a
        }
        break;
    }
    T bs2 = bs, bt2 = bt, bu2 = bu;
    exponentiate<T, TT>(bs2, bt2, bu2, 2, n, a);
    bs = add_mod<T, TT>(bs, bs2, n);
    bt = add_mod<T, TT>(bt, bt2, n);
    bu = add_mod<T, TT>(bu, bu2, n);
    bu = add_mod<T, TT>(bu, 1, n);
    if (bs != n - 1 || bt != 1 || bu != a)
        return false; // composite
    return true; // might be prime
}

template <class T, class TT> static int loop(int s, uint16_t cid, const TT &seed, uint64_t count)
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

    while(count--)
    {
	    v = u.get_seed(1);
	    v = convert_seed_to_number(v);
        r = isprime<T, TT>(v);
        rl = islnrcprime<T, TT>(v);
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

int inner_loop(int s, uint16_t cid, uint128_t seed, uint64_t count)
{
        return loop<uint64_t, uint128_t>(s, cid, seed, count);
        // return loop < uint128_t, uint256_t > (s, cid, start, end);
}

static int inner_self_test_64(void)
{
    uint64_t s, r, t;
    bool b;

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

    printf("Power ...\n");
    r = power<uint64_t, uint128_t>(3, 0xaa55, 197);
    if (r != 0xa7)
        return -1;

    printf("Known primes ...\n");
    b = isprime<uint64_t, uint128_t>(200003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }
    b = isprime<uint64_t, uint128_t>(2000003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }
    b = isprime<uint64_t, uint128_t>(20000003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }
    b = isprime<uint64_t, uint128_t>(2000000000003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }
    b = isprime<uint64_t, uint128_t>(20000000000000003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }
    b = isprime<uint64_t, uint128_t>(200000000000000003ull);
    if (!b) { printf("expected prime failed\n"); return (-1); }

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
    b = islnrcprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc M(3) failed\n");
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
    b = islnrcprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc M(31) failed\n");
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
    b = islnrcprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc M(61) failed\n");
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
