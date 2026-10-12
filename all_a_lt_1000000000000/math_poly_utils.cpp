
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <x86intrin.h>

using namespace std;

typedef unsigned __int128 uint128_t;
typedef signed __int128 int128_t;

typedef struct
{
    uint64_t u; /* coefficient of x */
    uint64_t v; /* constant coefficient */
} uint64_poly;

/* All arguments are reduced modulo m */
static inline uint64_t addm(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t t = a + b;
    if (t < a)
    {
        return t - m;
    }
    return (t >= m) ? t - m : t;
}

static inline uint64_t mulm(uint64_t a, uint64_t b, uint64_t m)
{
    return mul_mod(a, b, m);
}

/* Square (u*x + v) modulo x^2 + a*x + c */
static inline uint64_poly uint64_polysqr(uint64_poly p, uint64_t a, uint64_t c, uint64_t m)
{
    uint64_t u2 = m - mulm(p.u, p.u, m);
    uint64_t v2 = mulm(p.v, p.v, m);
    uint64_t uv = mulm(p.u, p.v, m);

    uint64_poly r;
    r.u = addm(addm(uv, uv, m), mulm(a, u2, m), m);
    r.v = addm(v2, mulm(c, u2, m), m);
    return r;
}

/* Multiply (u*x + v) by x modulo x^2 + a*x + c */
static inline uint64_poly uint64_polymul_x(uint64_poly p, uint64_t a, uint64_t c, uint64_t m)
{
    uint64_poly r;
    uint64_t t = m - p.u;
    r.u = addm(p.v, mulm(a, t, m), m);
    r.v = mulm(c, t, m);
    return r;
}

/* Multiply (u*x + v) by (s*x + t), modulo x^2 + a*x + c */
static uint64_poly uint64_polymul(uint64_poly p, uint64_poly q, uint64_t a, uint64_t c, uint64_t m)
{
    uint64_t us = m - mulm(p.u, q.u, m);
    uint64_t vt = mulm(p.v, q.v, m);
    uint64_t ut = mulm(p.u, q.v, m);
    uint64_t vs = mulm(p.v, q.u, m);

    uint64_poly r;
    r.u = addm(addm(ut, vs, m), mulm(a, us, m), m);
    r.v = addm(vt, mulm(c, us, m), m);
    return r;
}

/* Left-to-right binary exponentiation specialised to base x+0, modulo x^2 + a*x + c */
static uint64_poly uint64_polypow(uint64_poly b, uint64_t e, uint64_t a, uint64_t c, uint64_t m)
{
    uint64_poly r = b;
    if (e < 3)
    {
        if (e == 0)
        {
            return (uint64_poly){0, 1};
        }
        if (e == 1)
        {
            return r;
        }
        // e == 2
        return uint64_polysqr(r, a, c, m);
    }

    uint64_t bit = 1ull << (uint64_log_2(e) - 1);

    if (r.u == 1 && r.v == 0)
    {
        // exponentiaion base is x+0
        for (; bit != 0; bit >>= 1)
        {
            // square
            r = uint64_polysqr(r, a, c, m);
            if (e & bit)
            {
                // multiply by x
                r = uint64_polymul_x(r, a, c, m);
            }
        }
    }
    else
    {
        // exponentiation base is any
        for (; bit != 0; bit >>= 1)
        {
            // square
            r = uint64_polysqr(r, a, c, m);
            if (e & bit)
            {
                // slow multiply
                r = uint64_polymul(r, b, a, c, m);
            }
        }
    }

    return r;
}
