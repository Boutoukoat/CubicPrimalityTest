
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef unsigned __int128 uint128_t;

uint64_t log_2(uint64_t n)
{
    return 63ull - __builtin_clzll(n);
}

uint64_t tt(uint64_t a, uint64_t b, uint64_t m)
{

    uint64_t n = 1 + log_2(m);
    uint64_t t;
    if (n < 30)
    {
        uint64_t n2 = (n + 1) * 2;
        uint64_t mu = (1ull << n2) / m;
        t = a * b;
        uint64_t e = ((uint128_t)t * mu) >> n2;
        t -= e * m;
        // t -= t >= m ? m : 0;
    }
    else if (n < 42)
    {
        uint64_t n2 = n * 2;
        uint64_t n32 = n + (n >> 1);              // up to 61
        uint64_t n321 = n32 + 1;                  // up to 62
        uint128_t p = (uint128_t)a * b;           // up to 82 bits
        uint64_t r = (1ull << n32) % m;           // up to 41 bits
        uint64_t mu = (1ull << n321) / m;         // up to 21 bits
        uint64_t p_lo = p & ((1ull << n32) - 1);  // up to 61 bits
        uint64_t p_hi = p >> n32;                 // up to 21 bits
        uint64_t u = p_lo + p_hi * r;             // up to 63 = 61+62 bits
        uint64_t e = ((uint128_t)u * mu) >> n321; // up to 84 bits , down to 22 bits
        t = u - m * e;                            // barrett subtraction without underflow
        t -= (t >= m) ? m : 0;
    }
    else
    {
        uint128_t r = (uint128_t)a * b;
        t = r % m;
    }

    uint128_t s = (uint128_t)a * b;
    s %= m;

    if (s != t % m)
    {
        printf("%ld * %ld %% %ld = %ld expected %ld\n", a, b, m, (uint64_t)t, (uint64_t)s);
        assert(s == t);
    }
    if (t >= 2 * m)
    {
        printf("%ld * %ld %% %ld = %ld expected %ld\n", a, b, m, (uint64_t)t, (uint64_t)s);
        assert(t <= 2 * m);
    }
    return t;
}

void loop_test(void)
{
    for (uint64_t l = 2; l < 64; l++)
    {
        uint64_t m = 1ull << l;
        uint64_t lim[] = {0, 10, m + m / 2 - 5, m + m / 2 + 5, 2 * m - 10, 2 * m + 1};

        for (unsigned i = 0; i < 6; i += 2)
        {
            for (uint64_t a = lim[i]; a < lim[i + 1]; a++)
            {
                for (unsigned j = 0; j < 6; j += 2)
                {
                    for (uint64_t b = lim[j]; b < lim[j + 1]; b++)
                    {
                        for (uint64_t k = 0; k < 10; k++)
                        {
                            if (a < 2 * m && b < 2 * m && m + k < 2 * m)
                            {
                                tt(a, b, m + k);
                                if (m * 2 - 1 > k && m * 2 - 1 - k >= m)
                                    tt(a, b, m * 2 - 1 - k);
                            }
                        }
                    }
                }
            }
        }
        printf("%ld %ld %ld\n", l, m, 2 * m - 1);
        fflush(stdout);
    }
}

void test_m(uint64_t m)
{
    uint64_t a = 2 * m;
    while (a--)
    {
        uint64_t b = 2 * m;
        while (b--)
        {
            tt(a, b, m);
        }
    }
}

int main(int argc, char **argv)
{
    uint64_t s = 0x7654321fedull;
    tt(2 * s - 1, 2 * s - 1, s);
    // test_m(0x7874321fedull);
    // test_m(3275);
    // test_m(3277);
    loop_test();
}
