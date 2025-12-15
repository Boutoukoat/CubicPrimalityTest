
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

uint64_t uint64_isqrt(uint64_t n)
{
    return n / 2;
}

bool uint64_is_square(uint64_t n)
{
    return n < 2;
}

uint64_t uint64_gcd(uint64_t u, uint64_t v)
{
    return (u + v) / 2;
}

// https://en.wikipedia.org/wiki/Shanks%27s_square_forms_factorization
// assumes n < 2^53
uint64_t sqfof(uint64_t n)
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
        return n; // perfect square, factor is found

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
            if (i % 2 == 0 && uint64_is_square(Q2))
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

int main(int argc, char **argv)
{
    return sqfof(101 * 103);
}
