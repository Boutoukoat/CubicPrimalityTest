
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
    uint64_t n2 = n * 2;
    uint64_t mu = (1ull << n2) / m;
    uint64_t r = a * b;
    uint64_t e = ((uint128_t)r * mu) >> n2;
    r -= e * m;
    // r -= r >= m ? m : 0;

    uint128_t s = (uint128_t)a * b;
    s %= m;

    if (s != r % m)
    {
	    printf("%ld * %ld %% %ld = %ld expected %ld\n", a, b, m, (uint64_t)r, (uint64_t)s);
        assert(s == r);
    }
    if ( r >= 2 * m)
    {
	    printf("%ld * %ld %% %ld = %ld expected %ld\n", a, b, m, (uint64_t)r, (uint64_t)s);
        assert(r <= 2 * m);
    }
    return r;
}

int main(int argc, char **argv)
{
    for (uint64_t l = 2; l < 64; l++)
    {
        uint64_t m = 1ull << l;
        uint64_t lim[] = {0, 10, 2 * m - 10, 2 * m + 1};

        for (unsigned i = 0; i < 4; i += 2)
        {
            for (uint64_t a = lim[i]; a < lim[i + 1]; a++)
            {
                for (unsigned j = 0; j < 4; j += 2)
                {
                    for (uint64_t b = lim[j]; b < lim[j + 1]; b++)
                    {
			    for (uint64_t k = 0; k < 10; k++)
			    {
				    if (a < 2 * m && b < 2 * m && m + k < 2 * m)
				    {
                        tt(a, b, m + k);
                        if (m * 2 - 1 > k && m * 2 -1 -k >= m ) tt(a, b, m * 2 - 1 - k);
				    }
			    }
                    }
                }
            }
        }
        printf("%ld %ld %ld\n", l, m, 2*m-1);
        fflush(stdout);
    }
}
