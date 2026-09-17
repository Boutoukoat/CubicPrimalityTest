
#include "math_utils.cpp"

#include "iostream"

/*

a[0] = 3
a[n] = 2^a[n - 1] - 3

*/

static int test_k1(bool verbose)
{
    printf("Misc tests\n");

    uint64_t f = 1;
    uint64_t k = 3;
    printf("a[0] = %ld\n", k);
    printf("a[n] = 2^a[n - 1] - %ld\n", k);
    printf("\n");
    while (1)
    {
        f += 2;
        if (uint64_small_factor(f) != 1 || !uint64_is_prime_mr(f))
        {
            continue;
        }

        uint64_t b = k;
        uint64_t a = (1l << b) - k;
        uint64_t cpt = 0;
        unsigned c0 = 0, c1 = 0;
        unsigned i = 2;
        while (i < 25000)
        {
            // wrong code
            b = uint64_pow2_mod(b, f - 1);
            b = (b < k) ? b + f - 1 - k : b - 3;

            a = uint64_pow2_mod(a, f);
            a = (a < k) ? a + f - k : a - 3;

            if (a == 0)
            {
                cpt += 1;
                if (cpt == 1)
                {
                    c0 = i;
                }
                else
                {
                    c1 = i;
                    printf("%lu divides a[i*%u+%u]\n", f, c1 - c0, c0);
                    break;
                }
            }

            if (cpt == 0 && i > 1500)
            {
                break;
            }
            i++;
        }
        if (cpt == 1)
        {
            printf("%lu divides a[%u]\n", f, c0);
        }
    }

    return 0;
}

#define LIM (200000000)
#define MAX_FACTOR 5000000
static bool nnn[LIM];

static void eliminate(unsigned c1, unsigned c0)
{
    for (unsigned j = c0; j < LIM; j += c1)
    {
        nnn[j] = false;
    }
}

static void uint64_all_divisors(const factor_v &v, vector<uint64_t> &d)
{
    unsigned i, j, k;
    vector<uint64_t> u;
    i = v.size();
    while (i--)
    {
        uint64_t f = v[i].prime;
        j = v[i].count;
        while (j--)
        {
            u.push_back(f);
        }
    }
    j = (1 << u.size());
    i = 1;
    while (i < j)
    {
        k = u.size();
        uint64_t f = 1;
        while (k--)
        {
            if ((1 << k) & i)
            {
                f *= u[k];
            }
        }
        d.push_back(f);
        i++;
    }
    sort(d.begin(), d.end());
    d.erase(unique(d.begin(), d.end()), d.end());
}

static int test_64(uint64_t k, int plus, bool verbose)
{
    uint64_t f = 1;
    factor_v v;
    vector<uint64_t> d;
    memset(nnn, true, LIM);
    char pplus[6];
    snprintf(pplus, 5, "%s%d", (plus >= 0 ? "+" : ""), plus);

    if (verbose)
    {
        printf("Sieve small factors of %lu*2^n%s\n", k, pplus);
    }

    while (f < MAX_FACTOR)
    {
        f += 2;
        if (uint64_small_factor(f) != 1 || !uint64_is_prime_mr(f))
        {
            continue;
        }
        unsigned j;
        uint64_t q, i;
        if ((f & 7) == 1 || (f & 7) == 7)
        {
            if (int64_kronecker((plus > 0 ? -(int64_t)k : k), f) != 1)
            {
                if (verbose)
                {
                    printf("%lu never divides %lu*2^n%s (jacobi)\n", f, k, pplus);
                }
                continue;
            }
        }
        uint64_t t = uint64_mod_inv(k, f);
        if (t == 0)
        {
            i = 0;
            q = f + 1;
        }
        else
        {
            if (plus > 0)
            {
                t = f - t;
            }
            i = 0;
            q = 1;
            while (q != t && i < f)
            {
                q *= 2;
                q -= (q >= f) ? f : 0;
                i += 1;
                if (q == 1)
                {
                    break;
                }
            }
        }
        if (q != t)
        {
            if (verbose)
            {
                printf("%lu never divides %lu*2^n%s (discrete log)\n", f, k, pplus);
            }
            continue;
        }

        v.clear();
        d.clear();
        uint64_all_factors(v, f - 1);
        uint64_all_divisors(v, d);
        uint64_t k2, k1 = 0;
        uint64_t a = (plus ? f - 2 : 0);
        uint64_t dj = d.size() - 1;
        for (j = 0; j < d.size() - 1; j++)
        {
            k2 = d[j];
            a = mul_add_mod(uint64_pow2_mod(k2 - k1, f), a + 1, f - 1, f);
            k1 = k2;
            if (a + 2 == f && plus > 0)
            {
                dj = k2;
                break;
            }
            if (a == 0 && plus < 0)
            {
                dj = k2;
                break;
            }
        }
        if (verbose)
        {
            printf("%lu divides %lu*2^(i*%lu+%lu)%s\n", f, k, dj, i, pplus);
        }
        eliminate(d[j], i);
    }

    printf("====\n");
    for (unsigned j = 0; j < LIM; j++)
    {
        if (nnn[j])
        {
            printf("%lu*2^%u %s\n", k, j, pplus);
        }
    }

    return 0;
}

int main(int argc, char **argv)
{
    // const uint64_t k = 2996863034895ul;
    int64_t k = 5502699578414355ul;
    // const uint64_t k = 23669;
    int r = 1;
    bool err = false;
    bool verbose = false;

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-v"))
        {
            verbose = true;
            continue;
        }
        if (!strcmp(argv[i], "-k"))
        {
            k = strtol(argv[++i], 0, 10);
            continue;
        }
        if (!strcmp(argv[i], "-r"))
        {
            r = (int)strtol(argv[++i], 0, 10);
            continue;
        }
        err = true;
    }

    if (k == 0 || r == 0)
    {
        printf("k must be > 0 and |r| must be > 0\n");
        err = true;
    }
    if (k < 0)
    {
        printf("k must be positive\n");
        err = true;
    }
    if ((k & 1) == 0)
    {
        printf("k must be odd\n");
        err = true;
    }
    if (r < -1 || r > 1)
    {
        printf("|r| must be 1\n");
    }

    if (err)
    {
        printf("invalid command line\n");
        exit(1);
    }

    test_64(k, r, verbose);
    return 0;
}
