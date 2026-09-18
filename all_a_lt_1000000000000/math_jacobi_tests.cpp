
// from https://en.wikipedia.org/wiki/Jacobi_symbol
static int slow_jacobi(uint64_t a, uint64_t n)
{
    assert(n > 0 && n % 2 == 1);
    a %= n;
    unsigned t = 0;
    while (a)
    {
        while (a % 4 == 0)
        {
            a /= 4;
        }
        if (a % 2 == 0)
        {
            t ^= n;
            a /= 2;
        }
        t ^= (a & n) & 2;
        n %= a;
        n ^= a;
        a ^= n;
        n ^= a;
    }

    if (n != 1)
    {
        return 0;
    }

    return ((t ^ (t >> 1)) & 2) ? -1 : 1;
}

static int minus_one_kronecker(uint64_t n)
{
    assert(n > 0 && n % 2 == 1);
    return (n & 2) ? -1 : 1;
}

static int slow_kronecker128(int128_t a, int128_t b)
{
    int k = 1;

    /* (a|0) */
    if (b == 0)
        return (a == (int128_t)1 || a == (int128_t)-1) ? 1 : 0;

    // make b positive
    // K(a,-b) = K(a,b)*sgn(a)
    uint128_t ub;
    if (b < 0)
    {
        ub = -b;
        k = (a < 0) ? -k : k;
    }
    else
    {
        ub = b;
    }

    /* Handle negative numerator. */
    uint128_t ua;
    if (a < 0)
    {
        ua = -a;
        k = ((ub & 3) == 3) ? -k : k;
    }
    else
    {
        ua = a;
    }

    /* Remove factors of two from denominator. */
    if ((ub & 1) == 0)
    {
        int s = uint128_tzcnt(ub);
        ub >>= s;

        if ((ua & 1) == 0)
            return 0; /* gcd(a,b)>1 */

        uint64_t amod8 = (uint64_t)(ua & 7);
        if (s & 1)
        {
            k = (amod8 == 3 || amod8 == 5) ? -k : k;
        }
    }

    ua %= ub;

    while (ua)
    {

        /* Remove factors of two from numerator. */
        int s = uint128_tzcnt(ua);
        ua >>= s;

        if (s & 1)
        {
            uint64_t bmod8 = (uint64_t)(ub & 7);
            k = (bmod8 == 3 || bmod8 == 5) ? -k : k;
        }

        /* Quadratic reciprocity. */
        k = ((ua & 3) == 3 && (ub & 3) == 3) ? -k : k;

        uint128_t t = ua;
        ua = ub % ua;
        ub = t;
    }

    return (ub == 1) ? k : 0;
}

static int self_test_jacobi_64(void)
{
    uint64_t s, t;
    int j, k;

    printf("Jacobi ...\n");
    s = 5;
    t = 11;
    j = uint64_jacobi(s, t);
    if (j != 1)
    {
        printf("failed Jacobi(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 33;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != 0)
    {
        printf("failed Jacobi(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 34;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != -1)
    {
        printf("failed Jacobi(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 35;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != 1)
    {
        printf("failed Jacobi(%lu, %lu)\n", s, t);
        return -1;
    }

    // all pairs (s,t) < 101
    for (t = 1; t < 101; t += 2)
    {
        for (s = 1; s <= t; s += 2)
        {
            j = uint64_jacobi(s, t);
            k = slow_jacobi(s, t);
            if (j != k)
            {
                printf("failed Jacobi(%lu, %lu)\n", s, t);
                return -1;
            }
        }
    }

    printf("Kronecker ...\n");
    s = 5;
    t = 11;
    j = int64_kronecker(s, t);
    if (j != 1)
    {
        printf("failed kronecker(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 33;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != 0)
    {
        printf("failed kronecker(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 34;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != -1)
    {
        printf("failed kronecker(%lu, %lu)\n", s, t);
        return -1;
    }
    s = 35;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != 1)
    {
        printf("failed kronecker(%lu, %lu)\n", s, t);
        return -1;
    }

    j = int64_kronecker(11, 101);
    if (j != -1)
        return -1;
    j = int64_kronecker(-11, 101);
    if (j != -1)
        return -1;
    j = int64_kronecker(13, 101);
    if (j != 1)
        return -1;
    j = int64_kronecker(-13, 101);
    if (j != 1)
        return -1;
    j = int64_kronecker(-1, 101);
    if (j != 1)
        return -1;
    j = int64_kronecker(0, 101);
    if (j != 0)
        return -1;
    j = int64_kronecker(1, 101);
    if (j != 1)
        return -1;
    j = int64_kronecker(1, 0);
    if (j != 1)
        return -1;
    j = int64_kronecker(2, 0);
    if (j != 0)
        return -1;
    j = int64_kronecker(13, -101);
    if (j != 1)
        return -1;
    j = int64_kronecker(-13, -101);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, -11);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, -9);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, -7);
    if (j != 1)
        return -1;
    j = int64_kronecker(-2, -5);
    if (j != 1)
        return -1;
    j = int64_kronecker(-2, -3);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, -1);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, 1);
    if (j != 1)
        return 1;
    j = int64_kronecker(-2, 3);
    if (j != 1)
        return 1;
    j = int64_kronecker(-2, 5);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, 7);
    if (j != -1)
        return -1;
    j = int64_kronecker(-2, 9);
    if (j != 1)
        return -1;
    j = int64_kronecker(-2, 11);
    if (j != 1)
        return -1;
    j = int64_kronecker(2, 9);
    if (j != 1)
        return -1;
    j = int64_kronecker(2, -9);
    if (j != 1)
        return -1;
    j = int64_kronecker(2, 11);
    if (j != -1)
        return -1;
    j = int64_kronecker(2, -11);
    if (j != -1)
        return -1;
    j = int64_kronecker(3, 11);
    if (j != 1)
        return -1;
    j = int64_kronecker(-3, 11);
    if (j != -1)
        return -1;
    j = int64_kronecker(3, 13);
    if (j != 1)
        return -1;
    j = int64_kronecker(-3, 13);
    if (j != 1)
        return -1;
    j = int64_kronecker(3, 15);
    if (j != 0)
        return -1;
    j = int64_kronecker(-3, 15);
    if (j != 0)
        return -1;

    // check more negative numbers
    for (int64_t tt = 201; tt < 301; tt += 2)
    {
        for (int64_t ss = -20; ss < 20; ss++)
        {
            if (ss)
            {
                int j = int64_kronecker(ss, tt);
                int i = ss < 0 ? uint64_jacobi(-ss, tt) * minus_one_kronecker(tt) : uint64_jacobi(ss, tt);
                if (i != j)
                {
                    printf("Inconsistent jacobi/kronecker result (64 bits)\n");
                    return -1;
                }
            }
        }
    }

    // check more negative numbers
    for (int128_t tt = -20; tt < 20; tt++)
    {
        if (tt)
        {
            for (int128_t ss = -20; ss < 20; ss++)
            {
                if (ss)
                {
                    int j = int128_kronecker(ss, tt);
                    int i = int64_kronecker(ss, tt);
                    if (i != j)
                    {
                        printf("Inconsistent jacobi/kronecker result (128 bits small numbers)\n");
                        return -1;
                    }
                }
            }
        }
    }

    // check more large and negative numbers
    for (unsigned u = 50; u < 70; u++)
    {
        int128_t uu = 1;
        uu <<= u;

        for (int128_t tt = uu - 20; tt < uu + 20; tt++)
        {
            if (tt)
            {
                for (int128_t ss = uu - 20; ss < uu + 20; ss++)
                {
                    if (ss)
                    {
                        int j = int128_kronecker(ss, tt);
                        int i = slow_kronecker128(ss, tt);
                        if (i != j)
                        {
                            printf("Inconsistent jacobi/kronecker result (128 bits mid numbers)\n");
                            printf("ss 0x%16.16lx%16.16lx\n", (uint64_t)(ss >> 64), (uint64_t)ss);
                            printf("tt 0x%16.16lx%16.16lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
                            printf("slow %d fast %d\n", i, j);
                            return -1;
                        }
                    }
                }
            }
        }
    }
    // pass
    return 0;
}
