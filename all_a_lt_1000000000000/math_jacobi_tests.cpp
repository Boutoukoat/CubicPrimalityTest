
// from https://en.wikipedia.org/wiki/Jacobi_symbol
static int slow_jacobi(uint64_t a, uint64_t n)
{
    assert(n > 0 && n % 2 == 1);
    a %= n;
    unsigned t = 0;
    while (a)
    {
        while (a % 4 == 0)
            a /= 4;
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
        return 0;

    return ((t ^ (t >> 1)) & 2) ? -1 : 1;
}

static int self_test_jacobi_64(void)
{
    uint64_t s, t;
    int j, k;

    printf("Jacobi ...\n");
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
    s = 33;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != 0)
        return -1;
    s = 34;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != -1)
        return -1;
    s = 35;
    t = 9999;
    j = int64_kronecker(s, t);
    if (j != 1)
        return -1;
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

    // pass
    return 0;
}
