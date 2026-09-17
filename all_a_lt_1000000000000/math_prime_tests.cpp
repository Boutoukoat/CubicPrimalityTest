

static int self_test_prime_64(void)
{
    bool b;
    uint64_t t;

    printf("Known primes ...\n");
    b = uint64_is_prime_mr(200003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime_mr(2000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime_mr(20000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime_mr(2000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime_mr(20000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime_mr(200000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }

    printf("Isprime (MR) ...\n");

    t = 17 * 19;
    b = uint64_is_prime_bpsw(t);
    if (b)
    {
        printf("isprime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }

    printf("Isprime (NIST) ...\n");

    t = 17 * 19;
    b = uint64_is_prime_nist(t);
    if (b)
    {
        printf("isprime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }

    printf("Isprime (BPSW) ...\n");

    t = 17 * 19;
    b = uint64_is_prime_bpsw(t);
    if (b)
    {
        printf("isprime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }

    for (t = 3; t < 1001; t += 2)
    {
        bool b1 = uint64_is_prime_mr(t);
        bool b2	= uint64_is_prime_nist(t);
        bool b3	= uint64_is_prime_bpsw(t);
        if (b1 != b2 || b1 != b3)
        {
            printf("MR/BPSW/NIST difference for prime %lu\n", t);
            return -1;
        }
    }

    for (t = 10000000000003; t < 10000000001001; t += 2)
    {
        bool b1 = uint64_is_prime_mr(t);
        bool b2	= uint64_is_prime_nist(t);
        bool b3	= uint64_is_prime_bpsw(t);
        if (b1 != b2 || b1 != b3)
        {
            printf("MR/BPSW/NIST difference for prime %lu\n", t);
            return -1;
        }
    }

    return 0;
}
