

static int self_test_prime_64(void)
{
    bool b;
    uint64_t t, s;

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
        printf("is_prime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_mr(t);
    if (!b)
    {
        printf("is_prime M(61) failed\n");
        return -1;
    }

    printf("Isprime (NIST) ...\n");

    t = 17 * 19;
    b = uint64_is_prime_nist(t);
    if (b)
    {
        printf("is_prime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_nist(t);
    if (!b)
    {
        printf("is_prime M(61) failed\n");
        return -1;
    }

    printf("Isprime (BPSW) ...\n");

    t = 17 * 19;
    b = uint64_is_prime_bpsw(t);
    if (b)
    {
        printf("is_prime 323 failed\n");
        return -1;
    }

    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime M(3) failed\n");
        return -1;
    }

    t = 11;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime 11 failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime_bpsw(t);
    if (!b)
    {
        printf("is_prime M(61) failed\n");
        return -1;
    }

    for (t = 3; t < 1001; t += 2)
    {
        bool b1 = uint64_is_prime_mr(t);
        bool b2 = uint64_is_prime_nist(t);
        bool b3 = uint64_is_prime_bpsw(t);
        if (b1 != b2 || b1 != b3)
        {
            printf("MR/BPSW/NIST difference for prime %lu\n", t);
            return -1;
        }
    }

    for (t = 10000000000003; t < 10000000001001; t += 2)
    {
        bool b1 = uint64_is_prime_mr(t);
        bool b2 = uint64_is_prime_nist(t);
        bool b3 = uint64_is_prime_bpsw(t);
        if (b1 != b2 || b1 != b3)
        {
            printf("MR/BPSW/NIST difference for prime %lu\n", t);
            return -1;
        }
    }

    // checks at modulus 2^63
    s = (uint64_t)(-1ull) >> 1;
    for (t = s | 1; t > s - 4001; t -= 2)
    {
        bool b1 = uint64_is_prime_mr(t);
        bool b2 = uint64_is_prime_nist(t);
        bool b3 = uint64_is_prime_bpsw(t);
        if (b1 != b2 || b1 != b3)
        {
            printf("MR/BPSW/NIST difference for 63 bits 0x%lx\n", t);
            return -1;
        }
    }

    return 0;
}
static int self_test_prime_128(void)
{
    uint128_t p_lo, p_hi, p;
    bool b;

    printf("Isprime 128 (NIST) ...\n");

    // prime
    p_lo = 0x0000000000000005ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000005ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000007ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000007ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000000000bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000000000bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000017ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000017ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000025ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000025ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000049ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000049ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000000004dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000000004dul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000089ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000089ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000000077ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000000077ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000115ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000115ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000000000ddul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000000000ddul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000000021dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000000021dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000000001e1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000000001e1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000425ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000425ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000000047bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000000047bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000000821ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000000821ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000000081dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000000081dul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000001051ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000001051ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000000ff7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000000ff7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000002047ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000002047ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000001f37ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000001f37ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000004087ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000004087ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000000040fdul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000000040fdul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000000808dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000000808dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000007f7ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000007f7ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000010111ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000010111ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000000fbfbul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000000fbfbul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000002011dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000002011dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000001fed3ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000001fed3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000040201ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000040201ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000040be5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000040be5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000080201ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000080201ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000000801ebul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000000801ebul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000100403ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000100403ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000100febul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000100febul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000200407ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000200407ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000001ffbf1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000001ffbf1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000040080bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000040080bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000003fdfd3ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000003fdfd3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000000800807ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000000800807ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000007f87e5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000007f87e5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000100100ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000100100ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000000fffff7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000000fffff7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000000200103dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000000200103dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000200afcdul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000200afcdul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000004002011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000004002011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000000401ffeful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000000401ffeful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000008002009ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000008002009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000008031fe5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000008031fe5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000010004015ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000010004015ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000001005ffaful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000001005ffaful  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000002000401ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000002000401ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000001fff3ff7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000001fff3ff7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000040008031ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000040008031ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000003ff7ffc7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000003ff7ffc7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000080008003ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000080008003ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000007fed7fedul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000007fed7fedul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000100010005ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000100010005ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000000fff1fff1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000000fff1fff1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000200010049ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000200010049ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000001fffefe4dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000001fffefe4dul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000400020011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000400020011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000040037ffe3ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000040037ffe3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000800020003ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000000800020003ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000080001fffdul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000080001fffdul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000001000040009ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000001000040009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000ffff7fff1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000000ffff7fff1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000200004003dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000000200004003dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000020002bff97ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000020002bff97ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000004000080005ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000004000080005ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000040009fffebul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000040009fffebul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000008000080003ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000008000080003ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000800027fff9ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000000800027fff9ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000010000100045ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000010000100045ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000100003fffebul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00000100003fffebul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000002000010001ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000002000010001ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000020000afffcdul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000020000afffcdul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000040000200011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000040000200011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000040000ffff67ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000040000ffff67ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000080000200009ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000080000200009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000007ffff9fff79ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000007ffff9fff79ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000100000400003ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000100000400003ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000100002ffffd3ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000100002ffffd3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000020000040000dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000020000040000dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000200000bfffe5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0000200000bfffe5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000400000800013ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0000400000800013ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00003ffffcffff79ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00003ffffcffff79ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000080000080001ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x000080000080001ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00008000067ffd7bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00008000067ffd7bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0001000001000011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0001000001000011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0001000027ffff7ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0001000027ffff7ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0002000001000059ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0002000001000059ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000200001cffff97ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x000200001cffff97ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0004000002000005ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0004000002000005ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0003fffff7fffaabul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0003fffff7fffaabul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0008000002000001ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0008000002000001ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0007ffff81fffdb7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0007ffff81fffdb7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x00100000040000b1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x00100000040000b1ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0010000027ffffb5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0010000027ffffb5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0020000004000011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0020000004000011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x002000004bffff6ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x002000004bffff6ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0040000008000011ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0040000008000011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x003fffffaffffb95ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x003fffffaffffb95ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0080000008000015ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x0080000008000015ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x007ffffda7ffff8bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x007ffffda7ffff8bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x010000001000000bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x010000001000000bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00fffffc9fffff55ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x00fffffc9fffff55ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x020000001000002bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x020000001000002bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x01fffff98ffffd8dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x01fffff98ffffd8dul  \n");
        return -1;
    }

    // prime
    p_lo = 0x040000002000005bul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x040000002000005bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x04000000ffffffdful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x04000000ffffffdful  \n");
        return -1;
    }

    // prime
    p_lo = 0x080000002000001ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x080000002000001ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x07ffffff9ffffff7ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x07ffffff9ffffff7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x1000000040000047ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x1000000040000047ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0ffffff7ffffff97ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x0ffffff7ffffff97ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x200000004000001ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x200000004000001ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x1ffffff13ffffe7ful;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x1ffffff13ffffe7ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x4000000080000023ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x4000000080000023ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x40000004fffffff5ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x40000004fffffff5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x800000008000000dul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000000ul 0x800000008000000dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x800000067ffffff1ul;
    p_hi = 0x0000000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000000ul 0x800000067ffffff1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000100000095ul;
    p_hi = 0x0000000000000001ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000001ul 0x0000000100000095ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000009ffffffb5ul;
    p_hi = 0x0000000000000001ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000001ul 0x00000009ffffffb5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000100000005ul;
    p_hi = 0x0000000000000002ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000002ul 0x0000000100000005ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000006ffffffabul;
    p_hi = 0x0000000000000002ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000002ul 0x00000006ffffffabul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000200000019ul;
    p_hi = 0x0000000000000004ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000004ul 0x0000000200000019ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000fffffff67ul;
    p_hi = 0x0000000000000004ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000004ul 0x0000000fffffff67ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000200000057ul;
    p_hi = 0x0000000000000008ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000008ul 0x0000000200000057ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000000dffffff1ful;
    p_hi = 0x0000000000000008ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000008ul 0x0000000dffffff1ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000040000005dul;
    p_hi = 0x0000000000000010ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000010ul 0x000000040000005dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffffbffffffbfful;
    p_hi = 0x000000000000000ful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000000000ful 0xffffffbffffffbfful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000400000001ul;
    p_hi = 0x0000000000000020ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000020ul 0x0000000400000001ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffff8bfffff783ul;
    p_hi = 0x000000000000001ful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000000001ful 0xffffff8bfffff783ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000800000011ul;
    p_hi = 0x0000000000000040ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000040ul 0x0000000800000011ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000000affffff995ul;
    p_hi = 0x0000000000000040ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000040ul 0x000000affffff995ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000000800000007ul;
    p_hi = 0x0000000000000080ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000080ul 0x0000000800000007ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffff07fffffc3ful;
    p_hi = 0x000000000000007ful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000000007ful 0xffffff07fffffc3ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000100000000ful;
    p_hi = 0x0000000000000100ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000100ul 0x000000100000000ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000019fffffff65ul;
    p_hi = 0x0000000000000100ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000100ul 0x0000019fffffff65ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000001000000001ul;
    p_hi = 0x0000000000000200ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000200ul 0x0000001000000001ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffffefffffffd3ul;
    p_hi = 0x00000000000001fful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000000001fful 0xffffffefffffffd3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000002000000017ul;
    p_hi = 0x0000000000000400ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000400ul 0x0000002000000017ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfffffdffffffff1ful;
    p_hi = 0x00000000000003fful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000000003fful 0xfffffdffffffff1ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000002000000013ul;
    p_hi = 0x0000000000000800ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000000800ul 0x0000002000000013ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfffffa9fffffff51ul;
    p_hi = 0x00000000000007fful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000000007fful 0xfffffa9fffffff51ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000400000002ful;
    p_hi = 0x0000000000001000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000001000ul 0x000000400000002ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfffff67ffffffec5ul;
    p_hi = 0x0000000000000ffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000000ffful 0xfffff67ffffffec5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000004000000023ul;
    p_hi = 0x0000000000002000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000002000ul 0x0000004000000023ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffef3ffffffbf5ul;
    p_hi = 0x0000000000001ffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000001ffful 0xffffef3ffffffbf5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000008000000055ul;
    p_hi = 0x0000000000004000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000004000ul 0x0000008000000055ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x000007ffffffff5ful;
    p_hi = 0x0000000000004000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000004000ul 0x000007ffffffff5ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x000000800000000ful;
    p_hi = 0x0000000000008000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000008000ul 0x000000800000000ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x0000007fffffff97ul;
    p_hi = 0x0000000000008000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000008000ul 0x0000007fffffff97ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000010000000065ul;
    p_hi = 0x0000000000010000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000010000ul 0x0000010000000065ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffb7fffffffae7ul;
    p_hi = 0x000000000000fffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000000fffful 0xffffb7fffffffae7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000001000000001dul;
    p_hi = 0x0000000000020000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000020000ul 0x000001000000001dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffff6cfffffff6d3ul;
    p_hi = 0x000000000001fffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000001fffful 0xffff6cfffffff6d3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000002000000000dul;
    p_hi = 0x0000000000040000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000040000ul 0x000002000000000dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000bfffffffdc9ul;
    p_hi = 0x0000000000040000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000040000ul 0x00000bfffffffdc9ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x00000200000000bdul;
    p_hi = 0x0000000000080000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000080000ul 0x00000200000000bdul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffffc9fffffffec5ul;
    p_hi = 0x000000000007fffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000007fffful 0xffffc9fffffffec5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x00000400000000dbul;
    p_hi = 0x0000000000100000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000100000ul 0x00000400000000dbul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00000fffffffff5bul;
    p_hi = 0x0000000000100000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000100000ul 0x00000fffffffff5bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000040000000001ul;
    p_hi = 0x0000000000200000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000200000ul 0x0000040000000001ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00001bfffffffec1ul;
    p_hi = 0x0000000000200000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000200000ul 0x00001bfffffffec1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000008000000005ful;
    p_hi = 0x0000000000400000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000400000ul 0x000008000000005ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffff1ffffffff98bul;
    p_hi = 0x00000000003ffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000003ffffful 0xffff1ffffffff98bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000008000000008dul;
    p_hi = 0x0000000000800000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000000800000ul 0x000008000000008dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfffca7fffffffe71ul;
    p_hi = 0x00000000007ffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000007ffffful 0xfffca7fffffffe71ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000100000000009ul;
    p_hi = 0x0000000001000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000001000000ul 0x0000100000000009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffff5fffffffff89ul;
    p_hi = 0x0000000000fffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000000fffffful 0xffff5fffffffff89ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000100000000017ul;
    p_hi = 0x0000000002000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000002000000ul 0x0000100000000017ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00018ffffffffc15ul;
    p_hi = 0x0000000002000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000002000000ul 0x00018ffffffffc15ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000020000000000dul;
    p_hi = 0x0000000004000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000004000000ul 0x000020000000000dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x00007ffffffff353ul;
    p_hi = 0x0000000004000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000004000000ul 0x00007ffffffff353ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000020000000003dul;
    p_hi = 0x0000000008000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000008000000ul 0x000020000000003dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfff41ffffffffcc7ul;
    p_hi = 0x0000000007fffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000000007fffffful 0xfff41ffffffffcc7ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000400000000033ul;
    p_hi = 0x0000000010000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000010000000ul 0x0000400000000033ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfffe7ffffffffec5ul;
    p_hi = 0x000000000ffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000000ffffffful 0xfffe7ffffffffec5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000400000000023ul;
    p_hi = 0x0000000020000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000020000000ul 0x0000400000000023ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfff6bfffffffff97ul;
    p_hi = 0x000000001ffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000001ffffffful 0xfff6bfffffffff97ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000080000000003bul;
    p_hi = 0x0000000040000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000040000000ul 0x000080000000003bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffc8fffffffffdc1ul;
    p_hi = 0x000000003ffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000003ffffffful 0xffc8fffffffffdc1ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0000800000000057ul;
    p_hi = 0x0000000080000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000080000000ul 0x0000800000000057ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xff977ffffffff691ul;
    p_hi = 0x000000007ffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000007ffffffful 0xff977ffffffff691ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0001000000000033ul;
    p_hi = 0x0000000100000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000100000000ul 0x0001000000000033ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffd9fffffffffb29ul;
    p_hi = 0x00000000fffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000000fffffffful 0xffd9fffffffffb29ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000100000000005bul;
    p_hi = 0x0000000200000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000200000000ul 0x000100000000005bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffcefffffffff019ul;
    p_hi = 0x00000001fffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000001fffffffful 0xffcefffffffff019ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0002000000000041ul;
    p_hi = 0x0000000400000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000400000000ul 0x0002000000000041ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xffe7ffffffffea2bul;
    p_hi = 0x00000003fffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000003fffffffful 0xffe7ffffffffea2bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x000200000000001bul;
    p_hi = 0x0000000800000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000000800000000ul 0x000200000000001bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xff29ffffffffee99ul;
    p_hi = 0x00000007fffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00000007fffffffful 0xff29ffffffffee99ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x00040000000000a1ul;
    p_hi = 0x0000001000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000001000000000ul 0x00040000000000a1ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x006ffffffffffa33ul;
    p_hi = 0x0000001000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000001000000000ul 0x006ffffffffffa33ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0004000000000049ul;
    p_hi = 0x0000002000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000002000000000ul 0x0004000000000049ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xff7bfffffffffdc9ul;
    p_hi = 0x0000001ffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000001ffffffffful 0xff7bfffffffffdc9ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0008000000000007ul;
    p_hi = 0x0000004000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000004000000000ul 0x0008000000000007ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfc9ffffffffff56bul;
    p_hi = 0x0000003ffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000003ffffffffful 0xfc9ffffffffff56bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x00080000000000e7ul;
    p_hi = 0x0000008000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000008000000000ul 0x00080000000000e7ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf897fffffffff56bul;
    p_hi = 0x0000007ffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000007ffffffffful 0xf897fffffffff56bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0010000000000009ul;
    p_hi = 0x0000010000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000010000000000ul 0x0010000000000009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfe5ffffffffffc25ul;
    p_hi = 0x000000fffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000000fffffffffful 0xfe5ffffffffffc25ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0010000000000001ul;
    p_hi = 0x0000020000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000020000000000ul 0x0010000000000001ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfa6fffffffffff15ul;
    p_hi = 0x000001fffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000001fffffffffful 0xfa6fffffffffff15ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0020000000000061ul;
    p_hi = 0x0000040000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000040000000000ul 0x0020000000000061ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf2bffffffffffdd5ul;
    p_hi = 0x000003fffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000003fffffffffful 0xf2bffffffffffdd5ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0020000000000039ul;
    p_hi = 0x0000080000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000080000000000ul 0x0020000000000039ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf81fffffffffbb0ful;
    p_hi = 0x000007fffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x000007fffffffffful 0xf81fffffffffbb0ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x004000000000001dul;
    p_hi = 0x0000100000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000100000000000ul 0x004000000000001dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x1f7fffffffffeb81ul;
    p_hi = 0x0000100000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0000100000000000ul 0x1f7fffffffffeb81ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0040000000000085ul;
    p_hi = 0x0000200000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000200000000000ul 0x0040000000000085ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf03fffffffffff9dul;
    p_hi = 0x00001ffffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00001ffffffffffful 0xf03fffffffffff9dul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0080000000000007ul;
    p_hi = 0x0000400000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000400000000000ul 0x0080000000000007ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xe5ffffffffffff5bul;
    p_hi = 0x00003ffffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00003ffffffffffful 0xe5ffffffffffff5bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x008000000000003dul;
    p_hi = 0x0000800000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0000800000000000ul 0x008000000000003dul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf17fffffffffee99ul;
    p_hi = 0x00007ffffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x00007ffffffffffful 0xf17fffffffffee99ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x010000000000000bul;
    p_hi = 0x0001000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0001000000000000ul 0x010000000000000bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x4bfffffffffffe6bul;
    p_hi = 0x0001000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0001000000000000ul 0x4bfffffffffffe6bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0100000000000017ul;
    p_hi = 0x0002000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0002000000000000ul 0x0100000000000017ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xfeffffffffffffd3ul;
    p_hi = 0x0001fffffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0001fffffffffffful 0xfeffffffffffffd3ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0200000000000019ul;
    p_hi = 0x0004000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0004000000000000ul 0x0200000000000019ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xf7ffffffffffff8bul;
    p_hi = 0x0003fffffffffffful;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0003fffffffffffful 0xf7ffffffffffff8bul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0200000000000009ul;
    p_hi = 0x0008000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0008000000000000ul 0x0200000000000009ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x55fffffffffffc7ful;
    p_hi = 0x0008000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0008000000000000ul 0x55fffffffffffc7ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x040000000000000bul;
    p_hi = 0x0010000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0010000000000000ul 0x040000000000000bul  \n");
        return -1;
    }

    // not prime
    p_lo = 0xa7fffffffffff8b9ul;
    p_hi = 0x0010000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0010000000000000ul 0xa7fffffffffff8b9ul  \n");
        return -1;
    }

    // prime
    p_lo = 0x0400000000000023ul;
    p_hi = 0x0020000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0020000000000000ul 0x0400000000000023ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x33fffffffffff22ful;
    p_hi = 0x0020000000000001ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0020000000000001ul 0x33fffffffffff22ful  \n");
        return -1;
    }

    // prime
    p_lo = 0x0800000000000013ul;
    p_hi = 0x0040000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0040000000000000ul 0x0800000000000013ul  \n");
        return -1;
    }

    // not prime
    p_lo = 0x5fffffffffffe3dbul;
    p_hi = 0x0040000000000002ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x0040000000000002ul 0x5fffffffffffe3dbul  \n");
        return -1;
    }

    // prime
    p_lo = 0x080000000000001ful;
    p_hi = 0x0080000000000000ul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (!b)
    {
        printf("  0x0080000000000000ul 0x080000000000001ful  \n");
        return -1;
    }

    // not prime
    p_lo = 0x97fffffffffff8e9ul;
    p_hi = 0x007ffffffffffffdul;
    p = (p_hi << 64) + p_lo;
    b = uint128_is_prime_nist(p);
    if (b)
    {
        printf("  0x007ffffffffffffdul 0x97fffffffffff8e9ul  \n");
        return -1;
    }

    return 0;
}
