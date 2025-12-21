

static int self_test_barrett_64(void)
{
    uint64_t s, r, t;
    bool b;
    int j;

    printf("Barrett reduction sanity checks\n");
    barrett_t u;
    s = 101;
    barrett_precompute(&u, s);
    r = barrett_mul_mod(99, 38, u);
    if (r % s != 25)
        return -1;
    if (r > 2 * s)
        return -1;

    r = barrett_mul_mod(2 * s - 1, 2 * s - 1, u);
    if (r % s != 1)
        return -1;
    if (r > 2 * s)
        return -1;

    s = 0x7654321fedull;
    barrett_precompute(&u, s);
    r = barrett_mul_mod(2 * s - 1, 1, u);
    if (r % s != s - 1)
        return -1;
    if (r > 2 * s)
        return -1;

    r = barrett_mul_mod(2 * s - 1, 2 * s - 1, u);
    if (r % s != 1)
        return -1;
    if (r > 2 * s)
        return -1;

    printf("Barrett reduction validity\n");
    for (unsigned bits = 2; bits < 61; bits++)
    {
        uint64_t m_table[] = {1ull << bits,
                              (1ull << bits) + 5,
                              (3ull << (bits - 1)) - 5,
                              (3ull << (bits - 1)) + 5,
                              (2ull << bits) - 5,
                              (2ull << bits) - 1,
                              0,
                              0};
        for (unsigned im = 0; m_table[im]; im += 2)
        {
            for (uint64_t m = m_table[im]; m <= m_table[im + 1]; m += 1)
            {
                // verify modular multiplication input is <= 2 *m
                // verify modular multiplication output is <= 2 *m
                uint64_t l_table[] = {1,     2, 3,     m / 2 - 1, m / 2,     m / 2 + 1, m / 2 + 2, m - 2,
                                      m - 1, m, m + 1, m + 2,     m * 2 - 2, m * 2 - 1, m * 2,     0};
                barrett_precompute(&u, m);

                for (unsigned ia = 0; l_table[ia]; ia++)
                {
                    for (unsigned ib = 0; l_table[ib]; ib++)
                    {
                        uint64_t a = l_table[ia];
                        uint64_t b = l_table[ib];
                        uint128_t t = (uint128_t)a * b;
                        r = barrett_mul_mod(a, b, u);
                        assert(r <= 2 * m);
                        assert(t % m == r % m);
                    }
                }
            }
        }
        //	printf("Bits %u passed\n", bits);
    }

    // pass
    return 0;
}
