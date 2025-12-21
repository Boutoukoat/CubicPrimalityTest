
static int self_test_isqrt_64(void)
{
    printf("Integer square root\n");

    // uint64_isqrt() is limited to 63 bits
    uint64_t test[] = {0x12,         0x21,         0x1234,       0x4321,       0x12345,      0x54321, 0x12345678ul,
                       0x87654321ul, 0x7ffffffeul, 0x7ffffffful, 0xbffffffeul, 0xbffffffful, 0};

    for (unsigned j = 0; test[j] != 0; j++)
    {
        uint64_t s, s2, t;
        s = test[j];
        s2 = s * s;

        t = uint64_isqrt(s2);
        if (t != s)
        {
            return -1;
        }

        t = uint64_isqrt(s2 + 1);
        if (t != s)
        {
            return -1;
        }

        t = uint64_isqrt(s2 - 1);
        if (t + 1 != s)
        {
            return -1;
        }
    }
    return 0;
}
