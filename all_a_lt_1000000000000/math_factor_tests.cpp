
static int self_test_factor_64(void)
{
    uint64_t t;

    printf("Factors ...\n");

    t = uint64_small_factor(11 * 13);
    if (t != 11)
    {
        printf("small factor failed\n");
        return -1;
    }

    t = uint64_small_factor(101 * 103);
    if (t != 101)
    {
        printf("small factor failed\n");
        return -1;
    }

    t = uint64_small_factor(151 * 521);
    if (t != 151)
    {
        printf("small factor failed\n");
        return -1;
    }

    t = uint64_small_factor(157 * 521);
    if (t != 1)
    {
        printf("small factor failed\n");
        return -1;
    }

    t = uint64_sqfof_factor(101 * 103);
    if (!(t == 101 || t == 103))
    {
        printf("sqfof failed\n");
        return -1;
    }

    t = uint64_sqfof_factor(157 * 157);
    if (!(t == 157))
    {
        printf("sqfof failed\n");
        return -1;
    }

    t = uint64_brent_pollard_factor(101 * 103);
    if (!(t == 101 || t == 103))
    {
        printf("pollard failed\n");
        return -1;
    }

    t = uint64_brent_pollard_factor(157 * 157);
    if (!(t == 157))
    {
        printf("pollard failed\n");
        return -1;
    }

    return 0;
}
