
static int self_test_shift_64(void)
{
	uint64_t r, s, t;
    printf("Modular shift ...\n");

    s = 1;
    t = 65535;
    r = shift_mod(s, 16, t);
    if (r != 1)
    {
	    return -1;
    }
    r = shift_mod(s, 32, t);
    if (r != 1)
    {
	    return -1;
    }
    r = shift_mod(s, 48, t);
    if (r != 1)
    {
	    return -1;
    }

    s = 3ull << 47;
    t = (1ull << 60) - 1;
    r = shift_mod(s, 13, t);
    if (r != 3)
    {
	    return -1;
    }
    r = shift_mod(s, 23, t);
    if (r != 3072)
    {
	    return -1;
    }
    assert(r == 3072);
    r = shift_mod(s, 33, t);
    if (r != 3145728)
    {
	    return -1;
    }

    return 0;
}

static int self_test_shift_128(void)
{
    printf("Modular shift (128 bits) ...\n");

            uint128_t r, s, t;
    s = 1;
    t = 65535;
    r = uint128_shift_mod(s, 16, t);
    if (r != 1)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 32, t);
    if (r != 1)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 48, t);
    if (r != 1)
    {
            return -1;
    }

    s = 3ull << 47;
    t = (1ull << 60) - 1;
    r = uint128_shift_mod(s, 13, t);
    if (r != 3)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 23, t);
    if (r != 3072)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 33, t);
    if (r != 3145728)
    {
            return -1;
    }

    s = 3ull << 47;
    t = 1;
    t <<= 80;
    t -= 1;
    r = uint128_shift_mod(s, 33, t);
    if (r != 3)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 113, t);
    if (r != 3)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 193, t);
    if (r != 3)
    {
            return -1;
    }

    s = 5;
    s <<= 80;
    t = 1;
    t <<= 80;
    t -= 1;
    r = uint128_shift_mod(s, 0, t);
    if (r != 5)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 80, t);
    if (r != 5)
    {
            return -1;
    }
    r = uint128_shift_mod(s, 160, t);
    if (r != 5)
    {
            return -1;
    }


    return 0;
}
