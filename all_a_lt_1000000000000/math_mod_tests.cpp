
static int self_test_mod_64(void)
{
    uint64_t r, s, t, u;
    printf("Modular operations ...\n");

    t = 10101;
    s = 10103;
    r = square_mod(s, t);
    if (r != 4)
    {
	    return -1;
    }

    r = mul_mod(s, s, t);
    if (r != 4)
    {
	    return -1;
    }

    t = 10101;
    s = t;
    s <<= 24;
    s += 10103;
    r = square_mod(s, t);
    if (r != 4)
    {
	    return -1;
    }

    r = mul_mod(s, s, t);
    if (r != 4)
    {
	    return -1;
    }

    t = 10101;
    s = t;
    s <<= 24;
    s += 10103;
    u = 2;
    r = square_add_mod(s, u, t);
    if (r != 4 + u)
    {
	    return -1;
    }

    r = mul_add_mod(s, s, u, t);
    if (r != 4 + u)
    {
	    return -1;
    }

    return 0;
}

static int self_test_mod_128(void)
{
    uint128_t r, s, t, u;
    printf("Modular operations (128 bits) ...\n");

    t = 10101;
    s = 10103;
    r = uint128_square_mod(s, t);
    if (r != 4)
    {
	    return -1;
    }

    r = uint128_mul_mod(s, s, t);
    if (r != 4)
    {
	    return -1;
    }

    t = 10101;
    s = t;
    s <<= 100;
    s += 10103;
    r = uint128_square_mod(s, t);
    if (r != 4)
    {
	    return -1;
    }

    r = uint128_mul_mod(s, s, t);
    if (r != 4)
    {
	    return -1;
    }

    t = 10101;
    s = t;
    s <<= 110;
    s += 10103;
    u = 2;
    r = uint128_square_add_mod(s, u, t);
    if (r != 4 + u)
    {
            return -1;
    }

    r = uint128_mul_add_mod(s, s, u, t);
    if (r != 4 + u)
    {
            return -1;
    }

    return 0;
}





