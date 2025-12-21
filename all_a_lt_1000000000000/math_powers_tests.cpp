

static int self_test_powers_64(void)
{
    uint64_t s, r, t;
    bool b;

    printf("Perfect square ...\n");
    b = is_perfect_square(6);
    if (b)
        return -1;
    b = is_perfect_square(64);
    if (!b)
        return -1;
    b = is_perfect_square(27);
    if (b)
        return -1;
    b = is_perfect_square(0x1002001);
    if (!b)
        return -1;
    b = is_perfect_square(0x1002000);
    if (b)
        return -1;
    b = is_perfect_square(0x1002002);
    if (b)
        return -1;

    printf("Perfect cube ...\n");
    b = is_perfect_cube(6);
    if (b)
        return -1;
    b = is_perfect_cube(64);
    if (!b)
        return -1;
    b = is_perfect_cube(81);
    if (b)
        return -1;
    b = is_perfect_cube(0x1003003001);
    if (!b)
        return -1;
    b = is_perfect_cube(0x1003003000);
    if (b)
        return -1;
    b = is_perfect_cube(0x1003003002);
    if (b)
        return -1;

    printf("Perfect sursolid ...\n");
    b = is_perfect_sursolid(6);
    if (b)
        return -1;
    b = is_perfect_sursolid(64 * 16);
    if (!b)
        return -1;
    b = is_perfect_sursolid(81);
    if (b)
        return -1;
    b = is_perfect_sursolid(0x100500A00A005001ull);
    if (!b)
        return -1;
    b = is_perfect_sursolid(0x100500A00A005000ull);
    if (b)
        return -1;
    b = is_perfect_sursolid(0x100500A00A005002ull);
    if (b)
        return -1;

    printf("Semiprime numbers (aka biprimes)\n");
    unsigned semiprime[] = {6,   10,  14,  15,  21,  22,  26,  33,  34,  35,  38,  39,  46,  51,  55,  57,
                            58,  62,  65,  69,  74,  77,  82,  85,  86,  87,  91,  93,  94,  95,  106, 111,
                            115, 118, 119, 122, 123, 129, 133, 134, 141, 142, 143, 145, 146, 155, 158, 159,
                            161, 166, 177, 178, 183, 185, 187, 194, 201, 202, 203, 205, 0};
    unsigned psj = 3;
    for (unsigned psi = 0; semiprime[psi]; psi++)
    {
        factor_v f;
        while (psj < semiprime[psi])
        {
            f.clear();
            uint64_all_factors(f, psj);
            b = is_semiprime(f);
            if (b)
                return -1;
            psj++;
        }
        f.clear();
        uint64_all_factors(f, psj);
        b = is_semiprime(f);
        if (!b)
            return -1;
        psj++;
    }

    printf("Sphenic numbers (aka triprimes)\n");
    unsigned sphenic[] = {30,  42,  66,  70,  78,  102, 105, 110, 114, 130, 138, 154, 165, 170, 174, 182, 186, 190,
                          195, 222, 230, 231, 238, 246, 255, 258, 266, 273, 282, 285, 286, 290, 310, 318, 322, 345,
                          354, 357, 366, 370, 374, 385, 399, 402, 406, 410, 418, 426, 429, 430, 434, 435, 438, 0};
    psj = 3;
    for (unsigned psi = 0; sphenic[psi]; psi++)
    {
        factor_v f;
        while (psj < sphenic[psi])
        {
            f.clear();
            uint64_all_factors(f, psj);
            b = is_sphenic(f);
            if (b)
                return -1;
            psj++;
        }
        f.clear();
        uint64_all_factors(f, psj);
        b = is_sphenic(f);
        if (!b)
            return -1;
        psj++;
    }

    printf("Perfect power ...\n");
    unsigned perfect_powers[] = {1,   4,   8,   9,   16,  25,  27,  32,  36,  49,  64,   81,   100,  121, 125,
                                 128, 144, 169, 196, 216, 225, 243, 256, 289, 324, 343,  361,  400,  441, 484,
                                 512, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1000, 1024, 1089, 0};
    unsigned ppj = 1;
    for (unsigned ppi = 0; perfect_powers[ppi]; ppi++)
    {
        while (ppj < perfect_powers[ppi])
        {
            b = is_perfect_power(ppj);
            if (b)
                return -1;
            ppj++;
        }
        b = is_perfect_power(ppj);
        if (!b)
            return -1;
        ppj++;
    }
    b = is_perfect_power(0x100500A00A005001ull);
    if (!b)
        return -1;
    b = is_perfect_power(0x100500A00A005000ull);
    if (b)
        return -1;
    b = is_perfect_power(0x100500A00A005002ull);
    if (b)
        return -1;

    printf("Power ...\n");
    barrett_t bt;
    barrett_precompute(&bt, 197);
    r = pow_mod(2, 0xfedc, 197);
    s = pow2_mod(0xfedc, 197);
    t = barrett_pow_mod(2, 0xfedc, bt);
    if (r != 182 || r != s || r != t)
        return -1;
    r = pow_mod(2, 0x8765, 197);
    s = pow2_mod(0x8765, 197);
    t = barrett_pow_mod(2, 0x8765, bt);
    if (r != 103 || r != s || r != t)
        return -1;
    r = pow_mod(2, 0x81, 197);
    if (r != 153)
        return -1;
    r = pow_mod(2, 0x80, 197);
    s = pow2_mod(0x80, 197);
    if (r != 175 || r != s)
        return -1;
    r = pow_mod(2, 0x41, 197);
    if (r != 122)
        return -1;
    r = pow_mod(2, 0x40, 197);
    s = pow2_mod(0x40, 197);
    if (r != 61 || r != s)
        return -1;
    r = pow_mod(2, 0x22, 197);
    if (r != 155)
        return -1;
    r = pow_mod(2, 0x21, 197);
    if (r != 176)
        return -1;
    r = pow_mod(2, 0x20, 197);
    s = pow2_mod(0x20, 197);
    if (r != 88 || r != s)
        return -1;
    r = pow_mod(2, 0x1f, 197);
    if (r != 44)
        return -1;
    r = pow_mod(2, 0x1e, 197);
    if (r != 22)
        return -1;
    r = pow_mod(2, 0x1d, 197);
    s = pow2_mod(0x1d, 197);
    if (r != 11 || r != s)
        return -1;
    r = pow_mod(2, 0x1c, 197);
    s = pow2_mod(0x1c, 197);
    if (r != 104 || r != s)
        return -1;

    r = pow_mod(3, 0xaa55, 197);
    t = barrett_pow_mod(3, 0xaa55, bt);
    if (r != 0xa7 || r != t)
        return -1;

    barrett_precompute(&bt, 3725);
    r = pow_mod(3422, 252, 3725);
    t = barrett_pow_mod(3422, 252, bt);
    if (r != 1116 || r != t)
        return -1;

    // passed
    return 0;
}
