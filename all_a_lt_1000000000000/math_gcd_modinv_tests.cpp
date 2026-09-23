
static int self_test_gcd_modinv_64(void)
{
    uint64_t s, r, t;

    printf("Gcd ...\n");

    s = 12;
    t = 15;
    r = uint64_gcd(s, t);
    if (r != 3)
    {
        printf("gcd failed\n");
        return -1;
    }
    s = 12;
    t = 30;
    r = uint64_gcd(s, t);
    if (r != 6)
    {
        printf("gcd failed\n");
        return -1;
    }

    printf("Modinv ...\n");

    s = 11;
    t = 15;
    r = uint64_mod_inv(s, t);
    if (r != 11)
    {
        printf("modinv failed 1\n");
        return -1;
    }
    s = 12;
    t = 31;
    r = uint64_mod_inv(s, t);
    if (r != 13)
    {
        printf("modinv failed 2\n");
        return -1;
    }
    s = 1234567;
    t = 87654321;
    r = uint64_mod_inv(s, t);
    if (r != 75327931)
    {
        printf("modinv failed 3\n");
        return -1;
    }
    s = 11;
    t = 99;
    r = uint64_mod_inv(s, t);
    if (r != 0)
    {
        printf("modinv failed 4\n");
        return -1;
    }

    t = 1234567;
    s = 87654321 + t * t;
    r = uint64_mod_inv(s % t, t);
    if (r != 173611)
    {
        printf("modinv failed 5\n");
        return -1;
    }
    r = uint64_mod_inv(s, t);
    if (r != 173611)
    {
        printf("modinv failed 6\n");
        return -1;
    }

    t = 87654321;
    s = 1234567 + t * (t + 1);
    r = uint64_mod_inv(s % t, t);
    if (r != 75327931)
    {
        printf("modinv failed 7\n");
        return -1;
    }
    r = uint64_mod_inv(s, t);
    if (r != 75327931)
    {
        printf("modinv failed 8\n");
        return -1;
    }

    printf("Modinv pow2 ...\n");

    s = 1234567;
    r = uint64_mod_inv_pow2(s);
    if (r != 3481187596769554633ul)
    {
        printf("modinv pow2 failed\n");
        return -1;
    }

    r = uint64_mod_inv_pow2(r);
    if (r != 1234567)
    {
        printf("modinv pow2 failed\n");
        return -1;
    }

    // pass
    return 0;
}
