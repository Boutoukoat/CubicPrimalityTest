
#include "math_utils.cpp"

static int self_test_64(void)
{
    uint64_t s, r, t;
    bool b;
    int j;

    printf("Modular operations ...\n");
    s = 10103;
    t = 10101;
    r = square_mod(s, t);
    assert(r == 4);
    s = 10103;
    t = 10101;
    r = mul_mod(s, s, t);
    assert(r == 4);

    s = 1;
    t = 65535;
    r = shift_mod(s, 16, t);
    assert(r == 1);
    r = shift_mod(s, 32, t);
    assert(r == 1);
    r = shift_mod(s, 48, t);
    assert(r == 1);

    s = 3ull << 47;
    t = (1ull << 60) - 1;
    r = shift_mod(s, 13, t);
    assert(r == 3);
    r = shift_mod(s, 23, t);
    assert(r == 3072);
    r = shift_mod(s, 33, t);
    assert(r == 3145728);

    printf("Barrett reduction speed-up\n");
    barrett_t u;
    s = 101;
    barrett_precompute(&u, s);
    r = barrett_mul_mod(99, 38, u);
    if (r % s != 25)
        return -1;
    if (r > 2 * s)
        return -1;

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

    printf("Gcd ...\n");
    s = 12;
    t = 15;
    r = uint64_gcd(s, t);
    if (r != 3)
        return -1;
    s = 12;
    t = 30;
    r = uint64_gcd(s, t);
    if (r != 6)
        return -1;

    printf("Jacobi ...\n");
    s = 33;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != 0)
        return -1;
    s = 34;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != -1)
        return -1;
    s = 35;
    t = 9999;
    j = uint64_jacobi(s, t);
    if (j != 1)
        return -1;

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

    printf("Modinv ...\n");
    s = 11;
    t = 15;
    r = uint64_mod_inv(s, t);
    if (r != 11)
        return -1;
    s = 12;
    t = 31;
    r = uint64_mod_inv(s, t);
    if (r != 13)
        return -1;
    s = 1234567;
    t = 87654321;
    r = uint64_mod_inv(s, t);
    if (r != 75327931)
        return -1;
    s = 11;
    t = 99;
    r = uint64_mod_inv(s, t);
    if (r != 0)
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

    printf("Known primes ...\n");
    b = uint64_is_prime(200003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime(2000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime(20000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime(2000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime(20000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = uint64_is_prime(200000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }

    printf("Isprime ...\n");
    t = 1;
    t <<= 3;
    t -= 1;
    b = uint64_is_prime(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }

    t = 101;
    b = uint64_is_prime(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }

    t = 4493;
    b = uint64_is_prime(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = uint64_is_prime(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = uint64_is_prime(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }

    printf("Self-test completed\n");
    // pass
    return 0;
}

int main(int argc, char **argv)
{
    int rc = 0;
    rc = self_test_64();
    if (rc)
    {
        printf("failed\n");
        exit(1);
    }
    printf("All tests passed\n");
    return 0;
}
