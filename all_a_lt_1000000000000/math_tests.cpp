
#include "math_utils.cpp"

#include "math_barrett_tests.cpp"
#include "math_isqrt_tests.cpp"
#include "math_jacobi_tests.cpp"
#include "math_powers_tests.cpp"

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

    if (self_test_barrett_64() != 0)
    {
        printf("Barrett tests failed\n");
        return -1;
    }

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

    if (self_test_powers_64() != 0)
    {
        printf("Powers tests failed\n");
        return -1;
    }

    if (self_test_jacobi_64() != 0)
    {
        printf("Jacobi - Kronecker tests failed\n");
        return -1;
    }

    printf("Modular Power ...\n");
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

    if (self_test_isqrt_64() != 0)
    {
        printf("Integer square root failed\n");
        return -1;
    }

    printf("Factors ...\n");

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
