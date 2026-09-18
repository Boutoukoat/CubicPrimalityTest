
#include "math_utils.cpp"

#include "math_barrett_tests.cpp"
#include "math_isqrt_tests.cpp"
#include "math_jacobi_tests.cpp"
#include "math_mod_tests.cpp"
#include "math_phi_tests.cpp"
#include "math_powers_tests.cpp"
#include "math_prime_tests.cpp"
#include "math_shift_tests.cpp"

static int self_test_64(void)
{
    uint64_t s, r, t;
    bool b;
    int j;

    if (self_test_mod_64() != 0)
    {
        printf("Modular op failed\n");
        return -1;
    }

    if (self_test_mod_128() != 0)
    {
        printf("Modular op failed\n");
        return -1;
    }

    if (self_test_shift_64() != 0)
    {
        printf("Modular shift failed\n");
        return -1;
    }

    if (self_test_shift_128() != 0)
    {
        printf("modular shift failed\n");
        return -1;
    }

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
    s = uint64_pow2_mod(0xfedc, 197);
    t = barrett_pow_mod(2, 0xfedc, bt);
    if (r != 182 || r != s || r != t)
        return -1;
    r = pow_mod(2, 0x8765, 197);
    s = uint64_pow2_mod(0x8765, 197);
    t = barrett_pow_mod(2, 0x8765, bt);
    if (r != 103 || r != s || r != t)
        return -1;
    r = pow_mod(2, 0x81, 197);
    if (r != 153)
        return -1;
    r = pow_mod(2, 0x80, 197);
    s = uint64_pow2_mod(0x80, 197);
    if (r != 175 || r != s)
        return -1;
    r = pow_mod(2, 0x41, 197);
    if (r != 122)
        return -1;
    r = pow_mod(2, 0x40, 197);
    s = uint64_pow2_mod(0x40, 197);
    if (r != 61 || r != s)
        return -1;
    r = pow_mod(2, 0x22, 197);
    if (r != 155)
        return -1;
    r = pow_mod(2, 0x21, 197);
    if (r != 176)
        return -1;
    r = pow_mod(2, 0x20, 197);
    s = uint64_pow2_mod(0x20, 197);
    if (r != 88 || r != s)
        return -1;
    r = pow_mod(2, 0x1f, 197);
    if (r != 44)
        return -1;
    r = pow_mod(2, 0x1e, 197);
    if (r != 22)
        return -1;
    r = pow_mod(2, 0x1d, 197);
    s = uint64_pow2_mod(0x1d, 197);
    if (r != 11 || r != s)
        return -1;
    r = pow_mod(2, 0x1c, 197);
    s = uint64_pow2_mod(0x1c, 197);
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

    if (self_test_prime_64() != 0)
    {
        printf("Prime 64 bits failed\n");
        return -1;
    }

    if (self_test_prime_128() != 0)
    {
        printf("Prime 128 bits failed\n");
        return -1;
    }

    if (self_test_phi_64() != 0)
    {
        printf("Totient failed\n");
        return -1;
    }

    printf("Isqrt ...\n");

    if (self_test_isqrt_64() != 0)
    {
        printf("Integer square root failed\n");
        return -1;
    }

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
