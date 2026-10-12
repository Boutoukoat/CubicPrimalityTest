
#include "math_utils.cpp"
#include "math_poly_utils.cpp"

#include "math_barrett_tests.cpp"
#include "math_factor_tests.cpp"
#include "math_gcd_modinv_tests.cpp"
#include "math_isqrt_tests.cpp"
#include "math_jacobi_tests.cpp"
#include "math_mod_tests.cpp"
#include "math_phi_tests.cpp"
#include "math_powers_tests.cpp"
#include "math_prime_tests.cpp"
#include "math_shift_tests.cpp"
#include "math_poly_tests.cpp"

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

    if (self_test_gcd_modinv_64() != 0)
    {
        printf("gcd or modinv tests failed\n");
        return -1;
    }

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

    if (self_test_isqrt_64() != 0)
    {
        printf("Integer square root failed\n");
        return -1;
    }

    if (self_test_factor_64() != 0)
    {
        printf("Factor search failed\n");
        return -1;
    }

    if (self_test_poly_64() != 0)
    {
        printf("Polynomial tests failed\n");
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
