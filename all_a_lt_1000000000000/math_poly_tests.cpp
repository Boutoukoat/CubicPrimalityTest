
static int self_test_poly_64(void)
{
	printf("Quadratic polynomial test\n");
    // prime
    uint64_poly b = {1, 0};
    uint64_t n = 101;
    uint64_t a = 7;
    uint64_t inv3 = uint64_mod_inv(3, n);
    uint64_t c = mul_mod(mul_mod(a, a - 1, n), inv3, n);
    uint64_poly r = uint64_polypow(b, n + 1, a, c, n);
    if (r.u != 0 || r.v != c)
    {
	    return -1;
    }

    // known pseudoprime
    n = 1521734597;
    uint64_t k = 396738124;
    a = (7 + k * (k - 1)) % n;
    inv3 = uint64_mod_inv(3, n);
    c = mul_mod(mul_mod(a, a - 1, n), inv3, n);
    r = uint64_polypow(b, n + 1, a, c, n);
    if (r.u != 0 || r.v != c)
    {
	    return -1;
    }

    // composite
    n = 1521734597;
    k = 1;
    a = (7 + k * (k - 1)) % n;
    inv3 = uint64_mod_inv(3, n);
    c = mul_mod(mul_mod(a, a - 1, n), inv3, n);
    r = uint64_polypow(b, n + 1, a, c, n);
    if (r.u == 0 && r.v == c)
    {
	    return -1;
    }

    return 0;
}
