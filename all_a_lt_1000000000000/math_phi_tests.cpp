
static int self_test_phi_64(void)
{
    uint64_t t;
    printf("Euler Totient ...\n");

    t = 128;
    if (uint64_phi(t) != 64)
    {
        return -1;
    }

    t = 139;
    if (uint64_phi(t) != 138)
    {
        return -1;
    }

    t = 77;
    if (uint64_phi(t) != 60)
    {
        return -1;
    }

    t = 0xffff0001ul;
    if (uint64_phi(t) != 4272648192ul)
    {
        return -1;
    }

    return 0;
}
