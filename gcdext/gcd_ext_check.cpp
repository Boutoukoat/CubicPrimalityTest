// -----------------------------------------------------------------------
// verify details of the cubic test
//
// {gettime();l=10000;forprime(p=7,1000000000,if(p>l,print([l,gettime()]);l+=10000);
// forstep(Q=5,p-1,2,if(Q%3!=0,e=gcdext(p^3,Q^3);if((e[1]+e[2])%(p^3-1)==1,print([p,q,e])))));}
//
// -----------------------------------------------------------------------

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned __int128 uint128_t;
typedef signed __int128 int128_t;

// (u << 64 + v) mod n
static inline uint64_t uint64_long_mod(uint64_t u, uint64_t v, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(u), "1"(v), "r"(n));
    return r;
#else
    uint128_t t = ((uint128_t)u << 64) + v;
    return t % n;
#endif
}

// (u) mod n
static inline uint64_t uint128_long_mod(uint128_t u, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r = (uint64_t)(u >> 64), a = (uint64_t)u;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    return u % n;
#endif
}

static inline uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t tmp = (uint128_t)a * b;
    tmp %= n;
    return tmp;
#endif
}

static inline uint64_t square_mod(uint64_t a, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %2" : "=d"(r), "=a"(a) : "1"(a));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t tmp = (uint128_t)a * a;
    tmp %= n;
    return tmp;
#endif
}

// (u << s) mod n
static inline uint64_t shift_mod(uint64_t u, uint64_t s, uint64_t n)
{
#ifdef __x86_164__
    uint64_t r;
    asm("xor %0, %0\n shldq %b3, %1, %0\n shlxq %3, %1, %%rax\n divq %2"
        : "=&d"(r)
        : "r"(u), "r"(n), "c"(s)
        : "flags", "%rax");
    return r;
#else
    uint128_t t = (uint128_t)u;
    t <<= s;
    return t % n;
#endif
}

// count leading zeroed bits
static inline uint64_t uint64_lzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("lzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_clzll(a);
#endif
}

// count trailing zeroed bits
static inline uint64_t uint64_tzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("tzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_ctzll(a);
#endif
}

// simple floor(log_2) function
// log(1) = 0
// log(2) = 1
// log(3) = 1
// log(4) = 2 ...
static inline uint64_t uint64_log_2(uint64_t a)
{
    return 63 - uint64_lzcnt(a);
}

static int uint64_jacobi(uint64_t x, uint64_t y)
{
    // assert((y & 1) == 1);
    if (y == 1 || x == 1)
    {
        return 1;
    }

    if (x == 2)
    {
        // char j[4] = { -1,-1,1,1};
        // return j[((y - 3) >> 1) % 4];
        return ((y + 2) & 4) ? -1 : 1;
    }
    if (x == 3)
    {
        char j[6] = {0, (char)-1, (char)-1, 0, 1, 1};
        return j[((y - 3) >> 1) % 6];
    }
    if (x == 5)
    {
        char j[5] = {(char)-1, 0, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 5];
    }
    if (x == 7)
    {
        char j[14] = {1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 0, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 14];
    }
    if (x == 11)
    {
        char j[22] = {(char)-1, 1,        1,        1, 0, (char)-1, (char)-1, (char)-1, 1, (char)-1, (char)-1, 1,
                      (char)-1, (char)-1, (char)-1, 0, 1, 1,        1,        (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 22];
    }
    if (x == 13)
    {
        char j[13] = {1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 13];
    }
    if (x == 17)
    {
        char j[17] = {(char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        1, 0, 1,
                      1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 17];
    }
    if (x == 19)
    {
        char j[19] = {1,        1, (char)-1, 1,        (char)-1, (char)-1, 1,        1,        0,       (char)-1,
                      (char)-1, 1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 38;
        return t >= 19 ? -j[t - 19] : j[t];
    }
    if (x == 23)
    {
        char j[23] = {(char)-1, (char)-1, 1,        1, 1,        1,        1,        (char)-1,
                      1,        (char)-1, 0,        1, (char)-1, 1,        (char)-1, (char)-1,
                      (char)-1, (char)-1, (char)-1, 1, 1,        (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 46;
        return t >= 23 ? -j[t - 23] : j[t];
    }
    if (x == 29)
    {
        char j[29] = {(char)-1, 1, 1,        1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1,
                      1,        1, (char)-1, 0, (char)-1, 1, 1,        (char)-1, (char)-1, (char)-1,
                      (char)-1, 1, (char)-1, 1, 1,        1, (char)-1, 1,        1};
        return j[((y - 3) >> 1) % 29];
    }
    if (x == 31)
    {
        char j[31] = {1,        1, (char)-1, 1,        1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                      1,        1, (char)-1, 0,        1, (char)-1, (char)-1, (char)-1, 1,        1,        1,
                      (char)-1, 1, (char)-1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 62;
        return t >= 31 ? -j[t - 31] : j[t];
    }

    int t = 1;
    uint64_t a = x;
    uint64_t n = y;
    unsigned v = n & 7;
    unsigned c = (v == 3) || (v == 5);
    while (a)
    {
        v = __builtin_ctzll(a);
        a >>= v;
        t = (c & (v & 1)) ? -t : t;

        if (a < n)
        {
            uint64_t tmp = a;
            a = n;
            n = tmp;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            c = (v == 3) || (v == 5);
        }

        a -= n;
    }

    return (n == 1) ? t : 0;
}

// return smallest factor of n < 157, or 1 if none is found.
uint64_t small_factor(uint64_t n)
{
    if (n <= 152)
    {
        uint8_t stooopid_factor_table[] = {
            1, 1, 1, 1, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1,  2, 3, 2, 1, 2, 5, 2, 3, 2,  1, 2,
            1, 2, 3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 7, 2,  3, 2, 1, 2, 5, 2, 3, 2, 1,  2, 1,
            2, 3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 7, 2, 1, 2, 3,  2, 1, 2, 5, 2, 3, 2, 1, 2,  7, 2,
            3, 2, 5, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2,  1, 2, 5, 2, 3, 2, 7, 2, 11, 2, 3,
            2, 5, 2, 1, 2, 3, 2, 1, 2, 7, 2, 3, 2, 1, 2, 1, 2, 3, 2, 11, 2, 5, 2, 3, 2, 1, 2, 1, 2};
        return stooopid_factor_table[n];
    }
    if (!(n & 1))
        return 2;
    if ((uint64_t)(n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return 3;
    if ((uint64_t)(n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return 5;
    if ((uint64_t)(n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return 7;
    if ((uint64_t)(n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return 11;
    if ((uint64_t)(n * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return 13;
    if ((uint64_t)(n * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return 17;
    if ((uint64_t)(n * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return 19;
    if ((uint64_t)(n * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return 23;
    if ((uint64_t)(n * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return 29;
    if ((uint64_t)(n * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return 31;
    if (n < 37 * 37)
        return 1; // prime
    if ((uint64_t)(n * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return 37;
    if ((uint64_t)(n * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return 41;
    if ((uint64_t)(n * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return 43;
    if ((uint64_t)(n * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return 47;
    if ((uint64_t)(n * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return 53;
    if ((uint64_t)(n * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return 59;
    if ((uint64_t)(n * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return 61;
    if ((uint64_t)(n * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return 67;
    if ((uint64_t)(n * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return 71;
    if ((uint64_t)(n * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return 73;
    if ((uint64_t)(n * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return 79;
    if ((uint64_t)(n * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return 83;
    if ((uint64_t)(n * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return 89;
    if ((uint64_t)(n * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return 97;
    if (n < 101 * 101)
        return 1; // prime
    if ((uint64_t)(n * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return 101;
    if ((uint64_t)(n * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return 103;
    if ((uint64_t)(n * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return 107;
    if ((uint64_t)(n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return 109;
    if ((uint64_t)(n * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return 113;
    if ((uint64_t)(n * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return 127;
    if ((uint64_t)(n * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return 131;
    if ((uint64_t)(n * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return 137;
    if ((uint64_t)(n * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return 139;
    if ((uint64_t)(n * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return 149;
    if ((uint64_t)(n * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return 151;
    return 1; // no small factor
}

// modular exponentiation a^e mod m
uint64_t pow2_mod(uint64_t e, uint64_t m)
{
    uint64_t n = uint64_log_2(e);
    uint64_t s = (n >= 5) ? 5 : n;
    n -= s;
    uint64_t mask = e >> n;
    uint64_t result = shift_mod(1ull, mask, m);
    while (n >= 6)
    {
        n -= 6;
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        result = square_mod(result, m);
        mask = (e >> n) & 0x3f;
        result = shift_mod(result, mask, m);
    }
    while (n > 0)
    {
        n -= 1;
        result = square_mod(result, m);
        if ((e >> n) & 1)
        {
            result <<= 1;
            result -= (result >= m) ? m : 0;
        }
    }
    return result;
}

// modular exponentiation a^e mod m
uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t n = e;
    uint64_t s = a;
    uint64_t result = 1;
    while (n)
    {
        if (n & 1)
            result = mul_mod(result, s, m);
        s = square_mod(s, m);
        n >>= 1;
    }
    return result;
}

// MR strong test
bool witness(uint64_t n, uint64_t s, uint64_t d, uint64_t a)
{
    uint64_t x, y;
    if (n == a)
        return true;

    if (a == 2)
    {
        x = pow2_mod(d, n);
    }
    else
    {
        x = pow_mod(a, d, n);
    }

    while (s)
    {
        y = square_mod(x, n);
        if (y == 1 && x != 1 && x != n - 1)
            return false;
        x = y;
        --s;
    }
    if (y != 1)
        return false;
    return true;
}

// deterministic primality test for n < 2^64
bool uint64_is_prime(uint64_t n)
{
    uint64_t d = n / 2;
    uint64_t s = uint64_tzcnt(d);
    d >>= s++;

    if (n < 1373653)
        return witness(n, s, d, 2) && witness(n, s, d, 3);
    if (n < 9080191)
        return witness(n, s, d, 31) && witness(n, s, d, 73);
    if (n < 4759123141)
        return witness(n, s, d, 2) && witness(n, s, d, 7) && witness(n, s, d, 61);
    if (n < 1122004669633)
        return witness(n, s, d, 2) && witness(n, s, d, 13) && witness(n, s, d, 23) && witness(n, s, d, 1662803);
    if (n < 2152302898747)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11);
    if (n < 3474749660383)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13);
    if (n < 341550071728321)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17);
    if (n < 3825123056546413051)
        return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
               witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19) &&
               witness(n, s, d, 23);
    // n < 318665857834031151167461
    return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7) &&
           witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19) &&
           witness(n, s, d, 23) && witness(n, s, d, 29) && witness(n, s, d, 31) && witness(n, s, d, 37);
}

// binary gcd
uint64_t gcd(uint64_t u, uint64_t v)
{
    uint64_t t, k;

    if (u < v)
    {
        t = u;
        u = v;
        v = t;
    }
    if (v == 0)
        return u;

    k = uint64_tzcnt(u | v);
    u >>= uint64_tzcnt(u);

    do
    {
        v >>= uint64_tzcnt(v);

        if (u > v)
        {
            t = u;
            u = v;
            v = t;
        }
        v -= u;
    } while (v);
    return u << k;
}

// binary modular inverse 1/x mod m, with m odd and x < m, x and m coprime
uint64_t mod_inv(uint64_t x, uint64_t m)
{
    if (x < 2)
        return x;
    if (m < 3)
        return 0;
    uint64_t a = x, b = m, u = 1, v = 0;
    while (a != 0)
    {
        unsigned za = uint64_tzcnt(a);
        a >>= za;
        while (za--)
        {
            u += (u & 1) ? m : 0;
            u >>= 1;
        }
        if (a < b)
        {
            uint64_t t = a;
            uint64_t s = u;
            a = b;
            u = v;
            b = t;
            v = s;
        }
        a -= b;
        u = (u >= v) ? u - v : u + m - v;
    }
    return b == 1 ? v : 0;
}

// miserably slow textbook version of Extended Euclidean algorithm
void int128_gcdext(uint128_t a, uint128_t b, int128_t &ea, int128_t &eb, uint128_t &eg)
{
    int128_t r_prev = a;
    int128_t r_next = b;
    int128_t s_prev = 1;
    int128_t s_next = 0;
    int128_t t_prev = 0;
    int128_t t_next = 1;
    int128_t q, tmp;

    while( r_next != 0 )
    {
        q = r_prev / r_next;

	tmp = r_prev - q * r_next;
	r_prev = r_next;
	r_next = tmp;

	tmp = s_prev - q * s_next;
	s_prev = s_next;
	s_next = tmp;

	tmp = t_prev - q * t_next;
	t_prev = t_next;
	t_next = tmp;
    }
    eg = r_prev;
    ea = s_prev;
    eb = t_prev;
}

// conversion of a large number into a basis-10 string.
// returns the output string length.
unsigned uint128_sprint(char *ptr, uint128_t x)
{
    char buff[256];
    char *pb = buff;
    char *pt = ptr;

    if (x == 0)
    {
        *pb++ = '0';
    }
    while (x)
    {
        *pb++ = '0' + (char)(x % 10);
        x /= 10;
    }
    while (pb > buff)
    {
        *(pt++) = *(--pb);
    }
    *pt = 0;
    return pt - ptr;
}

void self_tests(void)
{
	int128_t ea, eb;
	uint128_t eg;
	int128_gcdext(12, 15, ea, eb, eg);
	assert(ea == -1);
	assert(eb == 1);
	assert(eg == 3);

	int128_gcdext(101*101*101, 103*103*103, ea, eb, eg);
	assert(ea == 407764);
	assert(eb == -384469);
	assert(eg == 1);

	char b[128];
	unsigned u;
	u = uint128_sprint(b, 123456);
	assert(strcmp(b, "123456") == 0);
	assert(u = strlen(b));
	assert(u = 10);
	u = uint128_sprint(b, 98765432109876543ull);
	assert(strcmp(b, "98765432109876543") == 0);
	assert(u = strlen(b));
	assert(u = 17);
}

int main(int argc, char **argv)
{
    uint64_t p_max = 1ull << 21;
    uint64_t l_max = 10000;

    self_tests();

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-p"))
        {
            p_max = strtoull(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-l"))
        {
            l_max = strtoull(argv[++i], 0, 0);
            continue;
        }
        printf("-p : maximum prime to search for < 10^12 approx\n");
        printf("-l : display interval e.g. 10000\n");
    }

    time_t t0 = time(NULL);
    uint64_t l = l_max;
    // iterate on multiples of 6k+/-1
    for (uint64_t p = 7, dp = 4; p < p_max; p += dp, dp = 6 - dp)
    {
        if (uint64_is_prime(p))
        {
		// display some progress from time to time
            if (p >= l)
            {
                char buff[256];
                char *ptr = buff;
                time_t t1 = time(NULL);
                *ptr++ = '[';
                ptr += uint128_sprint(ptr, p);
                *ptr++ = ',';
                ptr += uint128_sprint(ptr, t1 - t0);
                *ptr++ = ']';
                *ptr = 0;
                printf("%s\n", buff);
                fflush(stdout);
                t0 = t1;
		l += l_max;
            }

            uint128_t p3 = (uint128_t)p * p * p;
            for (uint64_t q = 5; q < p - 1; q += 2)
            {
                // q % 3 == 0
                if ((uint64_t)(q * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
                {
                    continue;
                }

                // extended euclidean algorithm
                // ep * p^3 + eq * q^3 = eg
                //
                // todo : investigate gcd(p^3, q^3) = gcd(p, q)^3
                //
                uint128_t q3 = (uint128_t)q * q * q;
                int128_t ep, eq, e;
		uint128_t eg;
                int128_gcdext(p3, q3, ep, eq, eg);

                e = ep + eq;
                if (e < 0)
		{
                    e = p3 - 1 - e;
		}
                e %= p3 - 1;
                if (e == 1)
                {
                    char buff[256];
                    char *ptr = buff;
                    ptr += uint128_sprint(ptr, p);
                    *ptr++ = ' ';
                    ptr += uint128_sprint(ptr, q);
                    *ptr++ = ' ';
                    *ptr++ = '[';
                    ptr += uint128_sprint(ptr, ep);
                    *ptr++ = ',';
                    ptr += uint128_sprint(ptr, eq);
                    *ptr++ = ',';
                    ptr += uint128_sprint(ptr, eg);
                    *ptr++ = ']';
		    *ptr = 0;
                    printf("%s\n", buff);
                    fflush(stdout);
                }
            }
        }
    }

    return (0);
}
