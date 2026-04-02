#include <algorithm>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>
#include <vector>

using namespace std;

#include "bison.gmp_expr.h"

typedef unsigned __int128 uint128_t;

static unsigned target_bit_length = 256;          /* arbitrary default length */
static bool verbose = false;                      /* arbitrary choice */
static bool hex_output = true;                    /* arbitrary choice */
static uint64_t max_prime_count = 128;            /* approximate number of primes to display */
static volatile uint64_t current_prime_count = 0; /* approximate number of primes to display */

static void worker(mpz_t p);

#define MT_QUEUE_SIZE 1024

struct mod_multithread_t
{
    pthread_mutex_t __attribute__((aligned(64))) display_mutex;
    pthread_cond_t __attribute__((aligned(64))) queue_cond;
    pthread_mutex_t __attribute__((aligned(64))) queue_mutex;
    volatile unsigned long __attribute__((aligned(64))) queue_head, queue_tail;
    sem_t __attribute__((aligned(64))) queue_sem;
    mpz_t queue[MT_QUEUE_SIZE];
};

void mt_initialize(mod_multithread_t *mt)
{
    mt->display_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->queue_cond = PTHREAD_COND_INITIALIZER;
    mt->queue_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->queue_head = 0;
    mt->queue_tail = 0;
    sem_init(&mt->queue_sem, 0, MT_QUEUE_SIZE);
    for (unsigned i = 0; i < MT_QUEUE_SIZE; i++)
    {
        mpz_init(mt->queue[i]);
    }
}
void mt_uninitialize(mod_multithread_t *mt)
{
    for (unsigned i = 0; i < MT_QUEUE_SIZE; i++)
    {
        mpz_clear(mt->queue[i]);
    }
}

static bool post_prime(mod_multithread_t *mt, mpz_t p)
{
    // push a prime to the queue
    int status = sem_trywait(&mt->queue_sem);
    if (status) 
    {
	    // cannot pass the blocking condition
	    return false;
    }
    pthread_mutex_lock(&mt->queue_mutex);
    mpz_set(mt->queue[mt->queue_head++ & (MT_QUEUE_SIZE - 1)], p);
    asm("" ::: "memory"); // prevent the compiler to optimize and move memory accesses across this line
    pthread_mutex_unlock(&mt->queue_mutex);
    pthread_cond_signal(&mt->queue_cond);
    // post successful
    return true;
}

static void get_prime(mod_multithread_t *mt, mpz_t p)
{
    // get a prime from the queue
    pthread_mutex_lock(&mt->queue_mutex);
    while (mt->queue_tail == mt->queue_head)
    {
        pthread_cond_wait(&mt->queue_cond, &mt->queue_mutex);
    }
    mpz_set(p, mt->queue[mt->queue_tail++ & (MT_QUEUE_SIZE - 1)]);
    asm("" ::: "memory"); // prevent the compiler to optimize and move memory accesses across this line
    pthread_mutex_unlock(&mt->queue_mutex);
    sem_post(&mt->queue_sem);
}

static void worker(mpz_t p);
static mod_multithread_t mt;

static void display_prime(mpz_t p)
{
    pthread_mutex_lock(&mt.display_mutex);
    if (hex_output)
    {
        gmp_printf("0x%Zx\n", p);
    }
    else
    {
        gmp_printf("%Zu\n", p);
    }
    current_prime_count += 1;
    pthread_mutex_unlock(&mt.display_mutex);
}

static void display_flush(void)
{
    pthread_mutex_lock(&mt.display_mutex);
    fflush(stdout);
    pthread_mutex_unlock(&mt.display_mutex);
}

// --------------------------------------------------------------------------------------------------
//
// Pass primes to another worker, best effort
//
// --------------------------------------------------------------------------------------------------

static void next_prime(mpz_t p)
{
    unsigned bit_length = mpz_sizeinbase(p, 2);

    if (verbose || bit_length == target_bit_length)
    {
        // prime is exact size, display it
        display_prime(p);
    }
    if (bit_length < target_bit_length)
    {
        // prime is not large, continue the calculation
        if (!post_prime(&mt, p))
        {
            // queue might be full, do it
            worker(p);
        }
    }
}

// --------------------------------------------------------------------------------------------------
//
// An array of small primes
//
// --------------------------------------------------------------------------------------------------

#define SMALL_PRIMES_COUNT 30000
static unsigned int small_primes[SMALL_PRIMES_COUNT];

static void generate_small_primes(void)
{
  mpz_t sp;
    mpz_init_set_ui(sp, 3);
    for(unsigned i = 0; i < SMALL_PRIMES_COUNT; i++)
    {
            small_primes[i] = mpz_get_ui(sp);
            mpz_nextprime(sp, sp);
    }
    mpz_clear(sp);
}

// --------------------------------------------------------------------------------------------------
//
// Number theory part
//
// --------------------------------------------------------------------------------------------------

// x % (2^b -1)
static uint64_t mpz_mod_mersenne(mpz_t x, uint64_t b)
{
    uint64_t mask = (1ull << b) - 1;
    unsigned s = mpz_size(x);
    const mp_limb_t *array = mpz_limbs_read(x);
    uint128_t t, r = 0;

    if ((b & (b - 1)) == 0)
    {
        // an exact power of 2
        for (unsigned i = 0; i < s; i += 1)
        {
            r += array[i];
        }
    }
    else
    {
        // add, shift first, and reduce after
        unsigned c = 0;
        for (unsigned i = 0; i < b; i += 1)
        {
            t = 0;
            for (unsigned j = i; j < s; j += b)
            {
                t += array[j];
            }
            t = (t & mask) + (t >> b);
            t = t << c;
            r += t;
            r = (r & mask) + (r >> b);
            c += 64;
            while (c >= b)
            {
                c -= b;
            }
        }
    }
    // reduce mod 2^b - 1
    while (r >> 64)
    {
        r = (r & mask) + (r >> b);
    }
    return (uint64_t)r;
}

typedef enum sieve_e
{
    COMPOSITE_FOR_SURE,
    PRIME_FOR_SURE,
    UNDECIDED
} sieve_t;

// sieve small primes
static sieve_t uint64_composite_sieve(uint64_t a)
{
    if (a <= 152)
    {
        // return COMPOSITE_FOR_SURE for small numbers which are composite for sure , without checking further.
        sieve_t stooopid_prime_table[] = {
            PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE};
        return stooopid_prime_table[a];
    }

    // divisibility is based on Barrett modular reductions by constants, use modular multiplications
    if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return COMPOSITE_FOR_SURE; // divisible by 3
    if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return COMPOSITE_FOR_SURE; // divisible by 5
    if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return COMPOSITE_FOR_SURE; // divisible by 7
    if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return COMPOSITE_FOR_SURE; // divisible by 11
    if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return COMPOSITE_FOR_SURE; // divisible by 13
    if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return COMPOSITE_FOR_SURE; // divisible by 17
    if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return COMPOSITE_FOR_SURE; // divisible by 19
    if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return COMPOSITE_FOR_SURE; // divisible by 23
    if ((uint64_t)(a * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return COMPOSITE_FOR_SURE; // divisible by 29
    if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return COMPOSITE_FOR_SURE; // divisible by 31
    if (a < 37 * 37)
        return PRIME_FOR_SURE;
    if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return COMPOSITE_FOR_SURE; // divisible by 37
    if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return COMPOSITE_FOR_SURE; // divisible by 41
    if ((uint64_t)(a * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return COMPOSITE_FOR_SURE; // divisible by 43
    if ((uint64_t)(a * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return COMPOSITE_FOR_SURE; // divisible by 47
    if ((uint64_t)(a * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return COMPOSITE_FOR_SURE; // divisible by 53
    if ((uint64_t)(a * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return COMPOSITE_FOR_SURE; // divisible by 59
    if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return COMPOSITE_FOR_SURE; // divisible by 61
    if ((uint64_t)(a * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return COMPOSITE_FOR_SURE; // divisible by 67
    if ((uint64_t)(a * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return COMPOSITE_FOR_SURE; // divisible by 71
    if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return COMPOSITE_FOR_SURE; // divisible by 73
    if ((uint64_t)(a * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return COMPOSITE_FOR_SURE; // divisible by 79
    if ((uint64_t)(a * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return COMPOSITE_FOR_SURE; // divisible by 83
    if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return COMPOSITE_FOR_SURE; // divisible by 89
    if ((uint64_t)(a * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return COMPOSITE_FOR_SURE; // divisible by 97
    if (a < 101 * 101)
        return PRIME_FOR_SURE;
    if ((uint64_t)(a * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return COMPOSITE_FOR_SURE; // divisible by 101
    if ((uint64_t)(a * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return COMPOSITE_FOR_SURE; // divisible by 103
    if ((uint64_t)(a * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return COMPOSITE_FOR_SURE; // divisible by 107
    if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return COMPOSITE_FOR_SURE; // divisible by 109
    if ((uint64_t)(a * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return COMPOSITE_FOR_SURE; // divisible by 113
    if ((uint64_t)(a * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return COMPOSITE_FOR_SURE; // divisible by 127
    if ((uint64_t)(a * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return COMPOSITE_FOR_SURE; // divisible by 131
    if ((uint64_t)(a * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return COMPOSITE_FOR_SURE; // divisible by 137
    if ((uint64_t)(a * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return COMPOSITE_FOR_SURE; // divisible by 139
    if ((uint64_t)(a * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return COMPOSITE_FOR_SURE; // divisible by 149
    if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return COMPOSITE_FOR_SURE; // divisible by 151

    // no factor less than 157
    if (a < 157 * 157)
        return PRIME_FOR_SURE; // prime
    return UNDECIDED;
}

static sieve_t mpz_composite_sieve(mpz_t n)
{
    // detect small 64-bit numbers
    unsigned sz = mpz_sizeinbase(n, 2);
    if (sz <= 64)
    {
        uint64_t a = mpz_get_ui(n);
        sieve_t sv = uint64_composite_sieve(a);
        return sv;
    }
    else
    {
        // large number, do the modular reduction in 2 steps
        // step 1 :
        //   reduce by a multiple of small factors
        // step 2:
        //    divisibility is based on Barrett modular reductions by constants, use modular multiplications

        uint64_t a;
        // 2^60-1 is divisible by 3,5,7,11,13,31,41,61,151 ...
        a = mpz_mod_mersenne(n, 60);

        if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
            return COMPOSITE_FOR_SURE; // divisible by 3
        if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
            return COMPOSITE_FOR_SURE; // divisible by 5
        if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
            return COMPOSITE_FOR_SURE; // divisible by 7
        if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
            return COMPOSITE_FOR_SURE; // divisible by 11
        if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
            return COMPOSITE_FOR_SURE; // divisible by 13
        if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
            return COMPOSITE_FOR_SURE; // divisible by 31
        if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
            return COMPOSITE_FOR_SURE; // divisible by 41
        if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
            return COMPOSITE_FOR_SURE; // divisible by 61
        if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
            return COMPOSITE_FOR_SURE; // divisible by 151

        // 2^56-1 is divisible by 3, 5, 17, 29, 43, 113, 127, .....
        a = mpz_mod_mersenne(n, 56);
        if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
            return COMPOSITE_FOR_SURE; // divisible by 17
        if ((uint64_t)(a * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
            return COMPOSITE_FOR_SURE; // divisible by 29
        if ((uint64_t)(a * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
            return COMPOSITE_FOR_SURE; // divisible by 43
        if ((uint64_t)(a * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
            return COMPOSITE_FOR_SURE; // divisible by 113
        if ((uint64_t)(a * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
            return COMPOSITE_FOR_SURE; // divisible by 127

        // 2^36-1 is divisible by 3,5,7,19,37,73,109, ...
        a = mpz_mod_mersenne(n, 36);
        if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
            return COMPOSITE_FOR_SURE; // divisible by 19
        if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
            return COMPOSITE_FOR_SURE; // divisible by 37
        if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
            return COMPOSITE_FOR_SURE; // divisible by 73
        if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
            return COMPOSITE_FOR_SURE; // divisible by 109

        // 2^44-1 is divisible by 3, 5, 23, 89, .....
        a = mpz_mod_mersenne(n, 44);
        if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
            return COMPOSITE_FOR_SURE; // divisible by 23
        if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
            return COMPOSITE_FOR_SURE; // divisible by 89
    }

    // no trivial small factor detected
    return UNDECIDED; // might be prime
}

static void worker(mpz_t p)
{
    mpz_t tmp, n, r, two;
    mpz_inits(tmp, n, r, two, 0);
    mpz_set_ui(two, 2);
    vector<uint64_t> already_checked;

    // n = 1 + 2 * p
    while (current_prime_count < max_prime_count)
    {
	already_checked.push_back(1);
        // n = 1 + 2 * p
        mpz_mul_2exp(tmp, p, 1);
        mpz_add_ui(n, tmp, 1);

        if (mpz_sizeinbase(n, 2) > target_bit_length)
        {
            // number already too large
            break;
        }
        if (mpz_composite_sieve(n) == COMPOSITE_FOR_SURE)
        {
            // number has a small factor
            break;
        }
        if ((mpz_get_ui(p) & 3) == 3)
        {
            // Euler-lagrange deterministic primality test 
	    // It is well known that if p==3 mod 4 is prime, 
	    // n = 2p+1 is also prime if and only if n divides 2^p-1.

            // 2^p mod n == 1
            mpz_powm(r, two, p, n);
            if (mpz_cmp_ui(r, 1) == 0)
            {
                // n is prime for sure
                next_prime(n);
            }
        }
        else
        {
            // Sophie-Germain deterministic primality test
	    // If p==1 mod 4 is prime, n=2p+1 is also prime if and only if n divides 2^p + 1.

            // 2^p mod n == n-1
            mpz_powm(r, two, p, n);
            if (mpz_cmp(r, tmp) == 0)
            {
                // n is prime for sure
                next_prime(n);
            }
        }
        break;
    }

    // Henri Lifschitz deterministic primality test
    static const unsigned mersenne_prime_exponent[] = {2, 3, 5, 7, 13, 17, 19, 31, 61, 127, 0};
    for (unsigned i = 0; mersenne_prime_exponent[i] != 0 && current_prime_count < max_prime_count; i++)
    {
        // verify large prime p > 2 * q
        if (mpz_cmp_ui(p, 2 * mersenne_prime_exponent[i]) > 0)
        {
            // n = 1 + 2 * p * q
            mpz_mul_ui(tmp, p, 2 * mersenne_prime_exponent[i]);
            mpz_add_ui(n, tmp, 1);

            if (mpz_sizeinbase(n, 2) > target_bit_length)
            {
                // number already too large
                break;
            }

	if (std::find(already_checked.begin(), already_checked.end(), mersenne_prime_exponent[i]) != already_checked.end())
	{
		// already tested
		continue;
	}
	already_checked.push_back(mersenne_prime_exponent[i]);

            if (mpz_composite_sieve(n) == COMPOSITE_FOR_SURE)
            {
                // number has a small factor
                continue;
            }
            // verify 2^n-1 mod n == 1
            mpz_powm(r, two, tmp, n);
            if (mpz_cmp_ui(r, 1) == 0)
            {
                // n is prime for sure
                next_prime(n);
            }
        }
    }

    // Pocklington deterministic primality test
    for (unsigned i = 0; i < SMALL_PRIMES_COUNT && current_prime_count < max_prime_count; i++)
    {
        // verify large prime p > 2 * q
        uint64_t q2 = 2 * small_primes[i];
        if (mpz_cmp_ui(p, q2) > 0)
        {
            // n = 1 + 2 * p * q
            mpz_mul_ui(tmp, p, q2);
            mpz_add_ui(n, tmp, 1);

            unsigned len_n = mpz_sizeinbase(n, 2);
            if (len_n > target_bit_length)
            {
                // number already too large
                break;
            }
	if (std::find(already_checked.begin(), already_checked.end(), small_primes[i]) != already_checked.end())
	{
		// already tested
		continue;
	}
	// already_checked.push_back(small_primes[i]);

            if (mpz_composite_sieve(n) == COMPOSITE_FOR_SURE)
            {
                // number has a small factor
                continue;
            }
            // verify 2^(n-1) mod n == 1
            mpz_powm(r, two, tmp, n);
            if (mpz_cmp_ui(r, 1) == 0)
            {
                // verify gcd(2^(2*q)-1, n) == 1
                if (q2 < len_n)
                {
                    // 1 << (2*q)-1
                    mpz_mul_2exp(tmp, r, q2);
                    mpz_sub_ui(r, tmp, 1);
                }
                else
                {
                    // (2^(2*q) % n + n - 1
                    mpz_set_ui(tmp, q2);
                    mpz_powm(r, two, tmp, n);
                    mpz_add(r, r, n);
                    mpz_sub_ui(r, r, 1);
                }
                mpz_gcd(tmp, r, n);
                if (mpz_cmp_ui(tmp, 1) == 0)
                {
                    // n is prime for sure
                    next_prime(n);
                }
            }
        }
    }

    mpz_clears(tmp, n, r, two, 0);
}

static void *worker_thread(void *arg)
{
    mpz_t p;
    mpz_init(p);

    while (1)
    {
        get_prime(&mt, p);
        if (mpz_sgn(p) == 0)
            break; // 0 : no more work to do

        worker(p);
    }

    mpz_clear(p);
    return 0;
}

int main(int argc, char **argv)
{
	mpz_t seed;
    unsigned thread_count = 1;
    mpz_init_set_ui(seed, 7);

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-b") && i < argc - 1)
        {
            target_bit_length = strtoul(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-s") && i < argc - 1)
        {
            mpz_expression_parse(seed, argv[++i]);
            continue;
        }
        if (!strcmp(argv[i], "-t") && i < argc - 1)
        {
            thread_count = strtoul(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-c") && i < argc-1)
        {
            max_prime_count = strtoull(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-verbose"))
        {
            verbose = true;
            continue;
        }
        if (!strcmp(argv[i], "-hex"))
        {
            hex_output = true;
            continue;
        }
        if (!strcmp(argv[i], "-dec"))
        {
            hex_output = false;
            continue;
        }
        printf("\n");
        printf("Generate proven primes\n");
        printf("-hex : hexadecimal output 0xf8ab01........\n");
        printf("-dec : decimal output 1234567.....\n");
        printf("-verbose : display all intermediate primes\n");
        printf("-s xxxxxx : prime seed, e.g. 607\n");
        printf("-b xxxxxx : bit length of output primes, e.g. 384\n");
        printf("-c xxxxxx : count of output primes, e.g. 1000\n");
        printf("-t xxxxxx : number of worker threads e.g. 3\n");
        printf("\n");
        exit(1);
    }

    // check parameters

    // make sure the seed is a prime  (BPSW test and 2 extra unnecessary rounds of Miller Rabin)
    if (!mpz_probab_prime_p(seed, 2))
    {
	    printf("Input seed must be prime (option -s)\n");
	    exit(1);
    }

    if (mpz_sizeinbase(seed, 2) > target_bit_length)
    {
	    printf("Input seed is too large (conflict on -s, -b options)\n");
	    exit(1);
    }

    if (thread_count < 1)
    {
        thread_count = 1;
    }

    // fill an array of small primes
    generate_small_primes();

    mt_initialize(&mt);

    pthread_t tid[thread_count];

    for (unsigned i = 0; i < thread_count; i++)
    {
        pthread_create(&tid[i], 0, worker_thread, 0);
    }

    // todo : use Morrisson test to find a first proven prime k * 2^n - 1
    // todo : use Proth test to find a first proven prime k * 2^n + 1
    // todo : use Lucas-Lehmer test to find a first proven prime 2^n-1

    // push the first prime
    next_prime(seed);
    mpz_clear(seed);

    while (current_prime_count < max_prime_count)
    {
        sleep(1);
	display_flush();
    }

    mpz_t zero;
    mpz_init_set_ui(zero, 0);
    for (unsigned i = 0; i < thread_count; i++)
    {
        post_prime(&mt, zero);
    }
    mpz_clear(zero);

    for (unsigned i = 0; i < thread_count; i++)
    {
        pthread_join(tid[i], 0);
    }

    mt_uninitialize(&mt);

    exit(0);
}
