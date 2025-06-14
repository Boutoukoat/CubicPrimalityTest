#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "cubic_primality_alloc.h"
#include "cubic_primality_mt.h"

#define CACHE_LINE_SIZE 64
#define MT_QUEUE_SIZE 16
#define MIN_THREAD_COUNT 7

#define FLAG_S 0x0001
#define FLAG_T 0x0002
#define FLAG_U 0x0004
#define FLAG_S2 0x0010
#define FLAG_T2 0x0020
#define FLAG_U2 0x0040
#define FLAG_ST 0x0100
#define FLAG_TU 0x0200
#define FLAG_US 0x0400
#define FLAG_INIT 0x1000
#define FLAG_E 0x8000
#define FLAG_END 0x0000

struct mod_multithread_t
{
    mpz_t s __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t s2 __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t t __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t t2 __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t u __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t u2 __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t st __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t tu __attribute__((aligned(CACHE_LINE_SIZE)));
    mpz_t us __attribute__((aligned(CACHE_LINE_SIZE)));

    mod_precompute_t *p;
    uint64_t a;

    pthread_cond_t start_cond, done_cond;
    pthread_mutex_t start_mutex, done_mutex;
    volatile unsigned start_queue_head, start_queue_tail;
    volatile unsigned done_queue_head, done_queue_tail;
    unsigned start_queue[MT_QUEUE_SIZE];
    unsigned done_queue[MT_QUEUE_SIZE];
};

static void mod_multithread_init(mod_multithread_t *mt, mod_precompute_t *p)
{
    unsigned min_size = (p->n + 256) * 2;
    mt->p = p;
    mt->start_cond = PTHREAD_COND_INITIALIZER;
    mt->done_cond = PTHREAD_COND_INITIALIZER;
    mt->start_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->done_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->start_queue_head = 0;
    mt->start_queue_tail = 0;
    mt->done_queue_head = 0;
    mt->done_queue_tail = 0;
    mpz_init2(mt->s, min_size);
    mpz_init2(mt->t, min_size);
    mpz_init2(mt->u, min_size);
    mpz_init2(mt->s2, min_size);
    mpz_init2(mt->t2, min_size);
    mpz_init2(mt->u2, min_size);
    mpz_init2(mt->st, min_size);
    mpz_init2(mt->tu, min_size);
    mpz_init2(mt->us, min_size);
}

static void mod_multithread_end(mod_multithread_t *mt)
{
    mpz_clears(mt->s, mt->t, mt->u, 0);
    mpz_clears(mt->s2, mt->t2, mt->u2, mt->st, mt->tu, mt->us, 0);
}

static void mod_multithread_notify_ready(mod_multithread_t *mt, unsigned job)
{
    pthread_mutex_lock(&mt->start_mutex);
    // printf("Notify Ready %4.4x\n", job);
    mt->start_queue[mt->start_queue_head++ & (MT_QUEUE_SIZE - 1)] = job;
    pthread_mutex_unlock(&mt->start_mutex);
    pthread_cond_signal(&mt->start_cond);
}

static void mod_multithread_notify_done(mod_multithread_t *mt, unsigned job)
{
    pthread_mutex_lock(&mt->done_mutex);
    mt->done_queue[mt->done_queue_head++ & (MT_QUEUE_SIZE - 1)] = job;
    pthread_mutex_unlock(&mt->done_mutex);
    pthread_cond_signal(&mt->done_cond);
}

static unsigned mod_multithread_get_job(mod_multithread_t *mt)
{
    unsigned job = 0;
    pthread_mutex_lock(&mt->start_mutex);
    while (mt->start_queue_tail == mt->start_queue_head)
    {
        pthread_cond_wait(&mt->start_cond, &mt->start_mutex);
    }
    job = mt->start_queue[mt->start_queue_tail++ & (MT_QUEUE_SIZE - 1)];
    pthread_mutex_unlock(&mt->start_mutex);
    return job;
}

static unsigned mod_multithread_get_result(mod_multithread_t *mt)
{
    unsigned job = 0;
    pthread_mutex_lock(&mt->done_mutex);
    while (mt->done_queue_tail == mt->done_queue_head)
    {
        pthread_cond_wait(&mt->done_cond, &mt->done_mutex);
    }
    job = mt->done_queue[mt->done_queue_tail++ & (MT_QUEUE_SIZE - 1)];
    //    printf("Get result %4.4x\n", job);
    pthread_mutex_unlock(&mt->done_mutex);
    return job;
}

// The little actions done by the worker threads
// First, get a job assignment
// then do the computations (all assignment involve multiplications or squaring)
// and return a status when done

// Ideally this would be the actions of a usual state machine if implemented as single thread.

static void *mod_worker(void *arg)
{
    unsigned job;
    struct mod_multithread_t *mt = (struct mod_multithread_t *)arg;
    mpz_t tmp;
    unsigned min_size = (mt->p->n + 256) * 2;
    mpz_init2(tmp, min_size);

    // get a new job
    while ((job = mod_multithread_get_job(mt)) != FLAG_END)
    {
        switch (job)
        {
        case FLAG_S2:
            mpz_mul(tmp, mt->s, mt->s);
            mpz_mul_ui(mt->s2, tmp, mt->a);
            break;
        case FLAG_T2:
            mpz_mul(mt->t2, mt->t, mt->t);
            break;
        case FLAG_U2:
            mpz_mul(mt->u2, mt->u, mt->u);
            break;
        case FLAG_ST:
            mpz_mul(tmp, mt->s, mt->t);
            mpz_mul_ui(mt->st, tmp, 2 * mt->a);
            break;
        case FLAG_TU:
            mpz_mul(tmp, mt->t, mt->u);
            mpz_mul_2exp(mt->tu, tmp, 1);
            break;
        case FLAG_US:
            mpz_mul(tmp, mt->u, mt->s);
            mpz_mul_2exp(mt->us, tmp, 1);
            break;
        case FLAG_S | FLAG_E:
            mpz_add(tmp, mt->s2, mt->st);
            mpz_add(mt->s, tmp, mt->tu);
            mpz_mod_fast_reduce(mt->s, mt->p);
            break;
        case FLAG_T | FLAG_E:
            mpz_add(tmp, mt->us, mt->t2);
            mpz_add(tmp, tmp, mt->s2);
            mpz_mul_ui(mt->t, tmp, mt->a);
            mpz_add(tmp, mt->t, mt->st);
            mpz_add(mt->t, tmp, mt->u2);
            mpz_mod_fast_reduce(mt->t, mt->p);
            break;
        case FLAG_U | FLAG_E:
            mpz_add(tmp, mt->us, mt->t2);
            mpz_add(tmp, tmp, mt->s2);
            mpz_mul_ui(mt->u, tmp, mt->a);
            mpz_mod_fast_reduce(mt->u, mt->p);
            break;
        case FLAG_S:
            mpz_add(tmp, mt->s2, mt->us);
            mpz_add(mt->s, tmp, mt->t2);
            mpz_mod_fast_reduce(mt->s, mt->p);
            break;
        case FLAG_T:
            mpz_add(tmp, mt->s2, mt->st);
            mpz_add(mt->t, tmp, mt->tu);
            mpz_mod_fast_reduce(mt->t, mt->p);
            break;
        case FLAG_U:
            mpz_add(mt->u, mt->u2, mt->st);
            mpz_mod_fast_reduce(mt->u, mt->p);
            break;
        case FLAG_INIT:
            // little tweak to initiate the state machine (s, t, u are already computed)
            mod_multithread_notify_done(mt, FLAG_S);
            mod_multithread_notify_done(mt, FLAG_T);
            job = FLAG_U;
            break;
        default:
            // cannot happen
            assert(0);
            break;
        }

        // job is done, tell the main thread
        mod_multithread_notify_done(mt, job);
    }

    mpz_clear(tmp);

    return 0;
}

// A state machine to handle the dependencies over different variables.
// Ideally this would be a traditional DFA built over a transition table.
// But considering the low number of transitions, it is simpler to write it as flat code.
//
// (This could be done in 9 transitions. But there are extra transitions to avoid
// delays and to speed up processing)
//
void mpz_inner_multithread_exponentiate(mpz_t s, mpz_t t, mpz_t u, mpz_t e, uint64_t a, mod_precompute_t *p)
{
    struct mod_multithread_t *mt =
        (struct mod_multithread_t *)cubic_allocate_function(sizeof(struct mod_multithread_t));
    mod_multithread_init(mt, p);
    pthread_t tids[MIN_THREAD_COUNT];
    unsigned bit = mpz_sizeinbase(e, 2) - 1;
    unsigned exponent_bit = 0;
    unsigned flags = 0;
    unsigned result = 0;
    unsigned mask = 0;

    // pointer copies
    mpz_set(mt->s, s);
    mpz_set(mt->t, t);
    mpz_set(mt->u, u);
    mt->a = a;

    // kick-off the worker threads
    for (unsigned i = 0; i < MIN_THREAD_COUNT; i++)
    {
        pthread_create(&tids[i], 0, mod_worker, mt);
    }

    // kick-off the state machine to a hard-coded initialization transition
    mod_multithread_notify_ready(mt, FLAG_INIT);

    while (true)
    {
        mask = FLAG_US | FLAG_S2 | FLAG_T2 | FLAG_ST | FLAG_U2 | FLAG_TU;
        if ((flags & mask) == mask)
        {
                flags &= ~mask;
        }

       	result = mod_multithread_get_result(mt);
        // printf("Flags %4.4x Result %4.4x\n", flags, result);

        if ((flags & (FLAG_U | FLAG_T)) == (FLAG_U | FLAG_T) && (result & FLAG_S))
        {
            if (bit == 0)
                break;
            flags &= ~(FLAG_ST | FLAG_US | FLAG_S2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_S2);
            exponent_bit = mpz_tstbit(e, --bit);
            continue;
        }
        if ((flags & (FLAG_U | FLAG_S)) == (FLAG_U | FLAG_S) && (result & FLAG_T))
        {
            if (bit == 0)
                break;
            flags &= ~(FLAG_ST | FLAG_TU | FLAG_T2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_T2);
            exponent_bit = mpz_tstbit(e, --bit);
            continue;
        }
        if ((flags & (FLAG_T | FLAG_S)) == (FLAG_S | FLAG_T) && (result & FLAG_U))
        {
            if (bit == 0)
                break;
            flags &= ~(FLAG_US | FLAG_TU | FLAG_U2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_U2);
            exponent_bit = mpz_tstbit(e, --bit);
            continue;
        }
        if ((flags & (FLAG_T)) == (FLAG_T) && (result & FLAG_S))
        {
            flags |= FLAG_S;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_ST | FLAG_S2);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_S2);
            continue;
        }
        if ((flags & (FLAG_U)) == (FLAG_U) && (result & FLAG_S))
        {
            flags |= FLAG_S;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_US | FLAG_S2);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_S2);
            continue;
        }
        if ((flags & (FLAG_S)) == (FLAG_S) && (result & FLAG_T))
        {
            flags |= FLAG_T;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_ST | FLAG_T2);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_T2);
            continue;
        }
        if ((flags & (FLAG_U)) == (FLAG_U) && (result & FLAG_T))
        {
            flags |= FLAG_T;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_TU | FLAG_T2);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_T2);
            continue;
        }
        if ((flags & (FLAG_S)) == (FLAG_S) && (result & FLAG_U))
        {
            flags |= FLAG_U;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_US | FLAG_U2);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_U2);
            continue;
        }
        if ((flags & (FLAG_T)) == (FLAG_T) && (result & FLAG_U))
        {
            flags |= FLAG_U;
            if (bit == 0)
                continue;
            flags &= ~(FLAG_TU | FLAG_U2);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_U2);
            continue;
        }
        if (result & FLAG_S)
        {
            flags |= FLAG_S;
            if (bit == 0)
                continue;
            flags &= ~FLAG_S2;
            mod_multithread_notify_ready(mt, FLAG_S2);
            continue;
        }
        if (result & FLAG_T)
        {
            flags |= FLAG_T;
            if (bit == 0)
                continue;
            flags &= ~FLAG_T2;
            mod_multithread_notify_ready(mt, FLAG_T2);
            continue;
        }
        if (result & FLAG_U)
        {
            flags |= FLAG_U;
            if (bit == 0)
                continue;
            flags &= ~FLAG_U2;
            mod_multithread_notify_ready(mt, FLAG_U2);
            continue;
        }

        if (exponent_bit)
        {
            // exponent_bit == 1
            mask = FLAG_US | FLAG_S2 | FLAG_T2;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_U;
                mod_multithread_notify_ready(mt, FLAG_U | FLAG_E);
            }

            mask = FLAG_TU | FLAG_S2 | FLAG_ST;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_S;
                mod_multithread_notify_ready(mt, FLAG_S | FLAG_E);
            }

            mask = FLAG_US | FLAG_S2 | FLAG_T2 | FLAG_ST | FLAG_U2;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_T;
                mod_multithread_notify_ready(mt, FLAG_T | FLAG_E);
            }
        }
        else
        {
            // exponent_bit == 0
            mask = FLAG_S2 | FLAG_US | FLAG_T2;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_S;
                mod_multithread_notify_ready(mt, FLAG_S);
            }

            mask = FLAG_S2 | FLAG_TU | FLAG_ST;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_T;
                mod_multithread_notify_ready(mt, FLAG_T);
            }

            mask = FLAG_U2 | FLAG_ST;
            if ((flags & mask) != mask && ((flags | result) & mask) == mask)
            {
                flags &= ~FLAG_U;
                mod_multithread_notify_ready(mt, FLAG_U);
            }
        }

        // default case
        flags |= result;
    }

    // gracefully stop the worker loops
    for (unsigned i = 0; i < MIN_THREAD_COUNT; i++)
    {
        mod_multithread_notify_ready(mt, FLAG_END);
    }

    // wait for the worker thread terminations
    for (unsigned i = 0; i < MIN_THREAD_COUNT; i++)
    {
        pthread_join(tids[i], 0);
    }

    // pointer copies
    mpz_set(s, mt->s);
    mpz_set(t, mt->t);
    mpz_set(u, mt->u);

    mod_multithread_end(mt);
}
