#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "cubic_primality_alloc.h"
#include "cubic_primality_mt.h"

#define CACHE_LINE_SIZE 64
#define MT_QUEUE_SIZE 16
#define MIN_THREAD_COUNT 7

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
    unsigned start_queue_head, start_queue_tail;
    unsigned done_queue_head, done_queue_tail;
    unsigned start_queue[MT_QUEUE_SIZE];
    unsigned done_queue[MT_QUEUE_SIZE];
};

static void mod_multithread_init(mod_multithread_t *mt, mod_precompute_t *p)
{
    mt->p = p;
    mt->start_cond = PTHREAD_COND_INITIALIZER;
    mt->done_cond = PTHREAD_COND_INITIALIZER;
    mt->start_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->done_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt->start_queue_head = 0;
    mt->start_queue_tail = 0;
    mt->done_queue_head = 0;
    mt->done_queue_tail = 0;
    mpz_inits(mt->s, mt->t, mt->u, 0);
    mpz_inits(mt->s2, mt->t2, mt->u2, mt->st, mt->tu, mt->us, 0);
}

static void mod_multithread_end(mod_multithread_t *mt)
{
    mpz_clears(mt->s, mt->t, mt->u, 0);
    mpz_clears(mt->s2, mt->t2, mt->u2, mt->st, mt->tu, mt->us, 0);
}

static void mod_multithread_notify_ready(mod_multithread_t *mt, unsigned job)
{
    pthread_mutex_lock(&mt->start_mutex);
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
    pthread_mutex_unlock(&mt->done_mutex);
    return job;
}

// The little actions done by the worker threads
// First, get a job assignment
// then do the computations (all assignment involve multiplications or squaring)
// and return a status when done

// Ideally this would be the actions of a usual state machine if implemented as single thread.

#define FLAG_S      0x1
#define FLAG_T      0x2
#define FLAG_U      0x4
#define FLAG_S2     0x8
#define FLAG_T2    0x10
#define FLAG_U2    0x20
#define FLAG_ST    0x40
#define FLAG_TU    0x80
#define FLAG_US   0x100
#define FLAG_INIT 0x200
#define FLAG_E   0x8000
#define FLAG_END    0x0

static void *mod_worker(void *arg)
{
    unsigned job;
    struct mod_multithread_t *mt = (struct mod_multithread_t *)arg;

    // get a new job
    while ((job = mod_multithread_get_job(mt)) != FLAG_END)
    {
        switch (job)
        {
        case FLAG_S2:
            mpz_mul(mt->s2, mt->s, mt->s);
            mpz_mul_ui(mt->s2, mt->s2, mt->a);
            break;
        case FLAG_T2:
            mpz_mul(mt->t2, mt->t, mt->t);
            break;
        case FLAG_U2:
            mpz_mul(mt->u2, mt->u, mt->u);
            break;
        case FLAG_ST:
            mpz_mul(mt->st, mt->s, mt->t);
            mpz_mul_ui(mt->st, mt->st, 2 * mt->a);
            break;
        case FLAG_TU:
            mpz_mul(mt->tu, mt->t, mt->u);
            mpz_mul_2exp(mt->tu, mt->tu, 1);
            break;
        case FLAG_US:
            mpz_mul(mt->us, mt->u, mt->s);
            mpz_mul_2exp(mt->us, mt->us, 1);
            break;
        case FLAG_S | FLAG_E:
            mpz_add(mt->s, mt->s2, mt->st);
            mpz_add(mt->s, mt->s, mt->tu);
            mpz_mod_fast_reduce(mt->s, mt->p);
            break;
        case FLAG_T | FLAG_E:
            mpz_add(mt->t, mt->us, mt->t2);
            mpz_add(mt->t, mt->t, mt->s2);
            mpz_mul_ui(mt->t, mt->t, mt->a);
            mpz_add(mt->t, mt->t, mt->st);
            mpz_add(mt->t, mt->t, mt->u2);
            mpz_mod_fast_reduce(mt->t, mt->p);
            break;
        case FLAG_U | FLAG_E:
            mpz_add(mt->u, mt->us, mt->t2);
            mpz_add(mt->u, mt->u, mt->s2);
            mpz_mul_ui(mt->u, mt->u, mt->a);
            mpz_mod_fast_reduce(mt->u, mt->p);
            break;
        case FLAG_S:
            mpz_add(mt->s, mt->s2, mt->us);
            mpz_add(mt->s, mt->s, mt->t2);
            mpz_mod_fast_reduce(mt->s, mt->p);
            break;
        case FLAG_T:
            mpz_add(mt->t, mt->s2, mt->st);
            mpz_add(mt->t, mt->t, mt->tu);
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
        }

        // job is done, tell the main thread
        mod_multithread_notify_done(mt, job);
    }
    return 0;
}

// A state machine to handle the dependencies over different variables.
// Ideally this would be a traditional DFA built over a transition table.
// But considering the low number of transitions, it is simpler to write it as flat code.
//
// (This could be done in 9 transitions. But there are extra transitions to avoid
// delays and to speed up processing)
//
void inner_mutithread_exponentiate(mpz_t s, mpz_t t, mpz_t u, mpz_t e, uint64_t a, mod_precompute_t *p)
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
        result = mod_multithread_get_result(mt);
        if ((flags & (FLAG_U | FLAG_T)) == (FLAG_U | FLAG_T) && (result & FLAG_S))
        {
            flags |= FLAG_S;
            if (bit == 0)
                break;
            flags &= ~(FLAG_S2 | FLAG_ST | FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_S2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            exponent_bit = mpz_tstbit(e, --bit) ? FLAG_E : 0;
            continue;
        }
        if ((flags & (FLAG_U | FLAG_S)) == (FLAG_U | FLAG_S) && (result & FLAG_T))
        {
            flags |= FLAG_T;
            if (bit == 0)
                break;
            flags &= ~(FLAG_T2 | FLAG_ST | FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_ST);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_T2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            exponent_bit = mpz_tstbit(e, --bit) ? FLAG_E : 0;
            continue;
        }
        if ((flags & (FLAG_T | FLAG_S)) == (FLAG_S | FLAG_T) && (result & FLAG_U))
        {
            flags |= FLAG_U;
            if (bit == 0)
                break;
            flags &= ~(FLAG_U2 | FLAG_US | FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_US);
            mod_multithread_notify_ready(mt, FLAG_TU);
            mod_multithread_notify_ready(mt, FLAG_U2);
            flags &= ~(FLAG_S | FLAG_T | FLAG_U);
            exponent_bit = mpz_tstbit(e, --bit) ? FLAG_E : 0;
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
            flags &= ~(FLAG_ST | FLAG_U2);
            mod_multithread_notify_ready(mt, FLAG_ST);
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
	mask = FLAG_TU | FLAG_ST | FLAG_T2 | FLAG_US;
        if ((flags & mask) == mask && (result & FLAG_S2))
        {
            flags &= ~(FLAG_S | FLAG_T);
            mod_multithread_notify_ready(mt, FLAG_S | exponent_bit);
            mod_multithread_notify_ready(mt, FLAG_T | exponent_bit);
            flags |= FLAG_S2;
            continue;
        }
	mask = FLAG_S2 | FLAG_TU | FLAG_U2;
        if ((flags & mask) == mask && (result & FLAG_ST))
        {
            flags &= ~(FLAG_T | FLAG_U);
            mod_multithread_notify_ready(mt, FLAG_T | exponent_bit);
            mod_multithread_notify_ready(mt, FLAG_U | exponent_bit);
            flags |= FLAG_S2;
            continue;
        }
        if ((flags & (FLAG_S2 | FLAG_US)) == (FLAG_S2 | FLAG_US) && (result & FLAG_T2))
	{
            flags &= ~FLAG_S;
            mod_multithread_notify_ready(mt, FLAG_S | exponent_bit);
            flags |= FLAG_T2;
            continue;
	}
        if ((flags & (FLAG_T2 | FLAG_US)) == (FLAG_T2 | FLAG_US) && (result & FLAG_S2))
	{
            flags &= ~FLAG_S;
            mod_multithread_notify_ready(mt, FLAG_S | exponent_bit);
            flags |= FLAG_S2;
            continue;
	}
        if ((flags & (FLAG_S2 | FLAG_T2)) == (FLAG_S2 | FLAG_T2) && (result & FLAG_US))
	{
            flags &= ~FLAG_S;
            mod_multithread_notify_ready(mt, FLAG_S | exponent_bit);
            flags |= FLAG_US;
            continue;
	}
        if ((flags & (FLAG_S2 | FLAG_ST)) == (FLAG_S2 | FLAG_ST) && (result & FLAG_TU))
	{
            flags &= ~FLAG_T;
            mod_multithread_notify_ready(mt, FLAG_T | exponent_bit);
            flags |= FLAG_TU;
            continue;
	}
        if ((flags & (FLAG_S2 | FLAG_TU)) == (FLAG_S2 | FLAG_TU) && (result & FLAG_ST))
	{
            flags &= ~FLAG_T;
            mod_multithread_notify_ready(mt, FLAG_T | exponent_bit);
            flags |= FLAG_ST;
            continue;
	}
        if ((flags & (FLAG_ST | FLAG_TU)) == (FLAG_ST | FLAG_TU) && (result & FLAG_S2))
	{
            flags &= ~FLAG_T;
            mod_multithread_notify_ready(mt, FLAG_T | exponent_bit);
            flags |= FLAG_S2;
            continue;
	}
        if ((flags & (FLAG_ST)) == (FLAG_ST) && (result & FLAG_U2))
        {
            flags &= ~FLAG_U;
            mod_multithread_notify_ready(mt, FLAG_U | exponent_bit);
            flags |= FLAG_ST;
            continue;
        }
        if ((flags & (FLAG_U2)) == (FLAG_U2) && (result & FLAG_ST))
        {
            flags &= ~FLAG_U;
            mod_multithread_notify_ready(mt, FLAG_U | exponent_bit);
            flags |= FLAG_U2;
            continue;
        }
        flags |= result;
    }

    /*
    T[0][FLAG_S] = 1
    T[0][FLAG_T] = 2
    T[0][FLAG_U] = 3
    T[FLAG_S][FLAG_T] = 2 4
    T[FLAG_S][FLAG_U] = 3 6
    T[FLAG_T][FLAG_S] = 1 4
    T[FLAG_T][FLAG_U] = 3 5
    T[FLAG_U][FLAG_S] = 1 6
    T[FLAG_U][FLAG_T] = 2 5
    T[FLAG_U|FLAG_T][FLAG_S] = 1 4 6
    T[FLAG_S|FLAG_U][FLAG_T] = 2 4 5
    T[FLAG_S|FLAG_T][FLAG_U] = 3 6 5
    T[FLAG_S2|FLAG_US][FLAG_T2] = 10
    T[FLAG_T2|FLAG_US][FLAG_S2] = 10
    T[FLAG_T2|FLAG_S2][FLAG_US] = 10
    T[FLAG_S2|FLAG_ST][FLAG_TU] = 11
    T[FLAG_S2|FLAG_TU][FLAG_ST] = 11
    T[FLAG_TU|FLAG_ST][FLAG_S2] = 11
    T[FLAG_U2][FLAG_ST] = 12
    T[FLAG_ST][FLAG_U2] = 12
    T[FLAG_S2|FLAG_TU|FLAG_U2][FLAG_ST] = 11 12
    T[FLAG_TU|FLAG_ST|FLAG_T2|FLAG_US][FLAG_S2] = 10 11
    */

    // pointer copies
    mpz_set(s, mt->s);
    mpz_set(t, mt->t);
    mpz_set(u, mt->u);

    // gracefully stop the worker threads
    for (unsigned i = 0; i < MIN_THREAD_COUNT; i++)
    {
        mod_multithread_notify_ready(mt, FLAG_END);
    }

    // wait for the worker threads termination
    for (unsigned i = 0; i < MIN_THREAD_COUNT; i++)
    {
        pthread_join(tids[i], 0);
    }

    mod_multithread_end(mt);
}
