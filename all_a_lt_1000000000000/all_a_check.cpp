// verify details of the cubic test
//
// {gettime();forstep(n=25,1000000000000,2,
// if(n%10000==1,print([n,gettime()]));
// if(n%3!=0&&!ispseudoprime(n),p=factor(n)[1,1];Q=n/p;g=gcd(p^3-1,Q^3-1);R=(n-1)%g;
//  if(R>2,A=(n^2+n+1)%g;
//   if(A>3,
//    for(k=1,(n+1)/2,
//     if(k%3!=2,a=(7+k*(k-1))%n;
//      if(Mod(a,n)^R==1,B=Mod(Mod(x,n),x^3-a*x-a)^R;
//       if(B^2+B+1==-x^2+x+a,print([n,p,Q,g,R,A]);break(7)))))))));}
//
// -----------------------------------------------------------------------

#include <assert.h>

#include "math_utils.cpp"

#include "v_table.cpp"
#define V_COUNT (sizeof(V) / sizeof(V[0]))

#include <pthread.h>
#include <semaphore.h>

// compute Mod(Mod(x, n), x^3 -a*x -a)^e
// assume n is less than 60 or 61 bits
// assume a < n
// if e odd, assume s == 0 and u == 0 and t == 1
static void cubic_exponentiate(uint64_t &s, uint64_t &t, uint64_t &u, uint64_t e, uint64_t n, uint64_t a)
{
    int bit = uint64_log_2(e);
    uint64_t tmp;
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;
    while (bit--)
    {
        // start Square
        tmp = square_mod(s, n);
        s2 = (uint128_t)tmp * a;
        t2 = (uint128_t)t * t;
        u2 = (uint128_t)u * u;
        tmp = mul_mod(s, t, n);
        st = (uint128_t)tmp * a;
        tu = (uint128_t)t * u;
        us = (uint128_t)u * s;
        st <<= 1;
        tu <<= 1;
        us <<= 1;
        if (e & (1ull << bit))
        {
            // finish Square and multiply by x
            us = s2 + us + t2;
            tmp = uint128_long_mod(us, n);
            uu = (uint128_t)tmp * a;
            ss = s2 + st + tu;
            tt = uu + u2 + st;
        }
        else
        {
            // finish Square
            ss = s2 + us + t2;
            tt = s2 + st + tu;
            uu = st + u2;
        }
        s = uint128_long_mod(ss, n);
        t = uint128_long_mod(tt, n);
        u = uint128_long_mod(uu, n);
    }

    // make sure the result is reduced
    s = uint64_long_mod(0, s, n);
    t = uint64_long_mod(0, t, n);
    u = uint64_long_mod(0, u, n);
}

// compute Mod(Mod(x, n), x^3 -a*x -a)^e with Barrett reduction
static void barrett_cubic_exponentiate(uint64_t &s, uint64_t &t, uint64_t &u, uint64_t e, const barrett_t &bt,
                                       uint64_t a)
{
    int bit = uint64_log_2(e);
    uint64_t tmp;
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;
    while (bit--)
    {
        // start Square
        tmp = barrett_mul_mod(s, s, bt);
        s2 = (uint128_t)tmp * a;
        t2 = (uint128_t)t * t;
        u2 = (uint128_t)u * u;
        tmp = barrett_mul_mod(s, t, bt);
        st = (uint128_t)tmp * a;
        tu = (uint128_t)t * u;
        us = (uint128_t)u * s;
        st <<= 1;
        tu <<= 1;
        us <<= 1;
        if (e & (1ull << bit))
        {
            // finish Square and multiply by x
            us = s2 + us + t2;
            tmp = barrett_long_mod(us, bt);
            uu = (uint128_t)tmp * a;
            ss = s2 + st + tu;
            tt = uu + u2 + st;
        }
        else
        {
            // finish Square
            ss = s2 + us + t2;
            tt = s2 + st + tu;
            uu = st + u2;
        }
        // make sure the result s,t,u is <= 2 * m
        s = barrett_long_mod(ss, bt);
        t = barrett_long_mod(tt, bt);
        u = barrett_long_mod(uu, bt);
    }

    // make sure the result s,t,u is < m
    s = uint64_long_mod(0, s, bt.m);
    t = uint64_long_mod(0, t, bt.m);
    u = uint64_long_mod(0, u, bt.m);
}

// compute Mod(Mod(x, n), x^3 -a*x -a)^2 with Barrett reduction
static void barrett_cubic_square(uint64_t &s, uint64_t &t, uint64_t &u, const barrett_t &bt, uint64_t a)
{
    uint64_t tmp;
    uint128_t s2, t2, u2, st, tu, us, uu, ss, tt;

    // start Square
    tmp = barrett_mul_mod(s, s, bt);
    s2 = (uint128_t)tmp * a;
    t2 = (uint128_t)t * t;
    u2 = (uint128_t)u * u;
    tmp = barrett_mul_mod(s, t, bt);
    st = (uint128_t)tmp * a;
    tu = (uint128_t)t * u;
    us = (uint128_t)u * s;
    st <<= 1;
    tu <<= 1;
    us <<= 1;
    // finish Square
    ss = s2 + us + t2;
    tt = s2 + st + tu;
    uu = st + u2;
    // make sure the result s,t,u is < m
    s = uint128_long_mod(ss, bt.m);
    t = uint128_long_mod(tt, bt.m);
    u = uint128_long_mod(uu, bt.m);
}

// read file, cautious about corruped or truncated files
bool read_file(char *fn, uint64_t *n, uint64_t *m, uint64_t *s)
{
    char line[100];
    bool n_found = false;
    bool m_found = false;
    bool s_found = false;
    bool valid = true;
    uint64_t n_temp = -1, m_temp = -1, s_temp = -1;
    FILE *f = fopen(fn, "rt");
    if (f)
    {
        memset(line, 0, 100);
        while (fgets(line, 100, f))
        {
            // trim the white-spaces at end of line (noisy stuff from hand-written files)
            int l = strlen(line);
            while (l > 0 && isspace(line[l - 1]))
            {
                l--;
            }
            line[l] = 0;

            if (*line == 0 || *line == '#')
            {
                // empty line, or comment, ignore
                continue;
            }
            if (!memcmp(line, "n=", 2))
            {
                n_temp = strtoull(&line[2], 0, 0);
                n_found = true;
                continue;
            }
            if (!memcmp(line, "m=", 2))
            {
                m_temp = strtoull(&line[2], 0, 0);
                m_found = true;
                continue;
            }
            if (!memcmp(line, "s=", 2))
            {
                s_temp = strtoull(&line[2], 0, 0);
                s_found = true;
                continue;
            }

            // not our file, give up
            valid = false;
            break;
        }
        fclose(f);
    }
    else
    {
        valid = false;
    }

    if (valid && n_found == true && m_found == true && s_found == true)
    {
        // file is complete, without missing field
        *n = n_temp;
        *m = m_temp;
        *s = s_temp;
        return true;
    }
    return false;
}

// write file, cautious about worst cases in the life of a computer
bool write_file(char *fn, char *bak, uint64_t n, uint64_t m, uint64_t s)
{
    uint64_t n_temp = -1, m_temp = -1, s_temp = -1;
    // first verify the current file is a valid one
    if (read_file(fn, &n_temp, &m_temp, &s_temp))
    {
        // Create a backup file from the file just validated
        // Unix : atomic file rename.
        // In case of crash, operation is done, or not done at all, and worst case is to have duplicate files.
        rename(fn, bak);
    }

    // create or overwrite the file with new values,
    // and in case of crash and new file is incomplete or corrupted, there is a backup just done, pfew.
    FILE *f = fopen(fn, "wt");
    if (f)
    {
        // write current time
        char buf[26];
        struct tm time_m;
        time_t time_p;
        time(&time_p);
        localtime_r(&time_p, &time_m);
        asctime_r(&time_m, buf);
        fprintf(f, "# %s\n", buf);

        // write restart parameters
        fprintf(f, "n=%ld\n", n);
        fprintf(f, "m=%ld\n", m);
        fprintf(f, "s=%ld\n", s);
        fclose(f);

        // never paranoid enough against full disks, need to read back the file.
        // This also flushes the file in the case of remote file (read-after-write security on NFS ....)
        n_temp = -1;
        m_temp = -1;
        s_temp = -1;
        if (read_file(fn, &n_temp, &m_temp, &s_temp))
        {
            // return true if the file is correctly written
            return (n_temp == n && m_temp == m && s_temp == s);
        }
    }
    return false;
}

// -----------------------------------------------------------------
//
// useful debug counters
//
// -----------------------------------------------------------------

static uint64_t modexp_count = 0;
static uint64_t exponentiate_count = 0;
static uint64_t skip_squarefree_count = 0;
static uint64_t worked_squarefree_count = 0;
static uint64_t skip_semiprime_count = 0;
static uint64_t worked_semiprime_count = 0;
static uint64_t skip_composite_count = 0;
static uint64_t worked_composite_count = 0;

// -----------------------------------------------------------------
//
// utility function
//
// -----------------------------------------------------------------

void verify_all_a(uint64_t n, uint64_t R_fermat, uint64_t R_cubic)
{

    struct barrett_t bn;
    barrett_precompute(&bn, n);
    uint64_t a = 7;
    for (uint64_t d = 0; d < n; d += 2)
    {
        // a = 7 + k * (k-1) unrolled as  a = a + 2k
        a += d;
        a -= (a >= n) ? n : 0;

        // verify test Mod(a,n)^R_fermat == 1
        modexp_count += 1;
        if (barrett_pow_mod(a, R_fermat, bn) == 1)
        {
            exponentiate_count += 1;
            // run cubic test
            // B = Mod(x, n)
            uint64_t bs = 0;
            uint64_t bt = 1;
            uint64_t bu = 0;
            // B = Mod(B, x^3 -ax -a)^R_cubic
            barrett_cubic_exponentiate(bs, bt, bu, R_cubic, bn, a);
            // B2 = B
            uint64_t bs2 = bs;
            uint64_t bt2 = bt;
            uint64_t bu2 = bu;
            // B2 = Mod(B2, x^3 -ax -a)^2
            barrett_cubic_square(bs2, bt2, bu2, bn, a);
            // check B2+B+1 is NOT -x^2 + x + 1
            bs2 = uint64_add_mod(bs2, bs, n);
            if (bs2 == n - 1)
            {
                bt2 = uint64_add_mod(bt2, bt, n);
                if (bt2 == 1)
                {
                    bu2 = uint64_add_mod(bu2, bu + 1, n);
                    if (bu2 == a)
                    {
                        // Aoutch, n = p*Q could be prime or pseudoprime ??!!??!!
                        char buff[256];
                        char *ptr = buff;
                        *ptr++ = '[';
                        ptr += uint128_sprint(ptr, n);
                        *ptr++ = ',';
                        ptr += uint128_sprint(ptr, a);
                        *ptr++ = ',';
                        ptr += uint128_sprint(ptr, R_fermat);
                        *ptr++ = ',';
                        ptr += uint128_sprint(ptr, R_cubic);
                        *ptr++ = ']';
                        *ptr = 0;
                        printf("%s\n", buff);
                        fflush(stdout);

                        // got 1 counter-example, it is time to stop
                        exit(1);
                    }
                    else
                    {
                        // polynomial test :
                        // n is composite for sure
                    }
                }
                else
                {
                    // polynomial test :
                    // n is composite for sure
                }
            }
            else
            {
                // polynomial test :
                // n is composite for sure
            }
        }
        else
        {
            // modexp test :
            // n is composite for sure
        }
    }
}

// -----------------------------------------------------------------
//
// multithread wrapper
//
// -----------------------------------------------------------------

#define MT_THREAD_COUNT 4 // a power of 2
#define MT_QUEUE_SIZE (2 * MT_THREAD_COUNT)

struct ring_entry_t
{
    uint64_t __attribute__((aligned(64))) n;
    uint64_t R_fermat;
    uint64_t R_cubic;
};

struct mod_multithread_t
{
    pthread_cond_t __attribute__((aligned(64))) start_cond, done_cond;
    pthread_mutex_t __attribute__((aligned(64))) start_mutex, done_mutex;
    volatile unsigned long __attribute__((aligned(64))) start_queue_head, start_queue_tail;
    volatile unsigned long __attribute__((aligned(64))) done_queue_head, done_queue_tail;
    sem_t __attribute__((aligned(64))) start_sem, done_sem;
    struct ring_entry_t start_queue[MT_QUEUE_SIZE];
    struct ring_entry_t done_queue[MT_QUEUE_SIZE];

    uint64_t display;
    time_t d0;
};

static void mod_multithread_notify_ready(mod_multithread_t *mt, const struct ring_entry_t &job)
{
    sem_wait(&mt->start_sem);
    pthread_mutex_lock(&mt->start_mutex);
    // printf("Notify Ready %ld\n", job.n);
    mt->start_queue[mt->start_queue_head++ & (MT_QUEUE_SIZE - 1)] = job;
    pthread_mutex_unlock(&mt->start_mutex);
    pthread_cond_signal(&mt->start_cond);
}

static void mod_multithread_notify_done(mod_multithread_t *mt, const struct ring_entry_t &job)
{
    sem_wait(&mt->done_sem);
    pthread_mutex_lock(&mt->done_mutex);
    // printf("Job done %ld\n", job.n);
    mt->done_queue[mt->done_queue_head++ & (MT_QUEUE_SIZE - 1)] = job;
    pthread_mutex_unlock(&mt->done_mutex);
    pthread_cond_signal(&mt->done_cond);
    sem_post(&mt->done_sem);
}

static struct ring_entry_t mod_multithread_get_job(mod_multithread_t *mt)
{
    struct ring_entry_t job;
    pthread_mutex_lock(&mt->start_mutex);
    while (mt->start_queue_tail == mt->start_queue_head)
    {
        pthread_cond_wait(&mt->start_cond, &mt->start_mutex);
    }
    job = mt->start_queue[mt->start_queue_tail++ & (MT_QUEUE_SIZE - 1)];
    // printf("Get job %ld\n", job.n);
    pthread_mutex_unlock(&mt->start_mutex);
    sem_post(&mt->start_sem);
    return job;
}

static struct ring_entry_t mod_multithread_get_result(mod_multithread_t *mt)
{
    struct ring_entry_t job;
    pthread_mutex_lock(&mt->done_mutex);
    while (mt->done_queue_tail == mt->done_queue_head)
    {
        pthread_cond_wait(&mt->done_cond, &mt->done_mutex);
    }
    job = mt->done_queue[mt->done_queue_tail++ & (MT_QUEUE_SIZE - 1)];
    // printf("Get result %ld\n", job.n);
    pthread_mutex_unlock(&mt->done_mutex);
    sem_post(&mt->done_sem);
    return job;
}

static struct mod_multithread_t mt;
static pthread_t mt_tids[MT_THREAD_COUNT];

void *mt_worker(void *)
{
    struct ring_entry_t job;

    while (1)
    {
        // get a new job
        job = mod_multithread_get_job(&mt);
        if (job.n == 0)
            break;
        verify_all_a(job.n, job.R_fermat, job.R_cubic);
        // job is done, tell the main thread
        mod_multithread_notify_done(&mt, job);
    }
    return 0;
}

void drain_done_queue(mod_multithread_t *mt)
{
    // best effort, non blocking attempt to drain whatever is in the queue
    while (mt->done_queue_tail != mt->done_queue_head)
    {
        struct ring_entry_t job = mod_multithread_get_result(mt);

        // multiple thread progress is approximative and out of order
        if (++mt->display > 1000)
        {
            time_t d1 = time(NULL);
            printf("completed n = %ld, %ld\n", job.n, d1 - mt->d0);
            mt->display = 0;
            mt->d0 = d1;
        }
    }
}

void mt_initialize(void)
{
    mt.start_cond = PTHREAD_COND_INITIALIZER;
    mt.done_cond = PTHREAD_COND_INITIALIZER;
    mt.start_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt.done_mutex = PTHREAD_MUTEX_INITIALIZER;
    mt.start_queue_head = 0;
    mt.start_queue_tail = 0;
    mt.done_queue_head = 0;
    mt.done_queue_tail = 0;
    sem_init(&mt.start_sem, 0, MT_QUEUE_SIZE);
    sem_init(&mt.done_sem, 0, MT_QUEUE_SIZE);

    mt.display = 0;
    mt.d0 = time(NULL);

    // kick-off the worker threads
    for (unsigned i = 0; i < MT_THREAD_COUNT; i++)
    {
        pthread_create(&mt_tids[i], 0, mt_worker, &mt);
    }
}

void mt_debug_display(void)
{
    printf("MT jobs requested ............ : %20ld\n", mt.start_queue_head);
    printf("MT jobs completed ............ : %20ld\n", mt.done_queue_head);
    printf("\n");
}

void mt_terminate(void)
{
    struct ring_entry_t job;
    // gracefully stop the worker loops
    for (unsigned i = 0; i < MT_THREAD_COUNT; i++)
    {
        job.n = 0;
        job.R_fermat = 0;
        job.R_cubic = 0;
        mod_multithread_notify_ready(&mt, job);
    }

    // wait for the worker thread terminations
    for (unsigned i = 0; i < MT_THREAD_COUNT; i++)
    {
        pthread_join(mt_tids[i], 0);
    }

    drain_done_queue(&mt);

    // clean resources
    sem_destroy(&mt.start_sem);
    sem_destroy(&mt.done_sem);
}

void mt_verify_all_a(uint64_t n, uint64_t R_fermat, uint64_t R_cubic)
{
#if 1
    // best effort, check there is possibly something in the response queue
    drain_done_queue(&mt);

    // enqueue a new request
    struct ring_entry_t job;

    job.n = n;
    job.R_fermat = R_fermat;
    job.R_cubic = R_cubic;
    mod_multithread_notify_ready(&mt, job);
#else
    // debug : ino work for worker threads, run in the thread context
    verify_all_a(n, R_fermat, R_cubic);

#endif
}

// -----------------------------------------------------------------
//
// entry point
//
// -----------------------------------------------------------------
int main(int argc, char **argv)
{
    uint64_t n_max = 1000000000000ul;
    uint64_t n_start = 5ul;
    char temp_filename1[100] = "cubic_all_a.txt\0";
    char temp_filename2[100] = "cubic_all_a.bak\0";
    bool use_temp_files = true;
    uint64_t interval_seconds = 600; // 10 minutes

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-max"))
        {
            n_max = strtoull(argv[++i], 0, 0);
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-n"))
        {
            n_start = strtoull(argv[++i], 0, 0);
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-s"))
        {
            interval_seconds = strtoull(argv[++i], 0, 0);
            continue;
        }
        if (!strcmp(argv[i], "-f"))
        {
            if (read_file(argv[++i], &n_start, &n_max, &interval_seconds))
            {
                printf("Start from value n=%ld, to max=%ld\n", n_start, n_max);
            }
            use_temp_files = false;
            continue;
        }
        if (!strcmp(argv[i], "-t"))
        {
            sprintf(temp_filename1, "%s.txt", argv[++i]);
            sprintf(temp_filename2, "%s.txt", argv[i]);
            printf("Temp file names are %s %s\n", temp_filename1, temp_filename2);
            continue;
        }
        printf("-n ddd : starting number, default is %ld\n", n_start);
        printf("-max ddd : maximum value for n = p*Q, default is %ld\n", n_max);
        printf("-s ddd : interval in seconds between restart file updates, default is %ld\n", interval_seconds);
        printf("-f aaa : alternative filename to restart at some (p,Q,n_max)\n");
        printf("\n");
        printf("-t aaa : temp filename prefix %s %s\n", temp_filename1, temp_filename2);
        printf("\n");
        printf("test n = p * Q, with p prime\n");
        printf("verify cubic test for all k < (n+1)/2\n");
        printf("NO composite number can pass the test and be found to be prime ! Never, never, never .....\n");
        printf("\n");
        exit(1);
    }

    if (use_temp_files)
    {
        if (read_file(temp_filename1, &n_start, &n_max, &interval_seconds))
        {
            printf("Start from value n=%ld, to max=%ld\n", n_start, n_max);
        }
        else
        {
            if (read_file(temp_filename2, &n_start, &n_max, &interval_seconds))
            {
                printf("Start from value n=%ld, to max=%ld\n", n_start, n_max);
            }
        }
    }

    // round n to next composite, not multiple of 3, not multiple of 2
    n_start = (n_start < 25) ? 25 : (n_start | 1);
    while (n_start % 3 == 0 || uint64_is_prime(n_start))
    {
        n_start += 2;
    }

    printf("Start from value n=%ld, to max=%ld\n", n_start, n_max);

    time_t t0 = time(NULL);
    time_t d0 = time(NULL);
    uint64_t n = n_start;
    uint64_t dn = n_start % 6 == 1 ? 4 : 2;
    uint64_t display = 0;
    modexp_count = 0;
    exponentiate_count = 0;
    worked_squarefree_count = 0;
    skip_squarefree_count = 0;
    worked_composite_count = 0;
    skip_composite_count = 0;
    worked_semiprime_count = 0;
    skip_semiprime_count = 0;

    mt_initialize();

    // temporary vector of factors
    factor_v f;
    f.reserve(64); // preallocate the vector, no further reallocations

    // update the restart file
    if (!write_file(temp_filename1, temp_filename2, n, n_max, interval_seconds))
    {
        printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
    }

    while (n <= n_max)
    {
        // master thread progress (in order)
        if (++display > 10000)
        {
            time_t d1 = time(NULL);
            printf("start n = %ld, %ld\n", n, d1 - d0);
            display = 0;
            d0 = d1;
        }

        f.clear();
        uint64_all_factors(f, n);
        if (is_prime(f))
        {
            // ------------------------------------------------------------------------------
            // prime, we are not interested with.
            // ------------------------------------------------------------------------------
        }
        else if (is_semiprime(f))
        {
            // ------------------------------------------------------------------------------
            // 2 factors have multiplicity 1   p < q
            //
            // input number is squarefree, but processing is greatly simplified when there are
            // only 2 factors
            // ------------------------------------------------------------------------------
            uint64_t p = f[0].prime;
            uint64_t Q = f[1].prime;
            uint128_t p3 = (uint128_t)p * p * p;
            uint128_t Q3 = (uint128_t)Q * Q * Q;
            uint128_t g = uint128_gcd(p3 - 1, Q3 - 1);
            uint64_t R = g > (n - 1) ? (n - 1) : (n - 1) % g;
            uint128_t An = (uint128_t)n * (n + 1) + 1;
            uint64_t A = g > An ? An : An % g;
            if (!(R < 3 || (p != 5 && R == 4) || (p < V_COUNT && R < V[p]) || (R < p - 1 && p % 10 != 1) ||
                  (Q < V_COUNT && R < V[Q]) || (R < Q - 1 && Q % 10 != 1) || A < 4))
            {
                worked_semiprime_count += 1;
                mt_verify_all_a(n, R, R);
            }
            else
            {
                skip_semiprime_count += 1;
            }
        }
        else if (is_squarefree(f))
        {
            // ------------------------------------------------------------------------------
            // All factors have multiplicity 1
            // ------------------------------------------------------------------------------
            uint64_t Rmin = n - 1;
            uint128_t An = (uint128_t)n * (n + 1) + 1;
            uint64_t i;
            for (i = 0; i < f.size(); i++)
            {
                uint64_t p = f[i].prime;
                uint64_t Q = n / p;
                uint128_t p3 = (uint128_t)p * p * p;
                uint128_t Q3 = (uint128_t)Q * Q * Q;
                uint128_t g = uint128_gcd(p3 - 1, Q3 - 1);
                uint64_t R = g > (n - 1) ? (n - 1) : (n - 1) % g;
                uint64_t A = g > An ? An : An % g;
                if (R < 3 || (p != 5 && R == 4) || (p < V_COUNT && R < V[p]) || (R < p - 1 && p % 10 != 1) || A < 4)
                {
                    // early terminate this number
                    break;
                }
                if (R < Rmin)
                {
                    Rmin = R;
                }
            }
            if (i == f.size())
            {
                // no early termination
                worked_squarefree_count += 1;
                mt_verify_all_a(n, Rmin, Rmin);
            }
            else
            {
                skip_squarefree_count += 1;
            }
        }
        else
        {
            // ------------------------------------------------------------------------------
            // general case : n is composite
            // iterate over all prime proper factors p of n
            // ------------------------------------------------------------------------------
            uint64_t Rmin = n - 1;
            uint64_t i;
            for (i = 0; i < f.size(); i++)
            {
                uint64_t p = f[i].prime;
                uint64_t Q = n / p;
                uint128_t p3 = (uint128_t)p * p * p;
                uint128_t Q3 = (uint128_t)Q * Q * Q;
                uint128_t g = uint128_gcd(p3 - 1, Q3 - 1);
                uint64_t R = g > (n - 1) ? (n - 1) : (n - 1) % g;
                if (R < 3 || (p != 5 && R == 4) || (p < V_COUNT && R < V[p]) || (R < p - 1 && p % 10 != 1))
                {
                    break;
                }
                if (R < Rmin)
                {
                    Rmin = R;
                }
            }
            if (i != f.size())
            {
                skip_composite_count += 1;
            }
            else
            {
                worked_composite_count += 1;
                mt_verify_all_a(n, n - 1, Rmin); // TODO : is it 0 or n-1 ???
            }
        }

        // after 30 minutes, update the temp files.
        // In case of crash, just restarting the program without parameter would restart from the last saved values
        time_t t1 = time(NULL);
        if (t1 - t0 >= interval_seconds)
        {
            if (!write_file(temp_filename1, temp_filename2, n, n_max, interval_seconds))
            {
                printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
            }
            else
            {
                t0 = t1;
            }
        }

        // next n, not a multiple of 3, not a multiple of 2
        n += dn;
        dn = 6 - dn;
    }

    mt_terminate();

    // update the restart file after the last iterations
    if (!write_file(temp_filename1, temp_filename2, n, n_max, interval_seconds))
    {
        printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
    }

    // a shell script could look like
    //    ./all_a -n 50000000 -s 600
    //    while $? != 0
    //       do
    //       ./all_a
    //       done

    // debug : verifies visually the number of "modexps" and "cubic exponentiates" and other values
    printf("modexp count ................. : %20ld\n", modexp_count);
    printf("cubic exponentiate count ..... : %20ld\n", exponentiate_count);
    printf("worked semiprime count ....... : %20ld\n", worked_semiprime_count);
    printf("skip semiprime count ......... : %20ld\n", skip_semiprime_count);
    printf("worked squarefree count ...... : %20ld\n", worked_squarefree_count);
    printf("skip squarefree count ........ : %20ld\n", skip_squarefree_count);
    printf("worked composite count ....... : %20ld\n", worked_composite_count);
    printf("skip composite count ......... : %20ld\n", skip_composite_count);
    printf("\n");

    // debug : verifies visually the number of multithread jobs done
    mt_debug_display();

    return (0);
}

// -----------------------------------------------------------------------
