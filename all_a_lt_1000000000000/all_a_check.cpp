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
//       if(B^2+B+1==-x^2+x+a.print([n,p,Q,g,R,A]);break(7)))))))));}
//
// -----------------------------------------------------------------------

#include <omp.h>
#include <assert.h>

#include "math_utils.cpp"

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
}

// read file, cautious about corruped or truncated files
bool read_file(char *fn, uint64_t *n, uint64_t *i, uint64_t *m, uint64_t *s)
{
    char line[100];
    bool n_found = false;
    bool i_found = false;
    bool m_found = false;
    bool s_found = false;
    bool valid = true;
    uint64_t n_temp = -1, i_temp = -1, m_temp = -1, s_temp = -1;
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
            if (!memcmp(line, "i=", 2))
            {
                i_temp = strtoull(&line[2], 0, 0);
                i_found = true;
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

    if (valid && n_found == true && i_found == true && m_found == true && s_found == true)
    {
        // file is complete, without missing field
        *n = n_temp;
        *i = i_temp;
        *m = m_temp;
        *s = s_temp;
        return true;
    }
    return false;
}

// write file, cautious about worst cases in the life of a computer
bool write_file(char *fn, char *bak, uint64_t n, uint64_t i, uint64_t m, uint64_t s)
{
    uint64_t n_temp = -1, i_temp = -1, m_temp = -1, s_temp = -1;
    // first verify the current file is a valid one
    if (read_file(fn, &n_temp, &i_temp, &m_temp, &s_temp))
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
        fprintf(f, "i=%ld\n", i);
        fprintf(f, "m=%ld\n", m);
        fprintf(f, "s=%ld\n", s);
        fclose(f);

        // never paranoid enough against full disks, need to read back the file.
        // This also flushes the file in the case of remote file (read-after-write security on NFS ....)
        n_temp = -1;
        i_temp = -1;
        m_temp = -1;
        s_temp = -1;
        if (read_file(fn, &n_temp, &i_temp, &m_temp, &s_temp))
        {
            // return true if the file is correctly written
            return (n_temp == n && i_temp == i && m_temp == m && s_temp == s);
        }
    }
    return false;
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
    uint64_t i_start = 0ul;
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
        if (!strcmp(argv[i], "-i"))
        {
            i_start = strtoull(argv[++i], 0, 0);
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
            if (read_file(argv[++i], &n_start, &i_start, &n_max, &interval_seconds))
            {
                printf("Start from values n=%ld factor# i=%ld, to max=%ld\n", n_start, i_start, n_max);
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
        printf("-i ddd : starting factor, default is %ld\n", i_start);
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
        if (read_file(temp_filename1, &n_start, &i_start, &n_max, &interval_seconds))
        {
            printf("Start from values n=%ld factor# i=%ld, to max=%ld\n", n_start, i_start, n_max);
        }
        else
        {
            if (read_file(temp_filename2, &n_start, &i_start, &n_max, &interval_seconds))
            {
                printf("Start from values n=%ld factor# i=%ld, to max=%ld\n", n_start, i_start, n_max);
            }
        }
    }

    // round n to next composite, not multiple of 3, not multiple of 2
    n_start = (n_start < 25) ? 25 : (n_start | 1);
    while (n_start % 3 == 0 || uint64_is_prime(n_start))
    {
        n_start += 2;
    }

    printf("Start from values n=%ld factor# i=%ld, to max=%ld\n", n_start, i_start, n_max);

    time_t t0 = time(NULL);
    time_t d0 = time(NULL);
    uint64_t n = n_start;
    uint64_t dn = n_start % 6 == 1 ? 4 : 2;
    uint64_t display = 0;
    uint64_t modexp_count = 0;
    uint64_t exponentiate_count = 0;

    // update the restart file
    if (!write_file(temp_filename1, temp_filename2, n, i_start, n_max, interval_seconds))
    {
        printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
    }

    while (n <= n_max)
    {
        factor_v f;
        bool is_exact_power = uint64_all_factors(f, n);
        if (is_exact_power == false)
        {
            // iterate over the prime factors of n
            for (uint64_t i = i_start; i < f.size(); i++)
            {
                // single thread progress
                if (++display > 10000)
                {
                    time_t d1 = time(NULL);
                    printf("n = %ld, i = %ld, %ld\n", n, i, d1 - d0);
                    display = 0;
                    d0 = d1;
                }

                uint64_t p = f[i].prime;
                uint64_t Q = n / p;
                if (f.size() == 2 && p == f[1].prime && Q == f[0].prime)
                {
                    // avoid double computation in the case of semiprimes
                    continue;
                }
                uint128_t p3 = (uint128_t)p * p * p;
                uint128_t Q3 = (uint128_t)Q * Q * Q;
                uint128_t g = uint128_gcd(p3 - 1, Q3 - 1);
                uint64_t R = g > (n - 1) ? (n - 1) : (n - 1) % g;
                if (R > 2)
                {
                    uint128_t A = (uint128_t)n * (n + 1) + 1;
                    A = g > A ? A : A % g;
                    if (A > 3)
                    {
                        struct barrett_t bn;
                        barrett_precompute(&bn, n);
                        //
                        // TODO : start a child thread from here
                        //
                        uint64_t a = 7;
                        for (uint64_t d = 0; d < n; d += 2)
                        {
                            // a = 7 + k * (k-1) unrolled as  a = a + 2k
                            a += d;
                            a -= (a >= n) ? n : 0;

                            // verify test Mod(a,n)^R == 1
                            modexp_count += 1;
                            if (barrett_pow(a, R, bn) == 1)
                            {
                                exponentiate_count += 1;
                                // run cubic test
                                // B = Mod(x, n)
                                uint64_t bs = 0;
                                uint64_t bt = 1;
                                uint64_t bu = 0;
                                // B = Mod(B, x^3 -ax -a)^R
                                cubic_exponentiate(bs, bt, bu, R, n, a);
                                // B2 = B
                                uint64_t bs2 = bs;
                                uint64_t bt2 = bt;
                                uint64_t bu2 = bu;
                                // B2 = Mod(B2, x^3 -ax -a)^2
                                cubic_exponentiate(bs2, bt2, bu2, 2, n, a);
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
                                            ptr += uint128_sprint(ptr, p);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, Q);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, g);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, R);
                                            *ptr++ = ',';
                                            ptr += uint128_sprint(ptr, A);
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
                }

                // after 30 minutes, update the temp files.
                // In case of crash, just restarting the program without parameter would restart from last saved values
                time_t t1 = time(NULL);
                if (t1 - t0 >= interval_seconds)
                {
                    if (!write_file(temp_filename1, temp_filename2, n, i, n_max, interval_seconds))
                    {
                        printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
                    }
                    else
                    {
                        t0 = t1;
                    }
                }
            }
        }

        // next n, not a multiple of 3, not a multiple of 2, and composite
        do
        {
            n += dn;
            dn = 6 - dn;
        } while (uint64_small_factor(n) == 1 || uint64_is_prime(n));

        // restart next number from first proper factor
        i_start = 0;
    }

    // update the restart file after the last iterations
    if (!write_file(temp_filename1, temp_filename2, n, i_start, n_max, interval_seconds))
    {
        printf("Unable to write the file for restart after crash (%s)\n", temp_filename1);
    }

    // a shell script could look like
    //    ./all_a -n 50000000 -s 600
    //    while $? != 0
    //       do
    //       ./all_a
    //       done

    // debug : verifies visually the number of "modexps" and "cubic exponentiates"
    printf("modexp count ........... : %20ld\n", modexp_count);
    printf("exponentiate count ..... : %20ld\n", exponentiate_count);

    return (0);
}
// -----------------------------------------------------------------------
