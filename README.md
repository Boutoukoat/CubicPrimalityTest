# CubicPrimalityTest

WORK IN PROGRESS !!!

Test code to verify a cubic primality test, based on linear recurrences

The official paper describing the test and its proofs is work-in-progress

So far, no counterexample (pseudo prime) has been found. Exhaustive tests completed to 10^15 (around 49 bits). This code as-it will not work for numbers > 62 bits.

# Quick user's guide :

The test code implements a distributed architecture to spread the work over many cores of multiple computers.
Results are collected by the server.

```
$ make
```

On server computer 10.10.10.1 , run 1 server thread to distribute the work and get the results

```
$ ./lnrc -server -e 4
```

On clients on the same subnetwork, e.g. CPUs with 20 cores, run 20 threads

```
$ ./lnrc -t 20 -s 10.10.10.1
```

# Complete user's guide :

Later.





