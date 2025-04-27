# CubicPrimalityTest

WORK IN PROGRESS !!!

Test code to verify a cubic primality test, based on linear recurrences

The official paper describing the test and its proofs is work-in-progress

So far, no counterexample (pseudo prime) has been found. i

# Simple utility based of GMP library for large integers

```
$ make

$ make check

$ ./cubic 2^127-1 2^127+0 2^127+1
2^127-1 might be prime, time=       0.151 msecs.
2^127+0 is composite for sure, time=       0.000 msecs.
2^127+1 is composite for sure, time=       0.000 msecs.

$ ./cubic 0x988a04da39838a3757afef4ae6ed84b092aa0ee673067e52140862e5d27af3adfd1d65489e91b068df21f5de5e78fe4a8deb967201c7944b0a0eabc31bb0b824d3cb6293156c0c84bc48072952f08711da7a8786050335f82ec0bba57adf9c22aad36ba2f4919a3ccd8a4717799d90ffc82189f5425a3026de65b4c7e11e9beb
0x988a04da39838a3757afef4ae6ed84b092aa0ee673067e52140862e5d27af3adfd1d65489e91b068df21f5de5e78fe4a8deb967201c7944b0a0eabc31bb0b824d3cb6293156c0c84bc48072952f08711da7a8786050335f82ec0bba57adf9c22aad36ba2f4919a3ccd8a4717799d90ffc82189f5425a3026de65b4c7e11e9beb might be prime, time=       6.194 msecs.

$ ./cubic 2^11213-1
2^11213-1 might be prime, time=    3147.283 msecs.

```

So far, no counterexample (pseudo prime) has been found. Exhaustive tests completed to 10^15

# Distributed test : Quick user's guide :

The test code implements a distributed architecture to spread the work over many cores of multiple computers.
Results are collected by the server. Exhaustive tests completed to 10^15 (around 49 bits). This code as-it will not work for numbers > 62 bits.

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





