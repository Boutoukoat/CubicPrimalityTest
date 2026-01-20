# CubicPrimalityTest

WORK IN PROGRESS !!!

This is test code to verify a cubic primality test, based on linear recurrences

The official paper describing the test and its proofs is work-in-progress

So far, no counterexample (pseudo prime) has been found. 

# The maths

Let irreducible f = x^3-a\*x-b such that kronecker(4\*a^3-27\*b^2, n) == 1.

Ignore B == 1 where B = x^(n-1) (mod n, f). This ensures x not in Zn.

Check x^2\*(B^2+B+1) == a.

Linear Fermat PRP tests have x = a in Zn.

Quadratic PRP tests either have x in Zn or roots in complex numbers.

Cubic PRP tests with kronecker(discriminant. n) == 1 can be chosen to have all three roots not in Zn.

With b = a we can rapidly find a suitable B with the Mod(n,a)^((a-1)/3) != 1 test with parametric prime a = 7+k\*(k-1).

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

# Complete user's guide :

Later.





