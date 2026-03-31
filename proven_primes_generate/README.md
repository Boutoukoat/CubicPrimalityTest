# generate proven primes

This is test code to generate large proven primes, NOT FOR CRYPTO USE.

The original challenge was to build unit tests with large composite and prime numbers. As they would run every so often in a ci/cd environment, 
primality or compositeness tests done on these numbers should NOT fail the expected outcome because of pseudoprimes.

Large primes can be constructively generated from smaller primes.

# The maths

From a prime p and some number q, with p large and q small, 

r = 1 + 2 \* p \* q 

could be proved prime with a reasonable verification

e.g. well known Sophie Germain primes

let q = 1
let p = 1 mod 3 , p is a proven prime
if 2^p == n-1 mod n, the 1 + 2\*p\*q is prime

Recursively, it is possible to build larger and larger proven primes

# The maths upside-down

To verify one of the generated number is a prime number, use a p-1 algorithm

factor p1 = 1 + 2 * p2 * q2

- since q2 by construction is small this is easy to factor
- p2 is by construction a large prime
- primality of p1 can be confirmed by a p-1 algorithm

factor p2 = 1 + 2 * p3 * q3

- since q3 by construction is small this is easy to factor
- p3 is by construction a large prime
- primality of p2 can be confirmed by a p-1 algorithm

factor p3 = 1 + 2 * p4 * q4

e.g.

0x81200FCF9EA491CA7438F6458884F4D123014C4A23C2CF3B15AEA8CCFA984D6562AE080D7DCF241258384B5ED737DB59C30198034EA16380BC0346F6969F0F7E 
= 1 + 0x2 \* 0x7DAD \* 0x838364952C0F6148AC08A67B6A34A58AA57E96BACC13719D3679E363FD642D0ED5A31C42DC9CA435A9AAFC2F1A2B819622E7D9C18840C07F8CEC599F309B


0x838364952C0F6148AC08A67B6A34A58AA57E96BACC13719D3679E363FD642D0ED5A31C42DC9CA435A9AAFC2F1A2B819622E7D9C18840C07F8CEC599F309A
= 1 + 0x2 \* 0x43 \* 0xFB3FCB9EE55E568331BC775E678E9BC7EBEE053703F74BAE41DD65FFFB04283AE861A4CC2760BF71FB8B6F285090402A4A53C6386B77E2627C084F842F


and so on for p3, p4, p5 ..... until pn is small enough (e.g. < 2^64) to be confirmed as a prime.


# Simple utility based of GMP library for large integers

Accept on command line the bit size of the generated primes 


```
$ python generate.py -b 256
0x85a1b591a16576e86c0c5f0a8c27794e8e61eff458cb27df62b7bf0745cb2d53 is prime
0x82645dd06e85a0147856890d3df41e4f600f8cd5eab9c85a7733471416edb13b is prime
0x81819ffd586a7a2772ded4da463586b15549bc23fb4b1a0fa99195427875de1f is prime
0xdd61ccd6b298b6cf2ddfdd5ac8921b0c5680b9e2b67cf408397dc34c434dff0f is prime
0x8a4ae32d18d32dfa6d5a41fb86fc3f13a2056fdb2db7b98c6a8d16ca9c4b47bb is prime
0x893fa10b88c3e52ccabd90bb499ac8a86e8a11966ba3c7864f974202cf35ac2f is prime
0x856a52616639b2c673b6e6a46c2356577cbf3a633d6a2da69c7a270dbcf7ec2b is prime
0x81eab34cc5e73601cbf02736ac1f2b36978d48bd4f5e0e51a9b3ae1e196d11cf is prime
0x852b61150ce8211a701dcef017fff5d866622858f960c01b8b3d75cd1573dcef is prime
0x850194db3022495f57959a3cf6bcd8407fa0614ec24b692833afceb73b919f8b is prime
.....

```

# Limits

Python, underlying GMP and memory limits apply.





