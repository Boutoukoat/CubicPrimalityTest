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
$ python generate.py -b 512
0xccd61a80f0269cee32975719d65b90bdefd746aade75f9e08434169aef47d76ef83b5fcdf1e93c2df29d290bac2bfdbc2365b5e6a831888be16f89f0f195e9eb is prime
0xcb25499b2b6d7d6b2a96a97034bb2c139471d9a3e2b024df2504a03639506d45aac10f990e3966842742d5d92fb9353d9a8bab9631006706e5c7df9ff6e10f53 is prime
0x931353d39e8454fb6d2916ea1fbdab2c81b97852b4c1f0a7da827e430e2d26a5ce55462cbd5739a6d0c49ecb97e79625d27b93c7a977b2ddd4931a562a8e66a3 is prime
0xfed3b602f33185169a277ae59e9826cb180fbafbceb1a4e309b049549cd9fd3680805f20a60c1c4662d3e9491c8b50c5f37e92f5798cc77ea6056edfee64360f is prime
0xfaff86b94b65227588854a6ccec496694ab647806678a3f841f644e1d79d6e96adf01b7203b7dcce82df5ae98d7c448420826af2832144e3ef5f1d53f7edf2d3 is prime
0xf9c5209cd6c7cdb0f6ca05becfbab08d8b42a4e00fe8c57ca620097e3c93ab1af679618273e07655f295be26ace1b4825596f8adf0778a17be72384db75dc047 is prime
0xf88aba80622a78ec650ec110d0b0cab1cbcf023fb958e7010a49ce1aa189e79f3f02a792e4090fdd624c2163cc4724808aab86695dcdcf4b8d85534776cd8dbb is prime
.....

```

# Limits

Python, underlying GMP and memory limits apply.





