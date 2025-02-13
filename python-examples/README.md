# Assorted Python snippets

Some random crap for the purposes of proving I'm not totally incompetent.  None of the crypto should be used for cryptographic purposes; it was written for learning purposes or to provide intermediate values for testing hardware implementations.

## Contents

- `aes.py`.  An implementation of AES256-GCM [FIPS 197](https://csrc.nist.gov/pubs/fips/197/final), [SP 800-38D](https://csrc.nist.gov/pubs/sp/800/38/d/final).
- `ASCON.py`.  An implementation of [ASCON v1.2](https://ascon.iaik.tugraz.at/files/asconv12-nist.pdf).  Recently standardized (initial public draft) with minor changes in [SP 800-232](https://csrc.nist.gov/pubs/sp/800/232/ipd).
  - `ascon_ref.py` is a third-party reference implementation to test against.
  
- `dilithium/`.  A to-the-specs version of Dilithium passing KATs ([v3.1](https://pq-crystals.org/dilithium/data/dilithium-specification-round3-20210208.pdf)).  Now standardized as ML-DSA in [FIPS 204](https://csrc.nist.gov/pubs/fips/204/final) with slight changes.  Run `kat.py` to go through known answer tests.  [`CompactFIPS202.py` from [https://keccak.team/](https://keccak.team/) and `aes256_ctr_drbg.py` mostly from [here](https://github.com/GiacomoPope/dilithium-py/blob/main/src/dilithium_py/drbg/aes256_ctr_drbg.py).]
- `kyber/`.  A to-the-specs version of Kyber passing KATs ([v3.02](https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf)).  Now standardized as ML-KEM in [FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) with slight changes.    Run `kat.py` to go through known answer tests.   [`CompactFIPS202.py` from [https://keccak.team/](https://keccak.team/) and `aes256_ctr_drbg.py` mostly from [here](https://github.com/GiacomoPope/kyber-py/blob/main/src/kyber_py/drbg/aes256_ctr_drbg.py).]
- `berlekamp_massey.py`.  An implementation of the Berlekamp--Massey algorithm for finding the smallest LFSR producing a given output.  Sort of a simple, slow Euclidean algorithm, useful for solving the "key equation" for algebraic codes of Reed--Solomon type (alternant codes or subfield subcodes of generalized Reed--Solomon codes).  [*Shift-Register Synthesis and BCH Decoding*](https://crypto.stanford.edu/~mironov/cs359/massey.pdf).
- `euclid.py`.  A couple implementations of the (extended) Euclidean algorithm, including fast binary versions for integers and polynomials over $GF(2)$.
- `gf.py`.  A class for arithmetic over $GF(2^n)$.
- `fft.py`.  Some Fast Fourier transform implementations (iterative and recursive DIT/Tukey--Cooley) along with randomized testing.  Some notes on DFT/FFT in `DFT.pdf`.
- `reed_muller.py`.  Encoding/decoding Reed--Muller codes.  Some lazy inefficiency/overcounting in the decoding.
- `reed_solomon.py`.  Encoding/decoding Reed--Solomon codes.  Running the script goes through some encoding + errors + decoding for some particular parameters $`[n,k,d=n-k+1]_q = [255,32,224]_{8}`$.
- `linear_regerssion.py`.  A simple linear regression in one variable, i.e. output is `y=a*x+b` based on a list of input pairs `(xi,yi)`.
- `sorting.py`.  A bunch of sorting algorithms (heap, quick, bubble, insertion, merge, selection).  Run script to do a randomized timing test.
- `gao_mateer.py`.  An additive fast Fourier transform in characteristic 2.  Quickly evaluates polynomials over an entire subspace of $GF(2^n)$.  Run the script to test.  Cf. [*Additive Fast Fourier Transforms over Finite Fields*](https://www.math.clemson.edu/~sgao/papers/GM10.pdf).
- A few PQC KEMs for the purposes of going through the schemes; they don't meet specifications.  Run `kemtest.py` to test all three (not so interesting to watch).
  - `frodo.py`.
  - `mceliece_full_scheme.py`.
  - `saber_scheme.py`.
- `miller_rabin.py`.  Miller--Rabin primality testing.  Running the script spits out probable primes less than a hundred and a uniformly random 128-byte probable prime, both tasks with four trials.

## Requirements

Some of the crypto testing (namely `aes256_ctr_drbg` in the Kyber/Dilithium KAT testing) uses [`pycryptodome`](https://www.pycryptodome.org/) (using the `Cryptodome` package name).
