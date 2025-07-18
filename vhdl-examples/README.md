# Some VHDL examples

## Contents

- `/bubblesort`.  Parallelized bubblesort with tb and top level test (16 switches, 16 LEDs, sorts 4 4-bit switch inputs and displays on LEDs).  Two rows of offset compare/swap units with feedback on clock.  Sorts array of $N$ numbers in $N/2$ cycles.  Not practical for $N$ large or large integer width.
- `/lfsr`.  A linear-feedback shift register with tb and quick mathematical explanation of LFSR as a dynamical system.
- `/pipelined-crypto`.  Rounds are unrolled and accept new input each cycle.  No modes, just single block encryption or hashing (when applicable).
  - `/aes`.  Randomized tb tests encryption against third-party.
  - `/gift`.   TB checks the three test vectors from [here](https://github.com/giftcipher/gift/tree/master/implementations/test%20vectors).  A couple architectures (one/two rounds per cycle, round-reduced).
  - `/chacha20`.  Uses the "original" 64/64 counter/nonce split.  A couple architectures (one/two rounds per cycle, round-reduced).  Randomized tb tests encryption against third-party.
  - `/sha2`.  Nothing fancy, single block hash and tbs.  Not practical, one round per cycle on the full 256-bit state, all 64 rounds unrolled.
  - `/sha3`.  Nothing fancy, basically just the Keccak permuation, some tbs.  Not practical, one round per cycle on the full 1600-bit state, all 24 rounds unrolled.
- `/slow-crypto`.  Reusing the round function.
  - `/aes`.  TB iterates encryption with feedback, compares output to third-party reference.
  - `/gift`.  TB checks the three test vectors from [here](https://github.com/giftcipher/gift/tree/master/implementations/test%20vectors).  Two rounds per cycle.
- `/prince`.  TB checks against test vectors from [here](https://eprint.iacr.org/2012/529.pdf).  Just a function; haven't put it in hardware to see if it will actually fit in a single cycle at any reasonable frequency (PRINCE is designed as a low latency block cipher by/for NXP, used for flash encryption in their products).


## Requirements

The Makefiles, testbenches, and test-vector scripts use make (obviously), [GHDL](https://ghdl.github.io/ghdl/), and Python; change variables in Makefiles as appropriate.  Uses pycryptodome under `Cryptodome` package name for some test vector generation scripts.
