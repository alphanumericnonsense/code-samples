#-----------------------------------------------------------------------------------#
# KYBER VERSION 3.02
#
# SPECS:  https://pq-crystals.org/kyber/resources.shtml
#
# Reference implementation:
# https://github.com/pq-crystals/kyber
#
# Follows the specs and passes KAT from the reference implementation;
# Beware sample order for ref and specs differs in one spot which affects KAT.
#
# Another python implementation:
# https://github.com/jack4818/kyber-py
#
# Below is most of the work, the CPAPKE.
#-----------------------------------------------------------------------------------#

import CompactFIPS202 as FIPS # keccak/sha3/shake NEED THIS AVAILABLE IN CWD!
from random import random # to get seeds for keccak and random message testing
from copy import copy
from math import floor

# Global constants
N = 256
Q = 3329

ZETA = 17 # primitive 256th root of unity modulo Q
ZETAINV = 1175 # inverse of ZETA modulo Q

# table of powers of zeta for whenever
ZETA_POWERS = {i : ZETA**i % Q if i >=0 else ZETAINV**(-i) % Q for i in range(-N+1,N)}

# used in the iterative fast NTT, NTT_iter(), INTT_iter()
NTT_iter_pows = {i : ZETA_POWERS[128//(2**i)] for i in range(7)}
INTT_iter_pows = {i : ZETA_POWERS[-(128//(2**i))] for i in range(7)}

# Fixed length for bitstream... not sure how big to make this.
# For parse() to create integers mod Q from bytes, 19% chance of failure per trial.
# Should really be a stream, but the Python Keccak implementation has fixed
# output length parameter
# can now fix, fix parse/XOF???
XOF_OUTPUT_LENGTH = 4096 # need 3 bytes for 2 integers in [0,3329) plus failures???

# NOTE TODO:  barrett reduction?, clean up NTT?, XOF/parse as generator

#--------------------------------------------------------------------------------------#
# key generation, encryption, decryption,
# export keys, and parameter generation
#--------------------------------------------------------------------------------------#
def param_dict(level = 5):
    """
    return a dictionary of parameters for security level in [1,3,5]
    """
    ret_dict = {}
    if level == 1:
        ret_dict['k'] = 2 # rank of free module over R_Q
        ret_dict['eta1'] = 3 # parameter for random sampling
        ret_dict['eta2'] = 2 # another parameter for random sampling
        ret_dict['du'] = 10 # for compression/noise in the ciphertext
        ret_dict['dv'] = 4 # for compression/noise in the ciphertext
    elif level == 3:
        ret_dict['k'] = 3 # rank of free module over R_Q
        ret_dict['eta1'] = 2 # parameter for random sampling
        ret_dict['eta2'] = 2 # another parameter for random sampling
        ret_dict['du'] = 10 # for compression/noise in the ciphertext
        ret_dict['dv'] = 4 # for compression/noise in the ciphertext
    elif level == 5:
        ret_dict['k'] = 4 # rank of free module over R_Q
        ret_dict['eta1'] = 2 # parameter for random sampling
        ret_dict['eta2'] = 2 # another parameter for random sampling
        ret_dict['du'] = 11 # for compression/noise in the ciphertext
        ret_dict['dv'] = 5 # for compression/noise in the ciphertext
    else:
        print("Specify level 1, 3, or 5 - returning None")
        return None
    return ret_dict

def pke_key_gen(d = None, level = 5):
    """
    returns (sk, pk) at security level in [1,3,5]
    """
    params = param_dict(level)
    if d is None:
        d = bytearray([int(256*random()) for i in range(32)]) # get 32 random bytes
    Gd = G_hash(d) # SHAKE512
    rho, sigma = Gd[:32], Gd[32:] # ???
    n = 0 # counter for randomness
    A_hat = MATRIX_ZERO(params['k'])
    for i in range(params['k']):
        for j in range(params['k']):
            A_hat[i][j] = parse(XOF(rho, bytearray([j]), bytearray([i]))) ###???
    s = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        s[i] = CBD(PRF(sigma, bytearray([n]), 64*params['eta1']), params['eta1'])
        n += 1
    e = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        e[i] = CBD(PRF(sigma, bytearray([n]), 64*params['eta1']), params['eta1'])
        n += 1
    s_hat = NTT(s, params['k'])
    e_hat = NTT(e, params['k'])
    t_hat = vector_add(NTT_mat_vec_mult(A_hat, s_hat, params['k']), e_hat, params['k'])

    pk = bytearray(0)
    for i in range(params['k']):
        pk += encode(t_hat[i], 12)
    pk += rho

    sk = bytearray(0)
    for i in range(params['k']):
        sk += encode(s_hat[i], 12)
    return (pk, sk)

def pke_encrypt(pk, m, coins = None, level = 5):
    """
    params: m - plaintext message bytearray (32 bytes)
            pk - public key byte array
    return: concatenated, compressed bytearray ciphertext

    Coins are generated inside if not passed as a parameter as in specs.
    """
    params = param_dict(level)
    if coins == None:
        coins = bytearray([int(256*random()) for i in range(32)]) # get 32 random bytes
    t_part = pk[:12*32*params['k']] # bytes for t_hat
    rho = pk[12*32*params['k']:] # rest of pk, seed for A_hat

    t_hat = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        t_slice = t_part[12*32*i : 12*32*(i + 1)]
        t_hat[i] = decode(t_slice, 12)

    A_hat_trans = MATRIX_ZERO(params['k'])
    for i in range(params['k']):
        for j in range(params['k']):
            A_hat_trans[i][j] = parse(XOF(rho, bytearray([i]), bytearray([j])))

    n = 0 # counter for randomness

    r = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        r[i] = CBD(PRF(coins, bytearray([n]), 64*params['eta1']), params['eta1'])
        n += 1

    e1 = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        e1[i] = CBD(PRF(coins, bytearray([n]), 64*params['eta2']), params['eta2'])
        n += 1

    e2 = CBD(PRF(coins, bytearray([n]), 64*params['eta2']), params['eta2'])

    r_hat = NTT(r, params['k'])
    u = vector_add(INTT(NTT_mat_vec_mult(A_hat_trans, r_hat, params['k']), params['k']), e1, params['k'])
    v = ring_add(ring_add(INTT(NTT_dot(t_hat, r_hat, params['k']), params['k']), e2), decompress(decode(m, 1), 1))

    c1 = bytearray(0)
    for ui in u:
        c1 += encode(compress(ui, params['du']), params['du'])
    c2 = encode(compress(v, params['dv']), params['dv'])
    return c1 + c2

def pke_decrypt(sk, c, level = 5):
    """
    params: c - concatenated, compressed bytearray ciphertext
            sk - secret key bytearray
    return: m - plaintext message bytearray
    """
    params = param_dict(level)
    L = params['du']*32*params['k']
    c1 = c[:L]
    c2 = c[L:]
    u = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        c1_slice = c1[32*params['du']*i : 32*params['du']*(i + 1)]
        u[i] = decompress(decode(c1_slice, params['du']), params['du'])
    v = decompress(decode(c2, params['dv']), params['dv'])
    s_hat = VECTOR_ZERO(params['k'])
    for i in range(params['k']):
        sk_slice = sk[32*12*i : 32*12*(i + 1)]
        s_hat[i] = decode(sk_slice, 12)
    m_bits = compress(ring_sub(v,INTT(NTT_dot(s_hat, NTT(u, params['k']), params['k']), params['k'])), 1)
    m = encode(m_bits, 1)
    return m

def pke_export_keys(pk, sk, name, level = 5):
    """
    Export keys to "name-pke-sec" and "name-pke-pub"
    """
    with open(name + f"-pke{level}-sec", "wb") as f:
        f.write(sk)
    with open(name + f"-pke{level}-pub", "wb") as f:
        f.write(pk)

#--------------------------------------------------------------------------------------#
# Encodings/decodings, keccak pieces (hashes, PRF, XOF)
#--------------------------------------------------------------------------------------#
def G_hash(B):
    return FIPS.SHA3_512(B)

def H_hash(B):
    return FIPS.SHA3_256(B)

def XOF(B, B1, B2):
    """
    Want to return bytestream, but must fix output length at the moment.
    Used to generate almost uniform integers in [0,3329) from bytes.
    """
    return FIPS.SHAKE128(B + B1 + B2, XOF_OUTPUT_LENGTH) # ???

def PRF(B_32, B0, outlength):
    return FIPS.SHAKE256(B_32 + B0, outlength) # ???

def KDF(B, outlength):
    """
    return length outlength byte array obtained by SHAKE256ing B
    """
    return FIPS.SHAKE256(B, outlength)

def encode(f, l):
    """
    R_Q element to byte array of length 32*l

    params: f in R_Q
            integer l
    return: bytearray
    """
    #print(f)
    bitlist = []
    for i in range(N):
        bitlist += to_bits(f[i], l)
    B = bits_to_bytes(bitlist)
    #print(bitlist)#####
    return B

def to_bits(n, l):
    """
    bit representation of n of length l, little-endian
    """
    #print(n)
    bitlist = []
    for i in range(l):
        bitlist.append(n % 2)
        n = n//2
    return bitlist

def decode(B, l):
    """
    Byte array to degree 256 polynomial with entries in [0,2^l)
    """
    f = RING_ZERO()
    bitlist = bytes_to_bits(B)
    for i in range(N):
        f[i] = sum([bitlist[i*l + j]*(2**j) for j in range(l)])
    return f

def bits_to_bytes(bits):
    """
    little-endian
    """
    L = len(bits)
    out_array = bytearray(L//8)
    for i in range(0, L, 8):
        S = 0
        for j in range(8):
            print
            S += bits[i+j]*(2**j) # i is a multple of 8
        out_array[i//8] = S
    return out_array

def bytes_to_bits(B):
    """
    little-endian
    """
    bitlist = []
    for b in B:
        x = int(b)
        for i in range(8):
            bitlist.append(x % 2)
            x = x//2
    return bitlist

def parse(B):
    """
    produces uniform elements of NTT domain R_Q from a bytestream
    How long should bytestream be???????
    """
    i = 0
    j = 0
    a_hat = RING_ZERO()
    while j < N:
        d1 = B[i] + 256*(B[i+1] % 16)
        d2 = B[i+1]//16 + 16*B[i+2]
        if d1 < Q:
            a_hat[j] = d1
            j += 1
        if d2 < Q and j < N:
            a_hat[j] = d2
            j += 1
        i += 3
    return a_hat

#------------------------------------------#
# convenient initializers for zero objects
#------------------------------------------#

def RING_ZERO():
    """
    return: list of zeros representing zero in R_Q
    """
    return [0 for i in range(N)]

def VECTOR_ZERO(K):
    """
    return: list of RING_ZERO() of size K
    """
    return [RING_ZERO() for i in range(K)]

def MATRIX_ZERO(K):
    """
    return K-by-K matrix of RING_ZERO (list of rows)
    """
    return [VECTOR_ZERO(K) for i in range(K)]

#------------------------------------------------------------------------------------------------#
# some suporting functions: ring_add, ring_sub, vector_add, vector_sub, compress,
# decompress, nearestintup, CBD
#------------------------------------------------------------------------------------------------#

# not used
def barrett_reduce(a):
    """
    Barrett reduction, 1/Q approx. 5039/2^24
    Intermediate steps fit into 32 bits
    Verified for |a| < (Q-1)^2,
    i.e. handles reducing multiplication.
    """
    b = a - ((5039*(a>>10))>>14)*Q
    if b >= Q:
        return b - Q
    elif b >= 0:
        return b
    else:
        return b + Q

def ring_add(a,b):
    """
    params: a, b in R_Q
    return: c = a + b in R_Q
    """
    c = RING_ZERO()
    for i in range(N):
        c[i] = a[i] + b[i]
        c[i] = c[i] % Q
    return c

def ring_sub(a,b):
    """
    params: a, b in R_Q
    return: c = a - b in R_Q with coefficients in [0,Q)
    """
    c = RING_ZERO()
    for i in range(N):
        c[i] = (a[i] - b[i]) % Q
    return c

def vector_add(a, b, K):
    """
    params: a, b in R_Q^K
            K module dimension
    return: vector addition of a and b
    """
    c = VECTOR_ZERO(K)
    for i in range(K):
        c[i] = ring_add(a[i], b[i])
    return c

def vector_sub(a, b, K):
    """
    params: a, b in R_Q^K
            K module dimension
    return: vector subtraction a - b
    """
    c = VECTOR_ZERO(K)
    for i in range(K):
        c[i] = ring_sub(a[i], b[i])
    return c

def compress(x,d):
    """
    First half of some deterministitc rounding; used to compress ciphertext during encryption,
    "decompressed" back to R_Q before decryption

    params: integer d and x in R_Q
    return: 0 <= y[i] < 2^d
    """
    y = RING_ZERO()
    for i in range(N):
        y[i] = nearestintup(((2.0**d) * x[i])/Q) % 2**d
    return y

def decompress(x,d):
    """
    Decompresses ciphertext before decryption; second half of some deterministic error

    params: integer d and integer list x of length N.
    return: y in R_Q
    """
    y = RING_ZERO()
    for i in range(N):
        y[i] = nearestintup((Q * x[i])/2.0**d) % Q
    return y

def nearestintup(x):
    """
    Round real number to the nearest integer, rounding up if half an integer

    params: real x
    return: nearest integer to x (rounded up if necessary)
    """
    y = floor(x)
    if x - y < 0.5:
        return int(y)
    else:
        return int(y + 1)

def CBD(B, eta):
    """
    params: B is 64*eta bytes in a bytearray
            eta is a small integer
    """
    bitlist = bytes_to_bits(B)
    f = RING_ZERO()
    for i in range(N):
        a = sum([bitlist[2*i*eta + j] for j in range(eta)])
        b = sum([bitlist[2*i*eta + eta + j] for j in range(eta)])
        f[i] = (a - b) % Q
    return f

#------------------------------------------------------------------------------------------------#
# NTT stuff below...
#------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------#
# iterative NTT parts
#-----------------------------------------------------------#
def NTT_iter(f):
    """
    iterative fast NTT "in place" to be applied to even/odd parts
    OUTPUT IS NOW BIT REVERSED (decimation in frequencey without final bit-reverse correction).

    more appropriate for hardware
    """
    # decimation in frequency
    m = 256
    n = 128
    zeta = 17
    fhat = copy(f)
    # without the line below, usual FFT; powers come up in the NTT specific to Kyber
    fhat = [fhat[i]*ZETA_POWERS[i] % Q for i in range(128)] # precompose with diag(zeta^i)
    for l in range(7,0,-1):
        m = m//2
        mhalf = m//2
        zeta = zeta**2 % Q # zeta_m
        omega = 1 # holds powers of zeta_m
        for j in range(mhalf):
            for r in range(0,n,m):
                u = fhat[r + j]
                v = fhat[r + j + mhalf]
                fhat[r + j] = (u + v) % Q
                fhat[r + j + mhalf] = (u - v)*omega % Q
            omega = omega*zeta # increase power of zeta_m
    #fhat = [fhat[bit_reversal(i)] for i in range(N//2)]
    return fhat

def INTT_iter(fhat):
    """
    iterative fast inverse NTT "in place" to be applied to even/odd parts

    decimation in time; input is assumed bit-reversed, output in normal order
    the 1/128 is taken care of in parent INTT call

    more appropriate for hardware
    """
    f = copy(fhat)
    # = [f[bit_reversal(i)] for i in range(128)] # index bit-reversal
    m = 1
    for i in range(7):
        M = 2*m
        for j in range(0,128,M):
            omega = 1
            for k in range(m):
                a = f[k+j]
                b = omega*f[k+j+m] % Q
                f[k+j] = (a + b) % Q
                f[k+j+m] = (a - b) % Q
                omega = omega*INTT_iter_pows[i] % Q
        m = M
    # without the line below, usual FFT; powers come up in the NTT specific to Kyber
    f = [f[i]*ZETA_POWERS[-i] % Q for i in range(128)] # postcompose with diag(zeta^(-i))
    return f

#-----------------------------------------------------------------------------------------#
# NTT and inverse, applies some FFT to the even and odd parts separately.
# Has bit-reversed intermediate order as in the specs.
#-----------------------------------------------------------------------------------------#
def NTT(f, K): # can change to mode = "recursive", "iterative", "matrix"
    """
    Chinese remainder theorem/DFT, takes f in R_Q to fhat
    with fhat entries f mod x^2 - ZETA^(2*br(i)+1) = (fhat[2*i], fhat[2*i+1]).
    Output has bit-reversed indices (br(i) above).

    Can be applied to vectors.

    params: K, integer module dimension
            f in R_Q or R_Q^K
    """
    if len(f) == K: # vector
        fhat = []
        for l in range(K):
            fhat.append(NTT(f[l], K))
        return fhat
    elif len(f) == N: # R_Q
        feven = copy(f[::2])
        fodd = copy(f[1::2])
        fhateven = NTT_iter(feven)
        fhatodd = NTT_iter(fodd)
        fhat = [fhateven[i//2] if i % 2 == 0 else fhatodd[i//2] for i in range(N)]
        return fhat
    else:
        print("in NTT(), shouldn't happen!")
        return None

def INTT(fhat, K):
    """
    Invert the Chinese remainder theorem/Fourier transform
    Can be applied to vectors or ring elements.

    params: K, integer module dimension
            fhat in R_Q or R_Q^K
    """
    if len(fhat) == K: # vector
        f = []
        for l in range(K):
            f.append(INTT(fhat[l], K))
        return f
    elif len(fhat) == N: # R_Q
        fhateven = copy(fhat[::2])
        fhatodd = copy(fhat[1::2])
        feven = INTT_iter(fhateven)
        fodd = INTT_iter(fhatodd)
        # 3303 = 1/128 mod 3329
        f = [3303*feven[i//2] % Q if i % 2 == 0 else 3303*fodd[i//2] % Q for i in range(N)]
        return f
    else:
        print("in INTT(), shouldn't happen!")
        return None

#-----------------------------------#
# other NTT stuff
#-----------------------------------#
def NTT_mult(fhat, ghat):
    """
    Multiplies fhat and ghat in the NTT domain (pointwise)
    Assumes NTT domain index is bit reversed.

    return: hhat in NTT domain
    """
    hhat = RING_ZERO()
    for i in range(N//2):
        j = bit_reversal(i)
        hhat[2*i] = (fhat[2*i]*ghat[2*i] + ZETA**(2*j+1)*fhat[2*i+1]*ghat[2*i+1]) % Q
        hhat[2*i+1] = (fhat[2*i]*ghat[2*i+1] + fhat[2*i+1]*ghat[2*i]) % Q
    return hhat

def NTT_dot(ahat, bhat, K):
    """
    params: K, integer module dimension
            a, b are length K vectors over R_Q
    return: dot product S = <a,b> in R_Q
    """
    Shat = RING_ZERO()
    for i in range(K):
        Shat = ring_add(Shat, NTT_mult(ahat[i], bhat[i]))
    return Shat

def NTT_mat_vec_mult(Ahat, xhat, K):
    """
    params: K, integer module dimension
            R_Q matrix A, R_Q vector x
    return: y = A*x in R_Q^K
    """
    yhat = VECTOR_ZERO(K)
    for i in range(K):
        yhat[i] = NTT_dot(Ahat[i], xhat, K)
    return yhat

def bit_reversal(i):
    """
    7-bit bit-reversal of i, 0 <= i < 128
    to artificially reorder the NTT according to specs,
    somewhat silly in software implementations.
    E.g. bit_reversal(25) = 76, 0011001 <-> 1001100

    params: integer 0 <= i < 128
    return: bit-reversed integer 0 <= i < 128

    NOT USED, just tack onto NTT() and NTT_inv() if wanted
    """
    s = bin(i)[2:]
    s = '0'*(7-len(s)) + s
    result = 0
    for j in range(len(s)):
        result += int(s[j])*2**j
    return result

#------------------------------------------------------------------------------------------------#
# For testing
#------------------------------------------------------------------------------------------------#

def random_message():
    """
    Return a random message, 32-byte bytearray
    """
    return bytearray([int(256*random()) for i in range(32)])

def random_R_Q():
    return [int(3329*random()) for i in range(N)]
