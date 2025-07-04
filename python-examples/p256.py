#----------------------------------------------------------------------------#
# ECDSA with NIST p-256 with SHA2-256
# https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf
# https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf
#
# some SCA countermeasures:
#   - random projectivization before each point scaling
#   - Montgomery ladder
#   - multiplicative blinding of nonce inversion, inv(k) = blind * inv(k * blind)
#   - additive scalar blinding, sP = (r+s)P + (n-r)P
#   - same code for double/add in projective coords
#----------------------------------------------------------------------------#

import os
from Crypto.Hash import SHA256

# global constants
p256    = 0xffffffff_00000001_00000000_00000000_00000000_ffffffff_ffffffff_ffffffff
G_order = 0xffffffff_00000000_ffffffff_ffffffff_bce6faad_a7179e84_f3b9cac2_fc632551
G_x     = 0x6b17d1f2_e12c4247_f8bce6e5_63a440f2_77037d81_2deb33a0_f4a13945_d898c296
G_y     = 0x4fe342e2_fe1a7f9b_8ee7eb4a_7c0f9e16_2bce3357_6b315ece_cbb64068_37bf51f5
G_z     = 1
G       = (G_x, G_y, G_z)
coeff_a = 0xffffffff_00000001_00000000_00000000_00000000_ffffffff_ffffffff_fffffffc
coeff_b = 0x5ac635d8_aa3a93e7_b3ebbd55_769886bc_651d06b0_cc53b0f6_3bce3c3e_27d2604b

#--------------------------#
# hash and modular inverse
#--------------------------#

def hash(M):
    hash_obj = SHA256.new(data=M)
    return int.from_bytes(hash_obj.digest(), 'big') % G_order

def binary_EEA(a,b):
    """
    Return (u,v,g) with u*a+v*b = g = gcd(a,b)
    Throughout, a*si+b*ti = ri for i = 0,1.
    Assumes a, b >= 0.

    "Binary" in that it works with odd/even instead of costlier integer division.
    """
    k = 0
    s0 = 1
    s1 = 0
    t0 = 0
    t1 = 1

    # trivial cases
    if a == 0:
        return a
    elif b == 0:
        return a

    # pull out common power of 2
    while not ((a|b)&1):
        a >>= 1
        b >>= 1
        k += 1

    # odd starting values
    r0 = a
    r1 = b

    while r0 != r1:
        if r0&1 == 1:
            if r0 >= r1:
                r0 = r0 - r1
                s0 = s0 - s1
                t0 = t0 - t1
            else:
                r1 = r1 - r0
                s1 = s1 - s0
                t1 = t1 - t0
        else:
            r0 >>= 1
            if not ((s0|t0)&1):
                s0 >>= 1
                t0 >>= 1
            else:
                s0 = (s0+b)>>1
                t0 = (t0-a)>>1
    return (s0, t0, r0*(1<<k))

def modinv(a, n):
    u, v, g = binary_EEA(a, n)
    if g != 1:
        print("failure in modular inversion!")
        return None
    return u % n

#------------------------------------#
# curve helpers
#------------------------------------#

# taken from https://eprint.iacr.org/2015/1060.pdf
def weier_add(P1, P2):
    X1, Y1, Z1 = P1
    X2, Y2, Z2 = P2
    
    X3 = (X1*Y2 + X2*Y1)*(Y1*Y2 - coeff_a*(X1*Z2 + X2*Z1) - 3*coeff_b*Z1*Z2) - (Y1*Z2 + Y2*Z1)*(coeff_a*X1*X2 + 3*coeff_b*(X1*Z2 + X2*Z1) - (coeff_a**2)*(Z1*Z2))
    
    Y3 = (3*X1*X2 + coeff_a*Z1*Z2)*(coeff_a*X1*X2 + 3*coeff_b*(X1*Z2 + X2*Z1) - (coeff_a**2)*(Z1*Z2)) + (Y1*Y2 + coeff_a*(X1*Z2 + X2*Z1) + 3*coeff_b*Z1*Z2)*(Y1*Y2 - coeff_a*(X1*Z2 + X2*Z1) - 3*coeff_b*Z1*Z2)
    
    Z3 = (Y1*Z2 + Y2*Z1)*(Y1*Y2 + coeff_a*(X1*Z2 + X2*Z1) + 3*coeff_b*Z1*Z2) + (X1*Y2 + X2*Y1)*(3*X1*X2 + coeff_a*Z1*Z2)
    
    return (X3 % p256, Y3 % p256, Z3 % p256)

def weier_double(P):
    return weier_add(P,P)

def proj_random(P):
    Px, Py, Pz = P
    while True:
        z = int.from_bytes(os.urandom(32), 'big') % p256
        if z != 0:
            return ((Px*z) % p256, (Py*z) % p256, (Pz*z) % p256)

# montgomery ladder https://cr.yp.to/bib/2003/joye-ladder.pdf
def weier_scale_basic(s, P):
    R0 = (0,1,0)
    R1 = P
    sbin = f"{s:0256b}"
    for bit in sbin:
        if bit == '0':
            R1 = weier_add(R0, R1)
            R0 = weier_double(R0)
        else:
            R0 = weier_add(R0, R1)
            R1 = weier_double(R1)
    return R0
    
def weier_scale(s, P):
    # randomize coords before scaling
    P = proj_random(P)
    # additive scalar blinding
    r_found = False
    while not r_found:
        r = int.from_bytes(os.urandom(32), 'big')
        if r < G_order:
            r_found = True
    splusrP = weier_scale_basic((r + s) % G_order, P)
    negrP = weier_scale_basic(G_order - r, P)
    return weier_add(splusrP, negrP)
    

# (x,y,x) -> (x/z,y/z,1) if z != 0
def to_affine(P):
    Px, Py, Pz = P
    if Pz == 0:
        print("z-coord is zero!")
        return None
    Pzinv = modinv(Pz, p256)
    return ((Px*Pzinv) % p256, (Py*Pzinv) % p256, (Pz*Pzinv) % p256)

# check point is on curve
def on_curve(P):
    Px, Py, Pz = P
    if (Pz*Py**2 - (Px**3 + coeff_a*Px*Pz**2 + coeff_b*Pz**3)) % p256 == 0:
        return True
    return False

#-------------------------------------------#
# keygen, sign, verify
#-------------------------------------------#

def keygen():
    while True:
        d = int.from_bytes(os.urandom(32), 'big')
        if 0 < d < G_order:
            Q = weier_scale(d, G)
            return (d, Q)
    
def sign(M, d):
    while True:
        k_found = False
        while not k_found:
            k = int.from_bytes(os.urandom(32), 'big')
            if 0 < k < G_order:
                k_found = True
        R = to_affine(weier_scale(k, G))
        r = R[0] % G_order
        e = hash(M)
        # multiplicative blinding of nonce inversion
        blind_found = False
        while not blind_found:
            blind = int.from_bytes(os.urandom(32), 'big') % G_order
            if blind != 0:
                blind_found = True
        kblind = (k*blind) % G_order
        kblindinv = modinv(kblind, G_order)
        sblind = (kblindinv*(e + r*d)) % G_order
        # remove blind
        s = (blind*sblind) % G_order
        if (r != 0) and (s != 0):
            return (r,s) # r.to_bytes(32, 'big') + s.to_bytes(32, 'big')
    
def verify(M, sigma, Q):
    if not on_curve(Q):
        return False
    r = sigma[0] # int.from_bytes(sigma[0:32], 'big')
    s = sigma[1] # int.from_bytes(sigma[32:64], 'big')
    if (0 < r < G_order) and (0 < s < G_order):
        pass
    else:
        return False
    sinv = modinv(s, G_order)
    e = hash(M)
    u = (e*sinv) % G_order
    v = (r*sinv) % G_order
    A = weier_scale(u, G)
    B = weier_scale(v, Q)
    Rprime = weier_add(A, B)
    rprime = to_affine(Rprime)[0]
    if rprime == r:
        return True
    return False

def test():
    d, Q = keygen()
    print(f"secret key d:\n\t{d:064x}\n")
    print(f"public key Q (projective coordinates):\n\tQx = {Q[0]:064x}\n\tQy = {Q[1]:064x}\n\tQz = {Q[2]:064x}\n")
    M = os.urandom(32)
    print(f"random 32-byte message M:\n\t{M.hex()}\n")
    sigma = sign(M, d)
    print(f"signature r:\n\t{sigma[0]:064x}\nsignature s:\n\t{sigma[1]:064x}\n")
    verified = verify(M, sigma, Q)
    print(f"verified: {verified}\n")
    return verified
    
if __name__ == "__main__":
    passcount = 0
    for i in range(100):
        print(f"test {i+1}/{100}: keygen, sign, verify 32 random bytes...\n")
        x = test()
        if x:
            passcount += 1
    print(f"passed {passcount}/100 tests\n")
