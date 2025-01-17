#--------------------------------------------------------------------------------------------------#
# DILITHIUM VERSION 3.1
#
# SPECS: https://pq-crystals.org/dilithium/resources.shtml
#
# Follows the specs and passes KAT from reference implementation
# https://github.com/pq-crystals/dilithium
#
# Another python implementation:
# https://github.com/jack4818/dilithium-py
#--------------------------------------------------------------------------------------------------#

from random import random
from time import time
import math
import CompactFIPS202 as FIPS
import os
from packing import *

QQ = 8380417
ZETA = 1753
ZETA_POWERS = [1, 1753, 3073009, 6757063, 3602218, 4234153, 5801164, 3994671, 5010068, 8352605, 1528066, 5346675, 3415069, 2998219, 1356448, 6195333, 7778734, 1182243, 2508980, 6903432, 394148, 3747250, 7062739, 3105558, 5152541, 6695264, 4213992, 3980599, 5483103, 7921677, 348812, 8077412, 5178923, 2660408, 4183372, 586241, 5269599, 2387513, 3482206, 3363542, 4855975, 6400920, 7814814, 5767564, 3756790, 7025525, 4912752, 5365997, 3764867, 4423672, 2811291, 507927, 2071829, 3195676, 3901472, 860144, 7737789, 4829411, 1736313, 1665318, 2917338, 2039144, 4561790, 1900052, 3765607, 5720892, 5744944, 6006015, 2740543, 2192938, 5989328, 7009900, 2663378, 1009365, 1148858, 2647994, 7562881, 8291116, 2683270, 2358373, 2682288, 636927, 1937570, 2491325, 1095468, 1239911, 3035980, 508145, 2453983, 2678278, 1987814, 6764887, 556856, 4040196, 1011223, 4405932, 5234739, 8321269, 5258977, 527981, 3704823, 8111961, 7080401, 545376, 676590, 4423473, 2462444, 749577, 6663429, 7070156, 7727142, 2926054, 557458, 5095502, 7270901, 7655613, 3241972, 1254190, 2925816, 140244, 2815639, 8129971, 5130263, 1163598, 3345963, 7561656, 6143691, 1054478, 4808194, 6444997, 1277625, 2105286, 3182878, 6607829, 1787943, 8368538, 4317364, 822541, 482649, 8041997, 1759347, 141835, 5604662, 3123762, 3542485, 87208, 2028118, 1994046, 928749, 2296099, 2461387, 7277073, 1714295, 4969849, 4892034, 2569011, 3192354, 6458423, 8052569, 3531229, 5496691, 6600190, 5157610, 7200804, 2101410, 4768667, 4197502, 214880, 7946292, 1596822, 169688, 4148469, 6444618, 613238, 2312838, 6663603, 7375178, 6084020, 5396636, 7192532, 4361428, 2642980, 7153756, 3430436, 4795319, 635956, 235407, 2028038, 1853806, 6500539, 6458164, 7598542, 3761513, 6924527, 3852015, 6346610, 4793971, 6653329, 6125690, 3020393, 6705802, 5926272, 5418153, 3009748, 4805951, 2513018, 5601629, 6187330, 2129892, 4415111, 4564692, 6987258, 4874037, 4541938, 621164, 7826699, 1460718, 4611469, 5183169, 1723229, 3870317, 4908348, 6026202, 4606686, 5178987, 2772600, 8106357, 5637006, 1159875, 5199961, 6018354, 7609976, 7044481, 4620952, 5046034, 4357667, 4430364, 6161950, 7921254, 7987710, 7159240, 4663471, 4158088, 6545891, 2156050, 8368000, 3374250, 6866265, 2283733, 5925040, 3258457, 5011144, 1858416, 6201452, 1744507, 7648983, 8380416, 8378664, 5307408, 1623354, 4778199, 4146264, 2579253, 4385746, 3370349, 27812, 6852351, 3033742, 4965348, 5382198, 7023969, 2185084, 601683, 7198174, 5871437, 1476985, 7986269, 4633167, 1317678, 5274859, 3227876, 1685153, 4166425, 4399818, 2897314, 458740, 8031605, 303005, 3201494, 5720009, 4197045, 7794176, 3110818, 5992904, 4898211, 5016875, 3524442, 1979497, 565603, 2612853, 4623627, 1354892, 3467665, 3014420, 4615550, 3956745, 5569126, 7872490, 6308588, 5184741, 4478945, 7520273, 642628, 3551006, 6644104, 6715099, 5463079, 6341273, 3818627, 6480365, 4614810, 2659525, 2635473, 2374402, 5639874, 6187479, 2391089, 1370517, 5717039, 7371052, 7231559, 5732423, 817536, 89301, 5697147, 6022044, 5698129, 7743490, 6442847, 5889092, 7284949, 7140506, 5344437, 7872272, 5926434, 5702139, 6392603, 1615530, 7823561, 4340221, 7369194, 3974485, 3145678, 59148, 3121440, 7852436, 4675594, 268456, 1300016, 7835041, 7703827, 3956944, 5917973, 7630840, 1716988, 1310261, 653275, 5454363, 7822959, 3284915, 1109516, 724804, 5138445, 7126227, 5454601, 8240173, 5564778, 250446, 3250154, 7216819, 5034454, 818761, 2236726, 7325939, 3572223, 1935420, 7102792, 6275131, 5197539, 1772588, 6592474, 11879, 4063053, 7557876, 7897768, 338420, 6621070, 8238582, 2775755, 5256655, 4837932, 8293209, 6352299, 6386371, 7451668, 6084318, 5919030, 1103344, 6666122, 3410568, 3488383, 5811406, 5188063, 1921994, 327848, 4849188, 2883726, 1780227, 3222807, 1179613, 6279007, 3611750, 4182915, 8165537, 434125, 6783595, 8210729, 4231948, 1935799, 7767179, 6067579, 1716814, 1005239, 2296397, 2983781, 1187885, 4018989, 5737437, 1226661, 4949981, 3585098, 7744461, 8145010, 6352379, 6526611, 1879878, 1922253, 781875, 4618904, 1455890, 4528402, 2033807, 3586446, 1727088, 2254727, 5360024, 1674615, 2454145, 2962264, 5370669, 3574466, 5867399, 2778788, 2193087, 6250525, 3965306, 3815725, 1393159, 3506380, 3838479, 7759253, 553718, 6919699, 3768948, 3197248, 6657188, 4510100, 3472069, 2354215, 3773731, 3201430, 5607817, 274060, 2743411, 7220542, 3180456, 2362063, 770441, 1335936, 3759465, 3334383, 4022750, 3950053, 2218467, 459163, 392707, 1221177, 3716946, 4222329, 1834526, 6224367, 12417, 5006167, 1514152, 6096684, 2455377, 5121960, 3369273, 6522001, 2178965, 6635910, 731434]

def param_dict(level):
    ret_dict = {}
    if level == 2:
        # Parameters, Level 2
        ret_dict['DD'] = 13
        ret_dict['TAU'] = 39
        ret_dict['GAMMA1'] = 2**17
        ret_dict['GAMMA2'] = (QQ-1)//88 # 95232
        ret_dict['KK'] = 4
        ret_dict['LL'] = 4
        ret_dict['ETA'] = 2
        ret_dict['BETA'] = 78
        ret_dict['OMEGA'] = 80
    elif level == 3:
        # Parameters, Level 3
        ret_dict['DD'] = 13
        ret_dict['TAU'] = 49
        ret_dict['GAMMA1'] = 2**19
        ret_dict['GAMMA2'] = (QQ-1)//32 # 261888
        ret_dict['KK'] = 6
        ret_dict['LL'] = 5
        ret_dict['ETA'] = 4
        ret_dict['BETA'] = 196
        ret_dict['OMEGA'] = 55
    elif level == 5:
        # Parameters, Level 5
        ret_dict['DD'] = 13
        ret_dict['TAU'] = 60
        ret_dict['GAMMA1'] = 2**19
        ret_dict['GAMMA2'] = (QQ-1)//32 # 261888
        ret_dict['KK'] = 8
        ret_dict['LL'] = 7
        ret_dict['ETA'] = 2
        ret_dict['BETA'] = 120
        ret_dict['OMEGA'] = 75
    else:
        print("input level in [2,3,5]; returning None.")
    return ret_dict

#------------------------------#
# Scalar Helpers
#------------------------------#

def ModPM(x,alpha):
    """
    return x mod alpha in the range [-alpha/2, alpha/2], rounding up.
    """
    if alpha % 2 == 0:
        y = x % alpha
        if y > alpha/2:
            y = y - alpha
    elif alpha % 2 == 1:
        y = x % alpha
        if y > (alpha - 1)/2:
            y = y - alpha
    else:
        print("ModPM error, returning None")
        return None
    return int(y)

def Pow2Round(r,d):
    r = r % QQ
    r0 = ModPM(r, 2**d)
    return ((r - r0)>>d, r0)

def Decompose(r, alpha):
    r = r % QQ
    r0 = ModPM(r, alpha)
    if r - r0 == QQ - 1:
        r1 = 0
        r0 = r0 - 1
    else:
        r1 = (r - r0)//alpha
    return (r1, r0)

def HighBits(r, alpha):
    r1, r0 = Decompose(r, alpha)
    return int(r1)

def LowBits(r, alpha):
    r1, r0 = Decompose(r, alpha)
    return int(r0)

def MakeHint(z, r, alpha):
    r1 = HighBits(r, alpha)
    v1 = HighBits(r + z, alpha)
    if r1 == v1:
        return 0
    else:
        return 1

def UseHint(h, r, alpha):
    m = (QQ - 1)//alpha
    r1, r0 = Decompose(r, alpha)
    if h == 1:
        if r0 > 0:
            return (r1 + 1) % m
        else:
            return (r1 - 1) % m
    else:
        return int(r1)

#------------------------------#
# Ring helpers
#------------------------------#

def RqModPM(a, b, alpha):
    c = []
    for i in range(256):
        c.append(ModPM(a[i], b[i], alpha))
    return c

def RqPow2Round(a, d):
    b1 = []
    b0 = []
    for i in range(256):
        r1, r0 = Pow2Round(a[i], d)
        b1.append(r1)
        b0.append(r0)
    return (b1, b0)

def RqDecompose(a, alpha):
    b1 = []
    b0 = []
    for i in range(256):
        r1, r0 = Decompose(a[i], alpha)
        b1.append(r1)
        b0.append(r0)
    return (b1, b0)

def RqHighBits(a, alpha):
    b1, b0 = RqDecompose(a, alpha)
    return b1

def RqLowBits(a, alpha):
    b1, b0 = RqDecompose(a, alpha)
    return b0

def RqMakeHint(z, a, alpha):
    b = []
    for i in range(256):
        b.append(MakeHint(z[i], a[i], alpha))
    return b

def RqUseHint(h, a, alpha):
    b = []
    for i in range(256):
        b.append(UseHint(h[i], a[i], alpha))
    return b

#------------------------------#
# Vector helpers
#------------------------------#

def VecModPM(a, b, alpha, k):
    c = []
    for i in range(k):
        c.append(RqModPM(a[i], b[i], alpha))
    return c

def VecPow2Round(a, d, k):
    b1 = []
    b0 = []
    for i in range(k):
        r1, r0 = RqPow2Round(a[i], d)
        b1.append(r1)
        b0.append(r0)
    return (b1, b0)

def VecDecompose(a, alpha, k):
    b1 = []
    b0 = []
    for i in range(k):
        r1, r0 = RqDecompose(a[i], alpha)
        b1.append(r1)
        b0.append(r0)
    return (b1, b0)

def VecHighBits(a, alpha, k):
    b1, b0 = VecDecompose(a, alpha, k)
    return b1

def VecLowBits(a, alpha, k):
    b1, b0 = VecDecompose(a, alpha, k)
    return b0

def VecMakeHint(z, a, alpha, k):
    b = []
    for i in range(k):
        b.append(RqMakeHint(z[i], a[i], alpha))
    return b

def VecUseHint(h, a, alpha, k):
    b = []
    for i in range(k):
        b.append(RqUseHint(h[i], a[i], alpha))
    return b

#------------------------------#
# NTT
#------------------------------#

def FFT256(f, zeta_in):
    """
    cooley-tukey fft with root of unity zeta_in**2???
    """
    zeta = zeta_in
    fhat = f.copy()
    for m in [256>>i for i in range(8)]:
        zeta = zeta**2 % QQ
        omega = 1
        for j in range(m//2):
            for r in range(0,256,m):
                u = fhat[r + j]
                v  = fhat[r + j + m//2]
                fhat[r + j] = (u + v) % QQ
                fhat[r + j + m//2] = ((u - v)*omega ) % QQ
            omega = (omega*zeta) % QQ
    return fhat

def IFFT256(fhat, zeta_in):
    """
    gentleman-sande fft with root of unity zeta_in**2???,
    1/256 scaling mod 8380417=QQ
    """
    zeta = zeta_in
    f = fhat.copy()
    for m in [2**i for i in range(8)]:
        for j in range(0,256,2*m):
            omega = 1
            for k in range(m):
                a = f[k + j]
                b = omega*f[k + j + m]
                f[k + j] = (a + b) % QQ
                f[k + j + m] = (a - b) % QQ
                omega = (omega * ZETA_POWERS[512 - 256//m]) % QQ
    # scale by 1/256 = 8347681 mod QQ
    return [(8347681 * f[i]) % QQ for i in range(256)]

def NTT(f):
    """
    usual fft on ZETA^2 with prescaling
    """
    #print(type(f))
    # scale
    fscale = [(ZETA_POWERS[i]*f[i]) % QQ for i in range(256)]
    # fft
    return FFT256(fscale, ZETA)

def INTT(fhat):
    """
    usual ifft on ZETA^2 with postscaling
    """
    # ifft
    newf = IFFT256(fhat, ZETA)
    # scale
    return [(newf[i]*ZETA_POWERS[-i])%QQ for i in range(256)]

# coefficient-wise NTT/INTT for vectors and matrices
def VecNTT(v):
    #print(len(v), v[0])
    return [NTT(vi) for vi in v]

def VecINTT(vhat):
    return [INTT(vhati) for vhati in vhat]

def MatNTT(A):
    return [VecNTT(row) for row in A]

def MatINTT(Ahat):
    return [VecINTT(rowhat) for rowhat in Ahat]

#------------------------------#
# Ring and Linear algebra
#------------------------------#

def RqZero():
    """
    return zero poly
    """
    return [0 for i1 in range(256)]

def RqOne():
    """
    ring unity
    """
    return [1]+[0 for i1 in range(255)]

def VecZero(k):
    """
    return zero vector of length k
    """
    return [RqZero() for i in range(k)]

# not used
def RqMultSB(a,b):
    """
    schoolbook poly mult, x^256 = -1
    """
    c = RqZero()
    for i in range(256):
        for j in range(256):
            k = i + j
            if k < 256:
                c[k] += a[i]*b[j]
                c[k] = c[k] % QQ
            else:
                c[k-256] -= a[i]*b[j]
                c[k-256] = c[k-256] % QQ
    return c

# not used
def RqMultNTT(a,b):
    """
    multiplication in Rq using NTT
    """
    ahat = NTT(a)
    bhat = NTT(b)
    return INTT([ahat[i] * bhat[i] for i in range(256)])

def PWMult(ahat,bhat):
    """
    multiplication in NTT domain
    """
    return [(ahat[i]*bhat[i]) % QQ for i in range(256)]

def RqScale(c, r):
    """
    scale poly r by int c
    """
    return [(ri*c) % QQ for ri in r]

def VecScale(chat, vhat):
    """
    multiply a vector by c NTT DOMAIN!!!
    """
    return [PWMult(chat, vhati) for vhati in vhat]

def RqAdd(a,b):
    c = []
    for i in range(256):
        c.append((a[i] + b[i]) % QQ)
    return c

def RqSub(a,b):
    c = []
    for i in range(256):
        c.append((a[i] - b[i]) % QQ)
    return c

def VecAdd(u, v, k):
    w = []
    for i in range(k):
        w.append(RqAdd(u[i], v[i]))
    return w

def VecSub(u, v, k):
    w = []
    for i in range(k):
        w.append(RqSub(u[i], v[i]))
    return w

def MatVecMult(Ahat, vhat, k, l):
    """
    kxl matrix times length l vector, poly entries
    IN NTT DOMAIN!!!
    """
    what = []
    for i in range(k):
        s = RqZero()
        for j in range(l):
            s = RqAdd(s, PWMult(Ahat[i][j], vhat[j]))
        what.append(s)
    return what

#------------------------------#
# Norms
#------------------------------#

def InfNorm(r):
    return abs(ModPM(r, QQ))

def RqInfNorm(a):
    oldmax = 0
    for ai in a:
        norm = InfNorm(ai)
        #print(ai, norm)
        if norm > oldmax:
            oldmax = norm
    return oldmax

def VecInfNorm(v):
    oldmax = 0
    for vi in v:
        norm = RqInfNorm(vi)
        if norm > oldmax:
            oldmax = norm
    return oldmax

#------------------------------#
# Hash
#------------------------------#

def H_hash(B, outlen):
    return FIPS.SHAKE256(bytes(B), outlen)

#------------------------------#
# Expansion from seed,
# sampleinball,
# other
#------------------------------#

def VecCount1s(h, k):
    """
    count the number of 1 coefficients in the vector of polys h
    """
    count = 0
    for i in range(k):
        for j in range(256):
            if h[i][j] == 1:
                count += 1
    return count

def SampleInBall(seed, tau):
    """
    use seed for PRNG to generate output of certain form
    section 5.3, pg. 19 and section 2.3 figure 2 pg. 10.

    NOTE TODO needs to be fixed!
    """
    #print("sampleseed", ba2h(seed))
    c = RqZero()
    myRNG = FIPS.gen_SHAKE256(seed)
    signint = 0
    #signbytes = b''
    for i in range(8):
        x = next(myRNG)
        #signbytes += bytes([x])
        signint |= x << (8*i) # bytes([next(myRNG)])
    #signint = int.from_bytes(signbytes, 'little')
    #print("signbytes", ba2h(signbytes))
    #print(int.from_bytes(signbytes, 'little'))
    #print(int.from_bytes(signbytes, 'big'))
    #print("signint", signint)
    for i in range(256 - tau, 256):
        while True:
            x = next(myRNG)
            if 0 <= x <= i:
                j = x # sample [0,i]
                break;
        s = 1 - (2*(signint&1)) # sample +/- 1
        signint >>= 1
        c[i] = c[j]
        c[j] = s
    myRNG.close()
    return c

def ExpandA(rho, k, l):
    """
    output is in NTT domain!!!
    uses SHAKE128
    """
    #print("expandA", ba2h(rho))
    Ahat = [VecZero(l) for i in range(k)]
    for i in range(k):
        for j in range(l):
            XOF = FIPS.gen_SHAKE128(rho + (256*i+j).to_bytes(2,'little')) # ??? = two bytes of 256*i+j in little-endian byte order?
            k = 0
            while k < 256:
                b0 = next(XOF)
                b1 = next(XOF)
                b2 = next(XOF)
                cand = b0 + (b1<<8) + ((b2 & 127)<<16)
                if cand < QQ:
                    Ahat[i][j][k] = cand
                    k += 1
            XOF.close()
    return Ahat

def ExpandMask(rhop, kappa, gamma1, l):
    """
    uses SHAKE256, no rejection sampling.
    18 or 20 bits from byte stream for each entry in [0, 2*gamma1] then subtract from gamma1.
    """
    y = VecZero(l)
    if gamma1 == 2**17:
            outlen = 256*18//8
            e = 18
    elif gamma1 == 2**19:
            outlen = 256*20//8
            e = 20
    mask = (1<<e) - 1
    for i in range(l):
        B = FIPS.SHAKE256(rhop + (kappa+i).to_bytes(2,'little'), outlen)
        sampleint = int.from_bytes(B, 'little')
        for j in range(256):
            y[i][j] = gamma1 - (mask & (sampleint >> (e*j))) # gamma1 - bits_to_int(Bbits[j*e : (j+1)*e])
    return y

def ExpandS(rhop, eta, k, l):
    """
    uses SHAKE256m section 5.3 pg. 20.
    """
    s1s2 = [RqZero() for i in range(k + l)]
    for i in range(k + l):
        XOF = FIPS.gen_SHAKE256(rhop + i.to_bytes(2,'little')) # ??? is two bytes rep i in little-endian byte order?
        j = 0
        # rejection sampling for coefficients
        while j < 256:
            b = next(XOF)
            for cand in [b&15, b>>4]:
                if j <256:
                    if eta == 2:
                        if cand < 15:
                            cand = cand % 5
                            cand = eta - cand
                            s1s2[i][j] = cand
                            j += 1
                    elif eta == 4:
                        if cand < 9:
                            cand = eta - cand
                            s1s2[i][j] = cand
                            j += 1
    return (s1s2[:l], s1s2[l:l + k])

#------------------------------#
# KeyGen, Sign, Verify
# NOTE: Lots to fix!
#------------------------------#

def KeyGen(zeta = None, level = 5):
    """
    generate public and private key
    """
    params = param_dict(level)
    k = params['KK']
    l = params['LL']
    eta = params['ETA']
    d = params['DD']
    if zeta is None:
        zeta = os.urandom(32) # 256 bits
    rho_rhop_K = H_hash(zeta, 32+64+32)
    #print(ba2h(rho_rhop_K))
    rho = rho_rhop_K[:32]
    rhop = rho_rhop_K[32:32+64]
    K = rho_rhop_K[32+64:]

    Ahat = ExpandA(rho, k, l)
    #print(Ahat[3][3])
    s1, s2 = ExpandS(rhop, eta, k, l)
    #print(s1[0])
    #print(s2[0])
    t = VecAdd(VecINTT(MatVecMult(Ahat, VecNTT(s1), k, l)), s2, k) # A*s1 + s2
    t1, t0 = VecPow2Round(t, d, k)

    #print(t[3])
    #print(t1[3])
    #print(t0[3])
    #print(ba2h(Pack_t1(t1)))
    tr = H_hash(rho + Pack_t1(t1), 32)
    #print("\n", ba2h(tr))
    pk = rho + Pack_t1(t1) #(rho, t1)
    #print(ba2h(pk))
    sk = rho + K + tr + Pack_s(s1, eta) + Pack_s(s2, eta) + Pack_t0(t0) # (rho, K, tr, s1, s2, t0)
    #print("keygen sizes: ", len(rho), len(K), len(tr), len(Pack_s(s1, eta)), len(Pack_s(s2, eta)), len(Pack_t0(t0)))
    # print("\tKeyGen packing test:")
    # pups1 = Unpack_s(Pack_s(s1, eta),l,eta)
    # pups2 = Unpack_s(Pack_s(s2, eta),k,eta)
    # pupt0 = Unpack_t0(Pack_t0(t0),k)
    # pupt1 = Unpack_t1(Pack_t1(t1),k)
    # print("s1: ", pups1 == s1)
    # print("s2: ", pups2 == s2)
    # print("t0: ", pupt0 == t0)
    # print("t1: ", pupt1 == t1)
    #print("t0: ", t0)
    #print("pupt0: ", pupt0)
    # print("s1[0]: ", s1[0])
    # print("s2[0]: ", s2[0])
    # print("t0[0]: ", t0[0])
    # print("t1[0]: ", t1[0])
    # print("rho: ", ba2h(rho))
    return (pk, sk)

def Sign(sk, M, randomized = False, coins = None, level = 5):
    """
    sign M with sk

    estimated prob of rejection in first check: 0.22, ???, and ??? for levels 2,3,5.
    """
    params = param_dict(level)
    k = params['KK']
    l = params['LL']
    gamma1 = params['GAMMA1']
    gamma2 = params['GAMMA2']
    beta = params['BETA']
    omega = params['OMEGA']
    tau = params['TAU']
    eta = params['ETA']

    rho = sk[:32]
    #print(ba2h(rho))
    Ahat = ExpandA(rho, k, l)
    #print(Ahat[0][0])
    tr = sk[64:96]
    mu = H_hash(tr+M, 64)
    #print("mu", ba2h(mu))
    kappa = 0
    if randomized == False:
        K = sk[32:64]
        rhop = H_hash(K+mu, 64)
    else:
        if coins is None:
            rhop = os.urandom(64)
        else:
            rhop = coins
    if level in [2,5]:
        ls = 96
    elif level == 3:
        ls = 128
    #print("rhop", ba2h(rhop))
    s1 = Unpack_s(sk[96:96+ls*l], l, eta)
    s2 = Unpack_s(sk[96+ls*l:96+ls*l+ls*k], k, eta)
    t0 = Unpack_t0(sk[96+ls*l+ls*k:], k)
    s1hat = VecNTT(s1)
    s2hat = VecNTT(s2)
    t0hat = VecNTT(t0)

    success = False
    #reject = 0
    while not success:
        y = ExpandMask(rhop, kappa, gamma1, l)
        #print("y[0]", y[0])
        #print("OK")
        #print(Ahat)
        #print("y = ", y)
        w = VecINTT(MatVecMult(Ahat, VecNTT(y), k, l))
        w1 = VecHighBits(w, 2*gamma2, k)
        #print("w[0]", w[0])
        #print("w1[0]", w1[0])
        #print("w1", ba2h(Pack_w1(w1, level)))
        ctilde = H_hash(mu + Pack_w1(w1, level), 32)
        #print("ctilde", ba2h(ctilde))
        c = SampleInBall(ctilde, tau)
        #print(c)
        chat = NTT(c)
        z = VecAdd(y, VecINTT(VecScale(chat, s1hat)), l)
        r0 = VecLowBits(
                    VecSub(w, VecINTT(VecScale(chat, s2hat)), k),
                    2*gamma2, k
                    ) # w - c*s2 # NTT!
        if (VecInfNorm(z) >= gamma1 - beta) or (VecInfNorm(r0) >= gamma2 - beta):
            success = False
            # print("bounds (z,r0): ", gamma1 - beta, gamma2 - beta)
            # #print("|c*s1| = ", VecInfNorm(VecINTT(VecScale(chat, s1hat))))
            # #print("|y| = ", VecInfNorm(y))
            # print("|z| = ", VecInfNorm(z))
            # print("|r0| = ", VecInfNorm(r0))
            #print("\tSign rejection1")
            #reject += 1
        else:
            # h = makehint(-c*t0, w-c*s2+c*t0, 2*gamma2)
            arg1 = VecSub(VecZero(k), VecINTT(VecScale(chat, t0hat)), k)
            arg2 = VecSub(w, VecINTT(VecScale(chat, VecSub(s2hat, t0hat, k))), k)
            h = VecMakeHint(arg1, arg2, 2*gamma2, k)
            if (VecInfNorm(VecINTT(VecScale(chat, t0hat))) >= gamma2) or (VecCount1s(h, k) > omega):
                success = False
                #print("arg1 = ", arg1)
                #print("arg2 = ", arg2)
                # print("bounds (c*t0, ones): ", gamma2, omega)
                # #print("h = ", h)
                # print("|c*t0| = ", VecInfNorm(VecINTT(VecScale(chat, t0hat))))
                # print("ones: ", VecCount1s(h, k))
                #print("\tSign rejection2")
                #reject += 1
            else:
                success = True
                #print("number of rejections: ", reject)
        kappa += l
    sigma = ctilde + Pack_z(z, gamma1) + Pack_h(h, omega) # (ctilde, z, h)

    #print("\tSign packing test:")
    #pupz = Unpack_z(Pack_z(z,gamma1),l,gamma1)
    #print("z: ", [[d[i] % QQ for i in range(256)] for d in pupz] == z)
    #print("h: ", Unpack_h(Pack_h(h,omega),k,omega) == h)
    #print("z: ", z)
    #print("pupz: ", pupz)
    # print("s1[0]: ", s1[0])
    # print("s2[0]: ", s2[0])
    # print("t0[0]: ", t0[0])
    # print("t1[0]: ", t1[0])
    # print("z[0]: ", z[0])
    # print("h[0]: ", h[0])
    # print("rho: ", ba2h(rho))
    # print("sign w1[0][:10]: ", w1[0][:10])
    # print("sign ctilde: ", ba2h(ctilde))
    return sigma

def Verify(pk, M, sigma, level = 5):
    """
    verify sigma on M with pk
    """
    params = param_dict(level)
    k = params['KK']
    l = params['LL']
    gamma1 = params['GAMMA1']
    gamma2 = params['GAMMA2']
    beta = params['BETA']
    omega = params['OMEGA']
    tau = params['TAU']
    d = params['DD']

    rho =  pk[:32]
    t1 = Unpack_t1(pk[32:], k)
    Ahat = ExpandA(rho, k, l)
    mu = H_hash(H_hash(rho + Pack_t1(t1), 32) + M, 64)
    ctilde = sigma[:32]
    if level in [3,5]:
        lz = 640
    elif level == 2:
        lz = 576
    z = Unpack_z(sigma[32:32+lz*l], l, gamma1)
    h = Unpack_h(sigma[32+lz*l:], k, omega)
    c = SampleInBall(ctilde, tau)
    chat = NTT(c)
    # compute three lines below with NTT
    scalarhat = PWMult(chat, NTT(RqScale(2**d, RqOne()))) # 2^d*c in NTT domain
    zhat = VecNTT(z)
    xhat = VecSub(MatVecMult(Ahat, zhat, k, l), VecScale(scalarhat, VecNTT(t1)), k) # A*z-(2^d*c)*t1 in NTT domain
    w1prime = VecUseHint(h, VecINTT(xhat), 2*gamma2, k)

    # print("t1[0]: ", t1[0])
    # print("z[0]: ", z[0])
    # print("h[0]: ", h[0])
    # print("rho: ", ba2h(rho))
    # print("verify w1prime[0][:10]: ", w1prime[0][:10])
    # print("verify ctilde: ", ba2h(ctilde))
    # print("z check: ", VecInfNorm(z), " <? ", gamma1 - beta)
    # print("h check: ", VecCount1s(h, k), " <? ", omega)
    # print("ctilde check: ", ba2h(ctilde)," ==? ", ba2h(H_hash(mu + Pack_w1(w1prime, level), 32)))

    return (VecInfNorm(z) < gamma1 - beta) and (ctilde == H_hash(mu + Pack_w1(w1prime, level), 32)) and (VecCount1s(h, k) <= omega)
