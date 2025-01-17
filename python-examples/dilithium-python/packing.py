#--------------------------------------------------------------------------------------------------#
# DILITHIUM VERSION 3.1
#
# SPECS: https://pq-crystals.org/dilithium/resources.shtml
#
# Reference implementation
# https://github.com/pq-crystals/dilithium
#
# All the annoying bitpacking...
#--------------------------------------------------------------------------------------------------#

import math

QQ = 8380417

def VecZero(dim):
    return [[0 for i in range(256)] for j in range(dim)]

#------------------------------#
# bit/byte packing and inverse
#------------------------------#

def PackVec(v,nbits,nbytes):
    r = 0
    dim = len(v)
    for vi in reversed(v):
        for vij in reversed(vi):
            r = (r<<nbits) | vij
    #print(math.log2(r), dim, nbytes)
    return r.to_bytes(dim*nbytes, 'little')

def UnpackVec(B, dim, nbits):
    r = int.from_bytes(B, 'little')
    v = [[0 for j in range(256)] for i in range(dim)]
    mask = (1<<nbits) - 1
    for i in range(dim):
        for j in range(256):
            v[dim-1-i][256-1-j] = mask & (r>>256*nbits*(dim-1-i) + (256-1-j)*nbits)
    return v

def bytes_to_bits(B):
    """
    take byte stream to list of bits
    B = (b0b1b2b3b4b5b6b7)...
    """
    retbits = []
    for b in B:
        bbits = []
        for i in range(8):
            bbits.append(b&1)
            b >>= 1
        bbits.reverse()
        retbits += bbits
    return retbits

def bits_to_int(b):
    """
    list of bits to int
    """
    brev = reversed(b) # not sure if needed, specs not clear?
    S = 0
    for i, d in enumerate(brev):
        S += d*(2**i)
    return S

def int_to_bits(n, width):
    """
    return list of 0/1 of length width for binary rep of n, e.g.
    int_to_bit(6,5) = [0,0,1,1,0]
    """
    x = bin(n)[2:] # string
    y = "0"*(width-len(x)) + x # string
    return [int(i) for i in y]

def bits_to_bytes(L):
    """
    list of bits to bytearray, assume 8|len(L).
    """
    B = b''
    for i in range(len(L)//8):
        nextint = bits_to_int(L[8*i:8*(i+1)])
        B += nextint.to_bytes(1,'big')
    return B

def Pack_w1(w1, level):
    if level in [3,5]:
        #print(w1)
        return PackVec(w1, 4, 256*4//8)
    elif level == 2:
        return PackVec(w1, 6, 256*6//8)

def Unpack_w1(w1, k, level):
    if level == 2:
        return UnpackVec(w1, k, 6)
    elif level in [3,5]:
        return UnpackVec(w1, k, 4)

# def Pack_w1(w1, level = 5):
#     """
#     pg. 17, pg. 18 fig 5, fig 6.
#     """
#     B = b''
#     bits = []
#     if level == 2:
#         # coeffs in [0,15], 4 bits
#         # each byte is 2 coefficients
#         for w1i in w1:
#             for w1ij in w1i:
#                 x = int_to_bits(w1ij,4)
#                 bits += x
#                 if len(bits) == 8:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     elif level in [3,5]:
#         # coeffs in [0, 43], 6 bits
#         # every 3 bytes is 4 coefficients
#         for w1i in w1:
#             for w1ij in w1i:
#                 x = int_to_bits(w1ij,6)
#                 bits += x
#                 if len(bits) == 24:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     return B
#
# def Unpack_w1(w1, k, level):
#     up_w1 = VecZero(k)
#     if level == 2: # 4 bits
#         for i in range(k):
#             for j in range(128):
#                 # each byte is two coeffs
#                 x = w1[128*i+j]
#                 y1y2 = bytes_to_bits(bytearray([x]))
#                 y1 = bits_to_int(y1y2[:4])
#                 y2 = bits_to_int(y1y2[4:])
#                 up_w1[i][2*j] = y1
#                 up_w1[i][2*j+1] = y2
#     elif level in [3,5]: # 6 bits
#         for i in range(k):
#             for j in range(64):
#                 # every 3 byte is four coeffs
#                 x1x2x3 = w1[3*64*i+3*j:3*64*i+3*(j+1)]
#                 y1_to_y4 = bytes_to_bits(x1x2x3)
#                 y1 = bits_to_int(y1_to_y4[:6])
#                 y2 = bits_to_int(y1_to_y4[6:12])
#                 y3 = bits_to_int(y1_to_y4[12:18])
#                 y4 = bits_to_int(y1_to_y4[18:24])
#                 up_w1[i][4*j] = y1
#                 up_w1[i][4*j+1] = y2
#                 up_w1[i][4*j+2] = y3
#                 up_w1[i][4*j+3] = y4
#     return up_w1

def Pack_t1(t1):
    return PackVec(t1, 10, 256*10//8)

def Unpack_t1(t1, k):
    return UnpackVec(t1, k, 10)

def Pack_t0(t0):
    t0alt = [[2**12 - t0ij for t0ij in t0i] for t0i in t0]
    return PackVec(t0alt, 13, 256*13//8)

def Unpack_t0(t0, k):
    x = UnpackVec(t0, k, 13)
    return [[2**12 - xij for xij in xi] for xi in x]

# def Pack_t0(t0):
#     """
#     pg. 18, pg. 19 fig 8.
#     coefficients are 2^12-v, v in [0,2^13).
#     stores these v in 13 bits.
#     8 coeffs in 13 bytes
#     """
#     B = b''
#     bits = []
#     for t0i in t0:
#         for t0ij in t0i:
#             x = int_to_bits(2**12 - t0ij, 13)
#             bits += x
#             if len(bits) == 104:
#                 B += bits_to_bytes(bits)
#                 bits = []
#     return B
#
# def Unpack_t0(t0, k):
#     up_t0 = VecZero(k)
#     for i in range(k):
#         for j in range(32):
#             # every 13 bytes is 8 coeffs
#             x1_to_x13 = t0[13*32*i+13*j:13*32*i+13*(j+1)]
#             y1_to_y8 = bytes_to_bits(x1_to_x13)
#             y1 = bits_to_int(y1_to_y8[:13])
#             y2 = bits_to_int(y1_to_y8[13:26])
#             y3 = bits_to_int(y1_to_y8[26:39])
#             y4 = bits_to_int(y1_to_y8[39:52])
#             y5 = bits_to_int(y1_to_y8[52:65])
#             y6 = bits_to_int(y1_to_y8[65:78])
#             y7 = bits_to_int(y1_to_y8[78:91])
#             y8 = bits_to_int(y1_to_y8[91:])
#             up_t0[i][8*j] = 2**12 - y1
#             up_t0[i][8*j+1] = 2**12 - y2
#             up_t0[i][8*j+2] = 2**12 - y3
#             up_t0[i][8*j+3] = 2**12 - y4
#             up_t0[i][8*j+4] = 2**12 - y5
#             up_t0[i][8*j+5] = 2**12 - y6
#             up_t0[i][8*j+6] = 2**12 - y7
#             up_t0[i][8*j+7] = 2**12 - y8
#     return up_t0

def Pack_s(s, eta):
    salt = [[(eta-sij)%QQ for sij in si] for si in s]
    if eta == 2:
        return PackVec(salt, 3, 256*3//8)
    elif eta == 4:
        return PackVec(salt, 4, 256*4//8)

def Unpack_s(s, dim, eta):
    if eta == 2:
        x = UnpackVec(s, dim, 3)
    elif eta == 4:
        x = UnpackVec(s, dim, 4)
    return [[(eta-xij)%QQ for xij in xi] for xi in x]

# def Pack_s(s, eta):
#     """
#     pg. 18, pg. 19 fig 9.
#     each coeff is eta-c mod q for some c in [0,2*eta].
#     store these c in either 3 or 4 bits per coeff depending on eta.
#     """
#     #QQ = 8380417
#     B = b''
#     bits = []
#     if eta == 2: # levels 2,5
#         # 3 bits per coeff, 8 coeffs per 3 bytes
#         for si in s:
#             for sij in si:
#                 x = int_to_bits((eta - sij) % QQ, 3)
#                 bits += x
#                 if len(bits) == 24:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     elif eta == 4: # level 3
#         # 4 bits per coeff, 2 coeff per byte
#         for si in s:
#             for sij in si:
#                 x = int_to_bits((eta - sij) % QQ, 4)
#                 bits += x
#                 if len(bits) == 8:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     #print("pack_s", len(B), len(s))
#     return B
#
# def Unpack_s(s, dim, eta):
#     """
#     bytes to bits to ints (eta - c)
#     """
#     #print("unpack_s", len(s), dim)
#     up_s = VecZero(dim)
#     if eta == 2: # 3 bits, level 2,5
#         for i in range(dim):
#             for j in range(32):
#                 # every 3 bytes is 8 coeffs
#                 x1x2x3 = s[3*32*i+3*j:3*32*i+3*(j+1)]
#                 y1_to_y8 = bytes_to_bits(x1x2x3)
#                 y1 = bits_to_int(y1_to_y8[:3])
#                 y2 = bits_to_int(y1_to_y8[3:6])
#                 y3 = bits_to_int(y1_to_y8[6:9])
#                 y4 = bits_to_int(y1_to_y8[9:12])
#                 y5 = bits_to_int(y1_to_y8[12:15])
#                 y6 = bits_to_int(y1_to_y8[15:18])
#                 y7 = bits_to_int(y1_to_y8[18:21])
#                 y8 = bits_to_int(y1_to_y8[21:24])
#                 up_s[i][8*j] = eta - y1
#                 up_s[i][8*j+1] = eta - y2
#                 up_s[i][8*j+2] = eta - y3
#                 up_s[i][8*j+3] = eta - y4
#                 up_s[i][8*j+4] = eta - y5
#                 up_s[i][8*j+5] = eta - y6
#                 up_s[i][8*j+6] = eta - y7
#                 up_s[i][8*j+7] = eta - y8
#     elif eta == 4: # 4 bits, level 3
#         for i in range(dim):
#             for j in range(128):
#                 # every byte is 2 coeffs
#                 x = s[128*i+j]
#                 #print(x)
#                 y1y2 = bytes_to_bits(bytearray([x]))
#                 y1 = bits_to_int(y1y2[:4])
#                 y2 = bits_to_int(y1y2[4:])
#                 up_s[i][2*j] = eta - y1
#                 up_s[i][2*j+1] = eta - y2
#     return up_s


def Pack_z(z,gamma1):
    zalt = [[(gamma1-zij)%QQ for zij in zi] for zi in z]
    if gamma1 == 2**17:
        return PackVec(zalt, 18, 256*18//8)
    elif gamma1 == 2**19:
        return PackVec(zalt, 20, 256*20//8)

def Unpack_z(z, l, gamma1):
    if gamma1 == 2**17:
        x = UnpackVec(z, l, 18)
    elif gamma1 == 2**19:
        x = UnpackVec(z, l, 20)
    return [[(gamma1-xij)%QQ for xij in xi] for xi in x]

# def Pack_z(z, gamma1):
#     """
#     pg. 18, pg. 20 fig 10.
#     coeffs are gamma1-c mod q with c in [0, 2*gamma1).
#     """
#     #QQ = 8380417
#     B = b''
#     bits = []
#     if gamma1 == 2**17: # level 2
#         # 18 bits per coeff, 4 coeffs in 9 bytes.
#         for zi in z:
#             for zij in zi:
#                 x = int_to_bits((gamma1 - zij) % QQ, 18)
#                 bits += x
#                 if len(bits) == 72:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     elif gamma1 == 2**19: # level 3,5
#         # 20 bits per coeff, 2 coeffs in 5 bytes.
#         for zi in z:
#             for zij in zi:
#                 x = int_to_bits((gamma1 - zij) % QQ, 20)
#                 bits += x
#                 if len(bits) == 40:
#                     B += bits_to_bytes(bits)
#                     bits = []
#     return B
#
# def Unpack_z(z, l, gamma1):
#     up_z = VecZero(l)
#     if gamma1 == 2**17: # 18 bits
#         for i in range(l):
#             for j in range(64):
#                 # every 9 bytes is 4 coeffs
#                 x1_to_x9 = z[9*64*i+9*j:9*64*i+9*(j+1)]
#                 y1_to_y4 = bytes_to_bits(x1_to_x9)
#                 y1 = bits_to_int(y1_to_y4[:18])
#                 y2 = bits_to_int(y1_to_y4[18:36])
#                 y3 = bits_to_int(y1_to_y4[36:54])
#                 y4 = bits_to_int(y1_to_y4[54:72])
#                 up_z[i][4*j] = gamma1 - y1
#                 up_z[i][4*j+1] = gamma1 - y2
#                 up_z[i][4*j+2] = gamma1 - y3
#                 up_z[i][4*j+3] = gamma1 - y4
#     elif gamma1 == 2**19: # 20 bits
#         for i in range(l):
#             for j in range(128):
#                 # every 5 bytes is 2 coeffs
#                 x1_to_x5 = z[5*128*i+5*j:5*128*i+5*(j+1)]
#                 y1y2 = bytes_to_bits(x1_to_x5)
#                 y1 = bits_to_int(y1y2[:20])
#                 y2 = bits_to_int(y1y2[20:40])
#                 up_z[i][2*j] = gamma1 - y1
#                 up_z[i][2*j+1] = gamma1 - y2
#     return up_z

def Pack_h(h, omega):
    """
    pg. 21.
    store location (or zero) of at most omega nonzero coeffs, 1 byte each.
    bytes [omega, omega + k - 1] store poly boundaries j, 0 <= j <= omega, one byte each?
    """
    nzind = b'' # non-zero coeff indices
    boundaries = b'' # poly boundaries
    countbound = 0
    for hi in h:
        for i, hij in enumerate(hi):
            if hij == 1:
                nzind += i.to_bytes(1,'big')
                countbound += 1
        boundaries += countbound.to_bytes(1,'big')
    k = len(h)
    for i in range(omega - len(nzind)):
        # fill with zeros
        nzind += (0).to_bytes(1,'big')
    return nzind + boundaries

def Unpack_h(h, k, omega):
    up_h = VecZero(k)
    nzind = h[:omega]
    boundaries = h[omega:]
    for i in range(k):
        if i > 0:
            start = boundaries[i-1]
        else:
            start = 0
        for j in range(start,boundaries[i]):
            up_h[i][nzind[j]] = 1
    return up_h

def ba2h(list_val):
     result = ''.join('{:02x}'.format(x) for x in list_val)
     return(result)
