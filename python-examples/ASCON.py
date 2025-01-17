#-------------------------------------------------------#
# ASCON-128, ASCON-HASH, ASCON-XOF
#
# internal state is list of 5 64-bit integers
# internal key, nonce, tag are 128-bit integers
# rate for all is 64-bits
# bytes in, bytes out
# XOF output length specified in bytes
#
# major functions:
#
# ascon_encrypt(key, nonce, ad, pt) return (ct, tag)
# ascon_decrypt(key, nonce, ad, ct, tag) return pt or None
# ascon_hash(msg) return digest
# ascon_xof(iv, outlen) return digest
#-------------------------------------------------------#

def print_state(state, fmt="bin"):
    if fmt == "bin":
        print(f"x0:{format(state[0], '064b')}\nx1:{format(state[1], '064b')}\nx2:{format(state[2], '064b')}\nx3:{format(state[3], '064b')}\nx4:{format(state[4], '064b')}")
    elif fmt == "hex":
        print(f"x0:{format(state[0], '016x')}\nx1:{format(state[1], '016x')}\nx2:{format(state[2], '016x')}\nx3:{format(state[3], '016x')}\nx4:{format(state[4], '016x')}")
    elif fmt == "int":
        print(f"x0:{state[0]}\nx1:{state[1]}\nx2:{state[2]}\nx3:{state[3]}\nx4:{state[4]}")

debug = False

#----------------------------------#
# permutation
#----------------------------------#

# number of rounds
a_aead = 12
b_aead = 6

a_hash = 12
b_hash = 12

a_xof = 12
b_xof = 12

RC = (0x00000000000000f0,
      0x00000000000000e1,
      0x00000000000000d2,
      0x00000000000000c3,
      0x00000000000000b4,
      0x00000000000000a5,
      0x0000000000000096,
      0x0000000000000087,
      0x0000000000000078,
      0x0000000000000069,
      0x000000000000005a,
      0x000000000000004b)

SBOX = (0x4, 0xb, 0x1f, 0x14, 0x1a, 0x15, 0x9, 0x2,
        0x1b, 0x5, 0x8, 0x12, 0x1d, 0x3, 0x6, 0x1c,
        0x1e, 0x13, 0x7, 0xe, 0x0, 0xd, 0x11, 0x18,
        0x10, 0xc, 0x1, 0x19, 0x16, 0xa, 0xf, 0x17)

def get_slice(state, i):
    """
    right to left, i=0 is LSB, i=63 is MSB
    """
    return (((state[0] >> (63 - i)) & 1) << 4) ^ \
           (((state[1] >> (63 - i)) & 1) << 3) ^ \
           (((state[2] >> (63 - i)) & 1) << 2) ^ \
           (((state[3] >> (63 - i)) & 1) << 1) ^ \
           (((state[4] >> (63 - i)) & 1) << 0)

def right_circular_shift_64(n, s):
    return ((n & ((1 << s) - 1)) << (64 - s)) ^ (n >> s)

def linear0(x0):
    return x0 ^ right_circular_shift_64(x0, 19) ^ right_circular_shift_64(x0, 28)

def linear1(x1):
    return x1 ^ right_circular_shift_64(x1, 61) ^ right_circular_shift_64(x1, 39)

def linear2(x2):
    return x2 ^ right_circular_shift_64(x2, 1) ^ right_circular_shift_64(x2, 6)

def linear3(x3):
    return x3 ^ right_circular_shift_64(x3, 10) ^ right_circular_shift_64(x3, 17)

def linear4(x4):
    return x4 ^ right_circular_shift_64(x4, 7) ^ right_circular_shift_64(x4, 41)

def add_round_constant(state, rnd, num_rnds):
    offset = {6:6, 8:4, 12:0}
    state[2] = state[2] ^ RC[rnd + offset[num_rnds]]
    return state

def substitution(state):
    newstate = [0,0,0,0,0]
    for i in range(64):
        # print(get_slice(state, i))
        subslice = SBOX[get_slice(state, i)]
        for j in range(5):
            newstate[j] ^= ((subslice >> (4-j) ) & 1) << (63 - i)
    return newstate

def linear(state):
    return [linear0(state[0]), linear1(state[1]), linear2(state[2]), linear3(state[3]), linear4(state[4])]

def ascon_perm(state, num_rnds):
    for i in range(num_rnds):
        state = add_round_constant(state, i, num_rnds)
        # print("\tafter rc")
        # print_state(state, fmt="hex")
        state = substitution(state)
        # print("\tafter sub")
        # print_state(state, fmt="hex")
        state = linear(state)
        # print("\tafter lin")
        # print_state(state, fmt="hex")
    return state

#----------------------------------#
# initialization
#----------------------------------#
IV_AEAD = 0x80400c0600000000

def initialize(key, nonce):
    """
    key and nonce 16 bytes
    """
    # convert key, nonce bytes to ints
    key0 = int.from_bytes(key[:8], 'big')
    key1 = int.from_bytes(key[8:], 'big')
    nonce0 = int.from_bytes(nonce[:8], 'big')
    nonce1 = int.from_bytes(nonce[8:], 'big')

    state = [IV_AEAD, key0, key1, nonce0, nonce1]
    if debug == True:
        print("start")
        print_state(state, fmt="hex")
    state = ascon_perm(state, a_aead)
    state[3] ^= key0
    state[4] ^= key1
    return state

#----------------------------------#
# finalization
#----------------------------------#
def finalize(state, key):
    """
    output tag as bytes
    """
    # convert key bytes to ints
    key0 = int.from_bytes(key[:8], 'big')
    key1 = int.from_bytes(key[8:], 'big')

    state[1] ^= key0
    state[2] ^= key1
    state = ascon_perm(state, a_aead)
    t0 = state[3] ^ key0
    t1 = state[4] ^ key1
    tag = (t0 << 64) ^ t1
    return tag.to_bytes(16, 'big')

#----------------------------------#
# encryption/decryption
#----------------------------------#

def pad_bytes_64(B):
    """
    append 1 bit, pad to multiple of 64 bits,
    assumes byte-size input
    """
    B += bytes([128])
    numzerobytes = (8 - (len(B) % 8)) % 8
    B += bytes(numzerobytes)
    return B

def absorb_ad_block(state, adblock):
    """
    adblock is 64-bit integer
    """
    state[0] ^= adblock
    state = ascon_perm(state, b_aead)
    return state

def absorb_ad(state, adbytes):
    if len(adbytes) != 0:
        adbytes_padded = pad_bytes_64(adbytes)
    else:
        adbytes_padded = adbytes
    for i in range(len(adbytes_padded) // 8):
        state = absorb_ad_block(state, int.from_bytes(adbytes_padded[8*i : 8*(i+1)], 'big'))
    state[4] ^= 1
    return state

def pt_in_ct_out(state, ptbytes):
    ct = b""
    ptbytes_padded = pad_bytes_64(ptbytes)
    for i in range(len(ptbytes_padded)//8 - 1):
        state[0] ^= int.from_bytes(ptbytes_padded[8*i:8*(i+1)], 'big')
        ct += state[0].to_bytes(8, 'big')
        state = ascon_perm(state, b_aead)
    state[0] ^= int.from_bytes(ptbytes_padded[-8:], 'big')
    lastblocksize = len(ptbytes) % 8
    ct += (state[0] >> (64-8*lastblocksize)).to_bytes(lastblocksize, 'big')
    return (state, ct)

def ct_in_pt_out(state, ctbytes):
    pt = b""
    lastblocksize = len(ctbytes) % 8
    for i in range(len(ctbytes)//8):
        ctblockint = int.from_bytes(ctbytes[8*i:8*(i+1)], 'big')
        ptint = state[0] ^ ctblockint
        pt += ptint.to_bytes(8, 'big')
        state[0] = ctblockint
        state = ascon_perm(state, b_aead)
    if lastblocksize == 0:
        ctlastint = 0
    else:
        ctlastint = int.from_bytes(ctbytes[-lastblocksize:], 'big')
    ptlastint = ctlastint ^ (state[0] >> (64 - 8*lastblocksize))
    pttilde = ptlastint.to_bytes(lastblocksize, 'big')
    pt += pttilde
    state[0] ^= int.from_bytes(pad_bytes_64(pttilde), 'big')
    return (state, pt)

def ascon_encrypt(key, nonce, ad, pt):
    state = initialize(key, nonce)
    if debug == True:
        print("after init")
        print_state(state, fmt="hex")
    state = absorb_ad(state, ad)
    if debug == True:
        print("after ad")
        print_state(state, fmt="hex")
    state, ct = pt_in_ct_out(state, pt)
    if debug == True:
        print("after pt")
        print_state(state, fmt="hex")
    tag = finalize(state, key)
    if debug == True:
        print("end")
        print_state(state, fmt="hex")
    return (ct, tag)

def ascon_decrypt(key, nonce, ad, ct, tag):
    state = initialize(key, nonce)
    if debug == True:
        print("after init")
        print_state(state, fmt="hex")
    state = absorb_ad(state, ad)
    if debug == True:
        print("after ad")
        print_state(state, fmt="hex")
    state, pt = ct_in_pt_out(state, ct)
    if debug == True:
        print("after ct")
        print_state(state, fmt="hex")
    comptag = finalize(state, key)
    if comptag == tag:
        return pt
    else:
        return None

#----------------------------------#
# XOF/hash
#----------------------------------#
INIT_STATE_XOF = (0xb57e273b814cd416, 0x2b51042562ae2420, 0x66a3a7768ddf2218, 0x5aad0a7a8153650c, 0x4f3e0e32539493b6)
INIT_STATE_HASH = (0xee9398aadb67f03d, 0x8bb21831c60f1002, 0xb48a92db98d5da62, 0x43189921b8f8e3e8, 0x348fa5c9d525e140)

def absorb_message(initstate, msgbytes):
    state = initstate
    msgbytes_padded = pad_bytes_64(msgbytes)
    for i in range(len(msgbytes_padded)//8 - 1):
        msgblockint = int.from_bytes(msgbytes_padded[8*i:8*(i+1)], 'big')
        state[0] ^= msgblockint
        state = ascon_perm(state, b_xof)
    msgblockint = int.from_bytes(msgbytes_padded[-8:], 'big')
    state[0] ^= msgblockint
    return state

def squeeze_digest(state, outlen):
    """
    outlen in bytes
    """
    h = b""
    for i in range((outlen + 7)//8):
        state = ascon_perm(state, b_xof)
        h += state[0].to_bytes(8, 'big')
    return h[:outlen]

def ascon_hash(msg):
    state = [INIT_STATE_HASH[i] for i in range(5)]
    state = absorb_message(state, msg)
    digest = squeeze_digest(state, 32)
    return digest

def ascon_xof(iv, outlen):
    state = [INIT_STATE_XOF[i] for i in range(5)]
    state = absorb_message(state, iv)
    digest = squeeze_digest(state, outlen)
    return digest

import ascon_ref
#----------------------------------#
# test AEAD
#----------------------------------#

def enc_dec_test(key, nonce, pt, ad):
    print("key", key.hex())
    print("nonce", nonce.hex())
    print("pt", pt.hex())
    print("ad", ad.hex())

    ct, tag = ascon_encrypt(key, nonce, ad, pt)
    print("ct", ct.hex())
    print("tag", tag.hex())

    ct_tag = ascon_ref.ascon_encrypt(key, nonce, ad, pt, variant="Ascon-128")
    ctref = ct_tag[:-16]
    tagref = ct_tag[-16:]
    print("ctref", ctref.hex())
    print("tagref", tagref.hex())

    pt = ascon_decrypt(key, nonce, ad, ct, tag)
    print("pt recovered", pt.hex())

    if ct == ctref and tag == tagref:
        print("\nPASS\n")
    else:
        print("\nFAIL\n")
    return ct == ctref and tag == tagref

#----------------------------------#
# test hash/XOF
#----------------------------------#

def hash_test(msg):
    digest = ascon_hash(msg)
    print("message", msg.hex())
    print("digest", digest.hex())
    digestref = ascon_ref.ascon_hash(msg, variant="Ascon-Hash", hashlength=32)
    print("digestref", digestref.hex())

    if digest == digestref:
        print("\nPASS\n")
    else:
        print("\nFAIL\n")
    return digest == digestref

def xof_test(iv, outlen):
    digest = ascon_xof(iv, outlen)
    print("iv", iv.hex())
    print("digest", digest.hex())
    digestref = ascon_ref.ascon_hash(iv, variant="Ascon-Xof", hashlength=outlen)
    print("digestref", digestref.hex())
    if digest == digestref:
        print("\nPASS\n")
    else:
        print("\nFAIL\n")
    return digest == digestref

#----------------------------------#
# randomized testing
# against 3rd party
#----------------------------------#

if __name__ == "__main__":
    import os
    print("aead test")
    for i in range(16):
        for j in range(16):
            key = os.urandom(16)
            nonce = os.urandom(16)
            pt = os.urandom(i)
            ad = os.urandom(j)
            x = enc_dec_test(key, nonce, pt, ad)
            assert x

    print("hash test")
    for i in range(16):
        msg = os.urandom(i)
        x = hash_test(msg)
        assert x

    print("xof test")
    for i in range(16):
        for j in range(16):
            iv = os.urandom(i)
            x = xof_test(iv, j)
            assert x
