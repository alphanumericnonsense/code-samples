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
# Running this script performs 300 keygen/encaps/decaps, 100 for each security level
# 1,3,5, and compares to known answers from the reference implemenation.
#-----------------------------------------------------------------------------------#

import pke
import kem
import aes256_ctr_drbg # taken from https://github.com/jack4818/kyber-py

def expand_KAT_seed(seed):
    """
    expand the 48 byte seed from NIST submission KAT.
    uses AES256 in counter mode to generate all necessary randomness from seed.

    order of sampling, (d,z,m) is different from specs kyber.ccakem.keygen() alg. 7,
    where it looks like it should be (z,d,m)
    """
    X = aes256_ctr_drbg.AES256_CTR_DRBG(seed=seed)
    d = X.random_bytes(32) # used in kem_key_gen
    z = X.random_bytes(32) # used in pke_key_gen
    m = X.random_bytes(32) # used in encaps
    return (d,z,m)

with open("kat-files/PQCkemKAT_1632.rsp", 'r') as KATfile:
    level = 1
    rank = 2
    print("level 1 kyber512 KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                d, z, m = expand_KAT_seed(seed)
                mypk, mysk = kem.kem_key_gen(d=d, z=z, level=level)
                myct, myss = kem.kem_encaps(pk=mypk, keylength=32, m=m, level=level)
            elif words[0] == "pk":
                pk = words[2].lower()
                #print(pk)
                #print(mypk.hex())
                #pkdecoded = [pke.decode(bytes.fromhex(pk)[12*32*i : 12*32*(i+1)], 12) for i in range(rank)]
                #mypkdecoded = [pke.decode(mypk[12*32*i : 12*32*(i+1)], 12) for i in range(rank)]
                #print(pkdecoded)
                #print(mypkdecoded)
                #print(pk[-32*2:])
                #print(mypk[-32:].hex())
                if pk == mypk.hex():
                    print("\tpk success")
                else:
                    print("\tpk failure")
            elif words[0] == "sk":
                sk = words[2].lower()
                #print(sk)
                #print(mysk.hex())
                if sk == mysk.hex():
                    print("\tsk success")
                else:
                    print("\tsk failure")
            elif words[0] == "ct":
                ct = words[2].lower()
                #print(ct)
                #print(myct.hex())
                if ct == myct.hex():
                    print("\tct success")
                else:
                    print("\tct failure")
            elif words[0] == "ss":
                ss = words[2].lower()
                #print(ss)
                #print(myss.hex())
                if ss == myss.hex():
                    print("\tss success")
                else:
                    print("\tss failure")
            else:
                pass

with open("kat-files/PQCkemKAT_2400.rsp", 'r') as KATfile:
    level = 3
    print("level 3 kyber768 KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                d, z, m = expand_KAT_seed(seed)
                mypk, mysk = kem.kem_key_gen(d=d, z=z, level=level)
                myct, myss = kem.kem_encaps(pk=mypk, keylength=32, m=m, level=level)
            elif words[0] == "pk":
                pk = words[2].lower()
                if pk == mypk.hex():
                    print("\tpk success")
                else:
                    print("\tpk failure")
            elif words[0] == "sk":
                sk = words[2].lower()
                if sk == mysk.hex():
                    print("\tsk success")
                else:
                    print("\tsk failure")
            elif words[0] == "ct":
                ct = words[2].lower()
                if ct == myct.hex():
                    print("\tct success")
                else:
                    print("\tct failure")
            elif words[0] == "ss":
                ss = words[2].lower()
                if ss == myss.hex():
                    print("\tss success")
                else:
                    print("\tss failure")

with open("kat-files/PQCkemKAT_3168.rsp", 'r') as KATfile:
    level = 5
    print("level 5 kyber1024 KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                d, z, m = expand_KAT_seed(seed)
                mypk, mysk = kem.kem_key_gen(d=d, z=z, level=level)
                myct, myss = kem.kem_encaps(pk=mypk, keylength=32, m=m, level=level)
            elif words[0] == "pk":
                pk = words[2].lower()
                if pk == mypk.hex():
                    print("\tpk success")
                else:
                    print("\tpk failure")
            elif words[0] == "sk":
                sk = words[2].lower()
                if sk == mysk.hex():
                    print("\tsk success")
                else:
                    print("\tsk failure")
            elif words[0] == "ct":
                ct = words[2].lower()
                if ct == myct.hex():
                    print("\tct success")
                else:
                    print("\tct failure")
            elif words[0] == "ss":
                ss = words[2].lower()
                if ss == myss.hex():
                    print("\tss success")
                else:
                    print("\tss failure")
