#--------------------------------------------------------------------------------------------------#
# DILITHIUM VERSION 3.1
#
# SPECS: https://pq-crystals.org/dilithium/resources.shtml
#
# Run 300 keygen/sign and compare to KAT from reference implementation
# https://github.com/pq-crystals/dilithium
# 100 for each security level 2,3,5
# Only uses the deterministic signing.
#
# Another python implementation:
# https://github.com/jack4818/dilithium-py
#--------------------------------------------------------------------------------------------------#

import dilithium
import aes256_ctr_drbg # taken from https://github.com/jack4818/kyber-py
from test import ba2h

def expand_KAT_seed(seed):
    """
    expand the 48 byte seed from NIST submission KAT.
    uses AES256 in counter mode to generate all necessary randomness from seed.

    """
    X = aes256_ctr_drbg.AES256_CTR_DRBG(seed=seed)
    zeta = X.random_bytes(32) # used in KeyGen
    return zeta

with open("kat-files/PQCsignKAT_Dilithium2.rsp", 'r') as KATfile:
    level = 2
    print("level 2 dilithium-d KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                zeta = expand_KAT_seed(seed)
                mypk, mysk = dilithium.KeyGen(zeta = zeta, level = level)
            elif words[0] == "mlen":
                mlen = words[2]
            elif words[0] == "msg":
                msg = words[2].lower()
                mysig = dilithium.Sign(mysk, bytes.fromhex(msg), randomized = False, coins = None, level = level)
            elif words[0] == "pk":
                pk = words[2].lower()
                if pk == mypk.hex():
                    print("\tpk success")
                else:
                    print("\tpk failure")
                    #print(pk[:192])
                    #print(mypk.hex()[:192])
                    #print("/t", mypk.hex()), len(pk)//2)
                    #print(len(mypk.hex())//2)
                    #print([int(pk[i] == mypk.hex()[i]) for i in range(len(pk))])
            elif words[0] == "sk":
                sk = words[2].lower()
                if sk == mysk.hex():
                    print("\tsk success")
                else:
                    print("\tsk failure")
                    #print(sk[:192])
                    #print(mysk.hex()[:192])
                    #print("\t", len(sk)//2, len(mysk.hex())//2)
                    #print([int(sk[i] == mysk.hex()[i]) for i in range(len(sk))])
            elif words[0] == "smlen":
                smlen = words[2]
                mysmlen = len(mysig) + int(mlen)
                if int(smlen) == mysmlen:
                    print("\tsmlen success")
                else:
                    print("\tsmlen failure")
            elif words[0] == "sm":
                sm = words[2].lower()
                mysm = mysig + bytes.fromhex(msg)
                if sm == mysm.hex():
                    print("\tsm success")
                else:
                    print("\tsm failure")
                    #print(sm[:192])
                    #print(mysm.hex()[:192])
                    #print(len(sm))#, sm)
                    #print(len(mysm.hex()))#, mysm.hex())
                    #print([int(sm[i] == mysm.hex()[i]) for i in range(len(sm))])
            else:
                pass

with open("kat-files/PQCsignKAT_Dilithium3.rsp", 'r') as KATfile:
    level = 3
    print("level 3 dilithium-d KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                zeta = expand_KAT_seed(seed)
                mypk, mysk = dilithium.KeyGen(zeta = zeta, level = level)
            elif words[0] == "mlen":
                mlen = words[2]
            elif words[0] == "msg":
                msg = words[2].lower()
                mysig = dilithium.Sign(mysk, bytes.fromhex(msg), randomized = False, coins = None, level = level)
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
                    #print("\t", len(sk)//2, len(mysk.hex())//2)
            elif words[0] == "smlen":
                smlen = words[2]
                mysmlen = len(mysig) + int(mlen)
                if int(smlen) == mysmlen:
                    print("\tsmlen success")
                else:
                    print("\tsmlen failure")
            elif words[0] == "sm":
                sm = words[2].lower()
                mysm = mysig + bytes.fromhex(msg)
                if sm == mysm.hex():
                    print("\tsm success")
                else:
                    print("\tsm failure")
            else:
                pass

with open("kat-files/PQCsignKAT_Dilithium5.rsp", 'r') as KATfile:
    level = 5
    print("level 5 dilithium-d KAT")
    for line in KATfile:
        if line != "\n":
            words = line.split()
            if words[0] == "count":
                count = words[2]
                print(f"count {count}")
            elif words[0] == "seed":
                hexseed = words[2]
                seed = bytes.fromhex(hexseed)
                zeta = expand_KAT_seed(seed)
                mypk, mysk = dilithium.KeyGen(zeta = zeta, level = level)
            elif words[0] == "mlen":
                mlen = words[2]
            elif words[0] == "msg":
                msg = words[2].lower()
                mysig = dilithium.Sign(mysk, bytes.fromhex(msg), randomized = False, coins = None, level = level)
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
                    #print("\t", len(sk)//2, len(mysk.hex())//2)
            elif words[0] == "smlen":
                smlen = words[2]
                mysmlen = len(mysig) + int(mlen)
                if int(smlen) == mysmlen:
                    print("\tsmlen success")
                else:
                    print("\tsmlen failure")
            elif words[0] == "sm":
                sm = words[2].lower()
                mysm = mysig + bytes.fromhex(msg)
                if sm == mysm.hex():
                    print("\tsm success")
                else:
                    print("\tsm failure")
            else:
                pass
