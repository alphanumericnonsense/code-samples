from Cryptodome.Hash import SHAKE256
import os
import math
import base64

class slhsk:

    def __init__(self, n):
        self.seed = os.urandom(n)
        self.prf = os.urandom(n)

class slhpk:

    def __init__(self, n):
        self.seed = os.urandom(n)
        self.root = None

class slhadrs:
    """
    32-byte addressing
    """

    def __init__(self, other = None):
        self.data = bytes(32) # the ADRS byte string
        if other is not None:
            self.data = other.data

    def setLayerAddres(self, l):
        self.data = l.to_bytes(4, 'big') + self.data[4:32]

    def setTreeAddress(self, t):
        self.data = self.data[0:4] + t.to_bytes(12, 'big') + self.data[16:32]

    def setTypeAndClear(self, Y):
        self.data = self.data[0:16] + Y.to_bytes(4, 'big') + bytes(12)

    def setKeyPairAddresss(self, i):
        self.data = self.data[0:20]+ i.to_bytes(4, 'big') + self.data[24:32]

    def setChainAddress(self, i):
        self.data = self.data[0:24]+ i.to_bytes(4, 'big') + self.data[28:32]

    def setTreeHeight(self, i):
        self.data = self.data[0:24]+ i.to_bytes(4, 'big') + self.data[28:32]

    def setHashAddress(self, i):
        self.data = self.data[0:28]+ i.to_bytes(4, 'big')

    def setTreeIndex(self, i):
        self.data = self.data[0:28]+ i.to_bytes(4, 'big')

    def getKeyPairAddress(self):
        return int.from_bytes(self.data[20:24], 'big')

    def getTreeIndex(self):
        return int.from_bytes(self.data[28:32], 'big')

class slh():

    types = {"WOTS_HASH":0, "WOTS_PK":1, "TREE":2, "FORS_TREE":3, "FORS_ROOTS":4, "WOTS_PRF":5, "FORS_PRF":6}

    def __init__(self, n, h, d, hprime, a, k, lg_w, m):
        self.n = n
        self.h = h
        self.d = d
        self.hprime = hprime
        self.a = a
        self.k = k
        self.lg_w = lg_w
        self.m = m

        self.w = 2**lg_w
        self.len1 = int(math.ceil(8*n/lg_w))
        self.len2 = int(math.log2(self.len1*(self.w-1))//self.lg_w + 1)
        self.length = int(self.len1 + self.len2)

        self.keygen()

    def setSK(self, SK):
        self.SK = SK

    def setPK(self, PK):
        self.PK = PK

    def keygen(self):
        self.SK = slhsk(self.n)
        self.PK = slhpk(self.n)
        ADRS = slhadrs()
        ADRS.setLayerAddres(self.d - 1)
        self.PK.root = self.xmss_node(0, self.hprime, ADRS)

    def sign(self, M, ctx):
        print("signing...")
        if len(ctx) > 255:
            return None
        addrnd = os.urandom(self.n)
        Mprime = bytes(1) + len(ctx).to_bytes(1, 'big') + ctx + M

        ADRS = slhadrs()
        opt_rand = addrnd
        R = self.PRF_msg(opt_rand, Mprime)
        SIG = R
        #print(f"len(R) {len(R)}")
        digest = self.H_msg(R, Mprime)
        md = digest[0 : math.ceil(self.k*self.a/8)]
        tmp_idx_tree = digest[math.ceil(self.k*self.a/8) : math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8)]
        tmp_idx_leaf = digest[math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8) : math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8) + math.ceil(self.h/8/self.d)]
        idx_tree = int.from_bytes(tmp_idx_tree, 'big') % 2**(self.h-self.h//self.d)
        idx_leaf = int.from_bytes(tmp_idx_leaf, 'big') % 2**(self.h//self.d)
        ADRS.setTreeAddress(idx_tree)
        ADRS.setTypeAndClear(self.types["FORS_TREE"])
        ADRS.setKeyPairAddresss(idx_leaf)
        SIG_FORS = self.fors_sign(md, ADRS)
        #print(f"len(SIG_FORS) {len(SIG_FORS)}")
        SIG += SIG_FORS
        PK_FORS = self.fors_pkFromSig(SIG_FORS, md, ADRS)
        SIG_HT = self.ht_sign(PK_FORS, idx_tree, idx_leaf)
        #print(f"len(SIG_HT) {len(SIG_HT)}")
        SIG += SIG_HT
        return SIG

    def verify(self, M, SIG, ctx):
        print("verifying...")
        if len(ctx) > 255:
            print("context too long")
            return False
        Mprime = bytes(1) + len(ctx).to_bytes(1, 'big') + ctx + M

        if len(SIG) != self.n*(1 + self.k*(1 + self.a) + self.h + self.d*self.length):
            print("signature length incorrect")
            print(f"is {len(SIG)} should be {self.n*(1 + self.k*(1 + self.a) + self.h + self.d*self.length)}")
            return False
        ADRS = slhadrs()
        R = SIG[0 : self.n]
        SIG_FORS = SIG[self.n : self.n * (1+self.k*(1+self.a))]
        SIG_HT = SIG[self.n * (1+self.k*(1+self.a)) :  self.n * (1 + self.k*(1 + self.a) + self.h + self.d*self.length)]
        digest = self.H_msg(R, Mprime)
        md = digest[0 : math.ceil(self.k*self.a/8)]
        tmp_idx_tree = digest[math.ceil(self.k*self.a/8) : math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8)]
        tmp_idx_leaf = digest[math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8) : math.ceil(self.k*self.a/8) + math.ceil((self.h - self.h/d)/8) + math.ceil(self.h/8/self.d)]
        idx_tree = int.from_bytes(tmp_idx_tree, 'big') % 2**(self.h-self.h//self.d)
        idx_leaf = int.from_bytes(tmp_idx_leaf, 'big') % 2**(self.h//self.d)
        ADRS.setTreeAddress(idx_tree)
        ADRS.setTypeAndClear(self.types["FORS_TREE"])
        ADRS.setKeyPairAddresss(idx_leaf)
        PK_FORS = self.fors_pkFromSig(SIG_FORS, md, ADRS)
        return self.ht_verify(PK_FORS, SIG_HT, idx_tree, idx_leaf)

    def getPK(self):
        return self.PK

    def getSK(self):
        return self.SK

    def xmss_node(self, i, z, adrs):
        if z == 0:
            adrs.setTypeAndClear(self.types["WOTS_HASH"])
            adrs.setKeyPairAddresss(i)
            node = self.wots_pkGen(adrs)
        else:
            lnode = self.xmss_node(2*i, z - 1, adrs)
            rnode = self.xmss_node(2*i + 1, z - 1, adrs)
            adrs.setTypeAndClear(self.types["TREE"])
            adrs.setTreeHeight(z)
            adrs.setTreeIndex(i)
            node = self.H(adrs, lnode + rnode)
        return node

    def xmss_sign(self, M, idx, adrs):
        auth = b""
        for j in range(0, self.hprime):
            k = (idx//2**j) ^ 1
            auth += self.xmss_node(k, j, adrs)
        adrs.setTypeAndClear(self.types["WOTS_HASH"])
        adrs.setKeyPairAddresss(idx)
        sig = self.wots_sign(M, adrs)
        SIG_XMSS = sig + auth
        return SIG_XMSS

    def xmss_pkFromSig(self, idx, SIG_XMSS, M, adrs):
        adrs.setTypeAndClear(self.types["WOTS_HASH"])
        adrs.setKeyPairAddresss(idx)
        sig = SIG_XMSS[0 : self.length*self.n]
        auth = SIG_XMSS[self.length*self.n : self.n*(self.length + self.hprime)]
        node0 = self.wots_pkFromSig(sig,  M, adrs)
        adrs.setTypeAndClear(self.types["TREE"])
        adrs.setTreeIndex(idx)
        for k in range(0, self.hprime):
            adrs.setTreeHeight(k+1)
            if (idx//2**k) % 2 == 0:
                adrs.setTreeIndex(adrs.getTreeIndex()//2)
                node1 = self.H(adrs, node0 + auth[k*self.n : (k + 1)*self.n])
            else:
                adrs.setTreeIndex((adrs.getTreeIndex() - 1)//2)
                node1 = self.H(adrs, auth[k*self.n : (k + 1)*self.n] + node0)
            node0 = node1
        return node0

    def wots_pkGen(self, adrs):
        tmp = b""
        skADRS = slhadrs(adrs)
        skADRS.setTypeAndClear(self.types["WOTS_PRF"])
        skADRS.setKeyPairAddresss(adrs.getKeyPairAddress())
        for i in range(0, self.length):
            skADRS.setChainAddress(i)
            sk = self.PRF(skADRS)
            adrs.setChainAddress(i)
            tmp += self.chain(sk, 0, self.w - 1, adrs)
        wotspkADRS = slhadrs(adrs)
        wotspkADRS.setTypeAndClear(self.types["WOTS_PK"])
        wotspkADRS.setKeyPairAddresss(adrs.getKeyPairAddress())
        pk = self.T(wotspkADRS, tmp)
        return pk

    def wots_sign(self, M, adrs):
        csum = 0
        msg = self.base_2b(M, self.lg_w, self.len1)
        for i in range(0, self.len1):
            csum += self.w - 1 - msg[i]
        csum <<= (8 - ((self.len2*self.lg_w) % 8)) % 8
        msg += self.base_2b(csum.to_bytes(math.ceil(self.len2*self.lg_w/8), 'big'), self.lg_w, self.len2)
        skadrs = slhadrs(adrs)
        skadrs.setTypeAndClear(self.types["WOTS_PRF"])
        skadrs.setKeyPairAddresss(adrs.getKeyPairAddress())
        sig = b""
        for i in range(0, self.length):
            skadrs.setChainAddress(i)
            sk = self.PRF(skadrs)
            adrs.setChainAddress(i)
            sig += self.chain(sk, 0, msg[i], adrs)
        return sig

    def wots_pkFromSig(self, sig, M, adrs):
        csum = 0
        msg = self.base_2b(M, self.lg_w, self.len1)
        for i in range(0, self.len1):
            csum += self.w - 1 - msg[i]
        csum <<= (8 - (self.len2*self.lg_w % 8)) % 8
        msg += self.base_2b(csum.to_bytes(math.ceil(self.len2*self.lg_w/8), 'big'), self.lg_w, self.len2)
        tmp = b""
        for i in range(0, self.length):
            adrs.setChainAddress(i)
            tmp += self.chain(sig[i*self.n : (i + 1)*self.n], msg[i], self.w - 1 - msg[i], adrs)
        wotspkadrs = slhadrs(adrs)
        wotspkadrs.setTypeAndClear(self.types["WOTS_PK"])
        wotspkadrs.setKeyPairAddresss(adrs.getKeyPairAddress())
        pk_sig = self.T(wotspkadrs, tmp)
        return pk_sig

    def chain(self, X, i, s, adrs):
        tmp = X
        for j in range(i, i + s):
            adrs.setHashAddress(j)
            tmp = self.F(adrs, tmp)
        return tmp

    def fors_skGen(self, adrs, idx):
        skadrs = slhadrs(adrs)
        skadrs.setTypeAndClear(self.types["FORS_PRF"])
        skadrs.setKeyPairAddresss(adrs.getKeyPairAddress())
        skadrs.setTreeIndex(idx)
        return self.PRF(skadrs)

    def fors_sign(self, md, adrs):
        SIG_FORS = b""
        indices = self.base_2b(md, self.a, self.k)
        for i in range(self.k):
            SIG_FORS += self.fors_skGen(adrs, i*2**self.a + indices[i])
            AUTH = b""
            for j in range(self.a):
                s = (indices[i]//2) ^ 1
                AUTH += self.fors_node(i*2**(self.a - j) + s, j, adrs)
            SIG_FORS += AUTH
        return SIG_FORS

    def fors_node(self, i, z, adrs):
        if z == 0:
            sk = self.fors_skGen(adrs, i)
            adrs.setTreeHeight(0)
            adrs.setTreeIndex(i)
            node = self.F(adrs, sk)
        else:
            lnode = self.fors_node(2*i, z - 1, adrs)
            rnode = self.fors_node(2*i + 1, z - 1, adrs)
            adrs.setTreeHeight(z)
            adrs.setTreeIndex(i)
            node = self.H(adrs, lnode + rnode)
        return node

    def fors_pkFromSig(self, SIG_FORS, md, adrs):
        indices = self.base_2b(md, self.a, self.k)
        root = b""
        for i in range(0, self.k):
            sk = SIG_FORS[i*(self.a + 1)*self.n : (i*(self.a + 1) + 1)*self.n]
            adrs.setTreeHeight(0)
            adrs.setTreeIndex(i*2**self.a + indices[i])
            node0 = self.F(adrs, sk)
            auth = SIG_FORS[(i*(self.a + 1) + 1)*self.n : (i + 1)*(self.a + 1)*self.n]
            for j in range(0, self.a):
                adrs.setTreeHeight(j + 1)
                if (indices[i]// 2**j) % 2 == 0:
                    adrs.setTreeIndex(adrs.getTreeIndex()//2)
                    node1 = self.H(adrs, node0 + auth[j*n : (j + 1)*n])
                else:
                    adrs.setTreeIndex((adrs.getTreeIndex() - 1) // 2)
                    node1 = self.H(adrs, auth[j*n : (j + 1)*n] + node0)
                node0 = node1
            root += node0
        forspkadrs = slhadrs(adrs)
        forspkadrs.setTypeAndClear(self.types["FORS_ROOTS"])
        forspkadrs.setKeyPairAddresss(adrs.getKeyPairAddress())
        pk = self.T(forspkadrs, root)
        return pk

    def ht_sign(self, M, idx_tree, idx_leaf):
        ADRS = slhadrs()
        ADRS.setTreeAddress(idx_tree)
        SIG_tmp = self.xmss_sign(M, idx_leaf, ADRS)
        SIG_HT = SIG_tmp
        root = self.xmss_pkFromSig(idx_leaf, SIG_tmp, M, ADRS)
        for j in range(1, self.d):
            idx_leaf = idx_tree % 2**self.hprime
            idx_tree >>= self.hprime
            ADRS.setLayerAddres(j)
            ADRS.setTreeAddress(idx_tree)
            SIG_tmp = self.xmss_sign(root, idx_leaf, ADRS)
            SIG_HT += SIG_tmp
            if j < d - 1:
                root = self.xmss_pkFromSig(idx_leaf, SIG_tmp, root, ADRS)
        return SIG_HT

    def ht_verify(self, M, SIG_HT, idx_tree, idx_leaf):
        ADRS = slhadrs()
        ADRS.setTreeAddress(idx_tree)
        SIG_tmp = SIG_HT[0 : (self.hprime + self.length)*self.n]
        node = self.xmss_pkFromSig(idx_leaf, SIG_tmp, M, ADRS)
        for j in range(1, self.d):
            idx_leaf = idx_tree % 2**self.hprime
            idx_tree >>= self.hprime
            ADRS.setLayerAddres(j)
            ADRS.setTreeAddress(idx_tree)
            SIG_tmp = SIG_HT[j*(self.hprime + self.length)*self.n : (j + 1)*(self.hprime + self.length)*self.n]
            node = self.xmss_pkFromSig(idx_leaf, SIG_tmp, node, ADRS)
            #print(f"node {j}: {node.hex()}")
        if node == self.PK.root:
            return True
        else:
            return False

    def base_2b(self, X, b, out_len):
        ii = 0
        bits = 0
        total = 0
        baseb = []
        for out in range(out_len):
            while bits < b:
                total = (total << 8) + int(X[ii])
                ii += 1
                bits += 8
            bits -= b
            baseb.append((total >> bits) % 2**b)
        return baseb

    # NOTE:  SHAKE256.read(n) gets n bytes of output, whereas output length is in BITS in specs

    def H(self, adrs, M_2):
        shake = SHAKE256.new()
        shake.update(self.PK.seed + adrs.data + M_2)
        return shake.read(self.n)

    def H_msg(self, R, M):
        shake = SHAKE256.new()
        shake.update(R + self.PK.seed + self.PK.root + M)
        return shake.read(self.m)

    def PRF(self, adrs):
        shake = SHAKE256.new()
        shake.update(self.PK.seed + adrs.data +  self.SK.seed)
        return shake.read(self.n)

    def PRF_msg(self, opt_rand, M):
        shake = SHAKE256.new()
        shake.update(self.SK.prf + opt_rand + M)
        return shake.read(self.n)

    def F(self, adrs, M_1):
        shake = SHAKE256.new()
        shake.update(self.PK.seed + adrs.data + M_1)
        return shake.read(self.n)

    def T(self, adrs, M):
        shake = SHAKE256.new()
        shake.update(self.PK.seed + adrs.data + M)
        return shake.read(self.n)

    def print_sig(self, sig):
        ptr = 0
        print("R:", sig[ptr : ptr + self.n].hex())
        print("")
        ptr += self.n
        print(f"SIG_FORS (k={self.k}, a={self.a}):")
        for i in range(self.k):
            print(f"\tFORS private key value {i}: {sig[ptr : ptr + self.n].hex()}")
            ptr += self.n
            print(f"\t\tFORS authpath[{i}]:")
            for j in range(self.a):
                print("\t\t\t", sig[ptr : ptr + self.n].hex())
                ptr += self.n
        print("")
        print(f"SIG_HT (d={self.d}):")
        for i in range(self.d):
            print(f"\tHT layer {i} XMSS sig:")
            print(f"\t\tSIG_WOTS+ (hprime={self.hprime}): ", sig[ptr : ptr + self.length*self.n].hex())
            ptr += self.length*self.n
            for j in range(self.hprime):
                print(f"\t\t\tXMSS authpath[{j}]: ", sig[ptr : ptr + self.n].hex())
                ptr += self.n

#-----------------------------------#
# test, parse sig
#-----------------------------------#

if __name__ == "__main__":
    # 128s params
    n = 16
    h = 63 # = hprime * d
    d = 7
    hprime = 9
    a = 12
    k = 14
    lg_w = 4
    m = 30
    myslh = slh(n, h, d, hprime, a, k, lg_w, m)
    #print(myslh.length)

    print("*128s params non-deterministic pure SLH_DSA signature*")
    print("")
    msg = os.urandom(32)
    ctx = os.urandom(32)
    print("message (32 bytes hex): ", msg.hex())
    print("context (32 bytes hex): ", ctx.hex())
    print("")
    sig = myslh.sign(msg, ctx)
    ver = myslh.verify(msg, sig, ctx)
    print("verified: ", ver)
    print("")
    print("parsing signature...")
    print("")
    myslh.print_sig(sig)


