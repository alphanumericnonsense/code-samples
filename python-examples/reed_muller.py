# For binary Reed-Muller codes, R(r,m), r <= m
# [n,k,d] = [2^m, sum_{i=0}^r binom(m,i), 2^{m-r}]
#
# works, need to invert generator matrix to get back to message?

import numpy as np
import math
import random

def binom(n,k):
    """
    return binomial coefficient n choose k using additive recursion
    """
    if k < 0 or n < 0:
        return 0
    elif k == 0:
        return 1
    else:
        return binom(n-1,k-1) + binom(n-1,k)

def getDim(r,m):
    """
    return dimension of R(r,m)
    """
    S = 0
    for i in range(r+1):
        S += binom(m, i)
    return S

def getGen(r,m):
    """
    return generator matrix G for R(r,m)
    """
    #print(r,m)
    if r == m:
        return np.identity(2**m)
    elif r == 0:
        return np.array([1 for i in range(2**m)])
    else:
        k = getDim(r,m)
        G = np.zeros((k,2**m))
        A = getGen(r, m-1)
        B = A

        kk = getDim(r-1,m-1)
        C = np.zeros((kk, 2**(m-1)))
        D = getGen(r-1, m-1)

        kkk = getDim(r, m-1)

        #print(k,kk,kkk)

        G[0:kkk, 0:2**(m-1)] = A
        G[0:kkk,2**(m-1):] = B
        G[kkk:,0:2**(m-1)] = C
        G[kkk:,2**(m-1):] = D
        return G

def encode(message,r,m):
    """
    return codeword for message, mG
    """
    return np.mod(np.matmul(np.array(message), getGen(r,m)),2)

def getBin(i,m):
    out = np.zeros(m)
    for j in range(m):
        out[m-j-1] = i&1
        i = i >> 1
    return out

def evalTransMono(b,v,m):
    """
    return evaluation vector of length 2^m,
    characteristic function of translated flat b + F_v
    """
    out = np.zeros(2**m)
    for i in range(2**m):
        coord = getBin(i,m)
        trans = np.mod(coord + b, 2)
        put1 = True
        for j in range(m):
            if trans[j] == 0 and v[j] == 1:
                put1 = False
                break
        if put1:
            out[i] = 1
    return out

def monomials(t,m):
    """
    return all weight t bitstrings of length m
    """
    biglist = []

    if t == m:
        biglist.append([1 for i in range(m)])
    elif t == 0 and m > 0:
        biglist.append([0 for i in range(m)])
    elif 0 < t < m:
        for x in monomials(t,m-1):
            biglist.append([0]+x)
        for x in monomials(t-1,m-1):
            biglist.append([1]+x)
    else:
        print("bad case?!", t, m)
    return biglist


def rmError(received,r,m):
    """
    return codeword closest to received (assuming < radius/2 errors).
    majority logic,
    inefficient as all 2^m translates are considered (even in the same flat)
    """
    k = getDim(r,m)
    out = np.zeros(2**m)
    received_tmp = np.zeros(2**m)
    for t in range(r,-1,-1):
        received_tmp = np.mod(received - out, 2)
        #print(t, received_tmp)
        L = monomials(t,m)
        for v in L:
            v = np.array(v)
            vcomp = 1-v
            basis_vec_img = evalTransMono(np.zeros(m), v, m)
            #print(v, vcomp, basis_vec_img)
            vote0 = 0
            vote1 = 0
            for j in range(2**m):
                b = getBin(j,m)
                A_b_vcomp = evalTransMono(b,vcomp,m)
                vote = np.mod(np.dot(received_tmp, A_b_vcomp), 2)
                #print(vote, A_b_vcomp)
                if vote == 0:
                    vote0 += 1
                else:
                    vote1 += 1
            if vote1 > vote0:
                out = np.mod(out + basis_vec_img, 2)
            #print(out)
    return out

def randomError(thresh, m):
    """
    return an error vector of weight <= thresh
    """
    t = math.floor((thresh+1)*random.random())
    return np.random.permutation(np.array([1 for i in range(t)]+[0 for j in range(2**m-t)]))

def randomMessage(r, m):
    """
    return a uniformly random message
    """
    k = getDim(r, m)
    return np.array([math.floor(2*random.random()) for i in range(k)])

def testEncDec(num_trials, r, m):
    """
    encode/decode random messages with errors within half the minimum distance
    """
    for i in range(num_trials):
        message = randomMessage(r, m)
        #print("message:", message)
        codeword = encode(message, r, m)
        #print("encoded:", codeword)
        error = randomError(2**(m-r-1) - 1, m)
        print("number of errors:", np.sum(error))
        received = np.mod(codeword + error, 2)
        #print("received:", received)
        corrected = rmError(received, r, m)
        #print("corrected:", corrected)
        if np.all(corrected == codeword):
            print("\tSUCCESS!")
        else:
            print("\tFAILURE!")

