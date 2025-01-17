import copy

def PolyMult(P, Q):
    R = [0 for k in range(len(P)+len(Q)-1)]
    for i in range(len(P)):
        for j in range(len(Q)):
            R[i+j] += P[i]*Q[j]
    return R

def PolyScale(P, c):
    return [c*p for p in P]

def PolyAdd(P, Q):
    M = max(len(P),len(Q))
    R = [0 for i in range(M)]
    for i in range(M):
        tmp = 0
        if i < len(P):
            tmp += P[i]
        if i < len(Q):
            tmp += Q[i]
        R[i] = tmp
    return R

def Monomial(d):
    return [0 for i in range(d)] + [1]

def BM(S):
    """
    Return connection poly for shortest LFSR producing sequence S,
    sum_i C_i*S_{N-i} = 0 for L <= N < n, C_0 = 1, [1/C = S?]
    """
    C = [1] # poly over field, connection poly to be returned
    B = [1] #poly over field
    x = 1 # integer
    L = 0 # integer, length of "connection poly"/recurrence relation/LFSR
    b = 1 # field elt
    N = 0 # integer, counts 0 to n - 1
    n = len(S)

    while N < n:
        d = 0 # field elt, discrepancy
        for i in range(L+1):
            d += C[i]*S[N-i]
        if d == 0:
            x += 1
        elif 2*L > N:
            C = PolyAdd(C, PolyScale(PolyMult(Monomial(x), B), -d/b))
            x += 1
        else:
            tmp = copy.copy(C)
            C = PolyAdd(C, PolyScale(PolyMult(Monomial(x), B), -d/b))
            B = tmp
            L = N + 1 - L
            b = d
            x = 1
        N += 1
    return C

def LFSR(C, S, n):
    """
    produce first n + len(S) elements of LFSR sequence with initial terms S,
    connection poly C, C_0 = 1
    """
    for i in range(n):
        tot = 0
        for i in range(1,len(C)):
            tot -= C[i]*S[len(S)-i]
        S.append(tot)
    return S
