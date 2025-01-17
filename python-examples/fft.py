import math

def bit_reversal(i, d):
    """
    d-bit bit-reversal of i, 0 <= i < 2^d
    """
    s = bin(i)[2:]
    #print(s)
    s = '0'*(d-len(s)) + s
    #print(s)
    result = 0
    for j in range(len(s)):
        result += int(s[j])*2**j
    return result

def fft_iter(f,n):
    zeta_n = math.cos(2*math.pi/n)+1j*math.sin(2*math.pi/n)
    d = round(math.log2(n))
    powers = [zeta_n**(n/2**(i+1)) for i in range(d)]
    g = [f[bit_reversal(i,d)] for i in range(len(f))]
    m = 1
    for i in range(d):
        M = 2*m
        for j in range(0,n,M):
            omega = 1
            for k in range(m):
                a = g[k+j]
                b = omega*g[k+j+m]
                g[k+j] = a + b
                g[k+j+m] = a - b
                omega = omega*powers[i]
        m = M
    return g

def ifft_iter(f,n):
    zeta_n = math.cos(2*math.pi/n)-1j*math.sin(2*math.pi/n)
    d = round(math.log2(n))
    powers = [zeta_n**(n/2**(i+1)) for i in range(d)]
    g = [f[bit_reversal(i,d)] for i in range(len(f))]
    m = 1
    for i in range(d):
        M = 2*m
        for j in range(0,n,M):
            omega = 1
            for k in range(m):
                a = g[k+j]
                b = omega*g[k+j+m]
                g[k+j] = a + b
                g[k+j+m] = a - b
                omega = omega*powers[i]
        m = M
    return [gi/n for gi in g]
