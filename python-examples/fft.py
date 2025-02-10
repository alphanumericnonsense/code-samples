import math
import random

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

def fft_recur(f,n):
    zeta_n = math.cos(2*math.pi/n)+1j*math.sin(2*math.pi/n)
    if n == 1:
        return f
    f_even = [f[2*i] for i in range(n//2)]
    f_odd = [f[2*i + 1] for i in range(n//2)]
    f_front = [fft_recur(f_even, n//2)[i] + (zeta_n**i)*fft_recur(f_odd, n//2)[i] for i in range(n//2)]
    f_back = [fft_recur(f_even, n//2)[i] - (zeta_n**i)*fft_recur(f_odd, n//2)[i] for i in range(n//2)]
    return f_front + f_back

def ifft_recur(f,n):
    """
    introduces 1/n one power of two at each depth
    """
    zeta_n = math.cos(2*math.pi/n)-1j*math.sin(2*math.pi/n)
    if n == 1:
        return [f[0]]
    f_even = [f[2*i] for i in range(n//2)]
    f_odd = [f[2*i + 1] for i in range(n//2)]
    f_front = [ifft_recur(f_even, n//2)[i]/2 + (zeta_n**i)*ifft_recur(f_odd, n//2)[i]/2 for i in range(n//2)]
    f_back = [ifft_recur(f_even, n//2)[i]/2 - (zeta_n**i)*ifft_recur(f_odd, n//2)[i]/2 for i in range(n//2)]
    return f_front + f_back

def fft_matmult(f,n):
    zeta_n = math.cos(2*math.pi/n)+1j*math.sin(2*math.pi/n)
    M = [[zeta_n**(i*j)for j in range(n)] for i in range(n)]
    return [sum([M[i][k]*f[k] for k in range(n)]) for i in range(n)]

def ifft_matmult(f,n):
    zeta_n = math.cos(2*math.pi/n)-1j*math.sin(2*math.pi/n)
    M = [[zeta_n**(i*j)for j in range(n)] for i in range(n)]
    return [sum([M[i][k]*f[k] for k in range(n)]) for i in range(n)]

def test(n):
    """
    n should be a power of two
    """
    f = [complex(random.random(), random.random()) for i in range(n)]
    fhatiter = fft_iter(f,n)
    fhatrecur = fft_recur(f,n)
    fiter = ifft_iter(fhatiter,n)
    frecur = ifft_recur(fhatrecur,n)
    frecur = ifft_recur(fhatrecur,n)
    print(f"input f length {n}")
    print(f)
    print("transform should be")
    print(fft_matmult(f,n))
    print("fhatiter")
    print(fhatiter)
    print("fhatrecur")
    print(fhatrecur)
    print("inverting...")
    print(fiter)
    print(frecur)
    print('\n')

if __name__ == "__main__":
    for i in [1,2,4,8]:
        test(i)
