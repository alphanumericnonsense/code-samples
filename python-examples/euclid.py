def binary_EEA(a,b):
    """
    Return (u,v,g) with u*a+v*b = g = gcd(a,b)
    Throughout, a*si+b*ti = ri for i = 0,1.
    Assumes a, b >= 0.

    "Binary" in that it works with odd/even instead of costlier integer division.

    Should also work for a, b representing polynomials over GF(2)
    after replacing addition/subtraction with XOR and reducing degree.
    """
    k = 0
    s0 = 1
    s1 = 0
    t0 = 0
    t1 = 1

    # trivial cases
    if a == 0:
        return a
    elif b == 0:
        return a

    # pull out common power of 2
    while not ((a|b)&1):
        a >>= 1
        b >>= 1
        k += 1

    # odd starting values
    r0 = a
    r1 = b

    while r0 != r1:
        if r0&1 == 1:
            if r0 >= r1:
                r0 = r0 - r1
                s0 = s0 - s1
                t0 = t0 - t1
            else:
                r1 = r1 - r0
                s1 = s1 - s0
                t1 = t1 - t0
        else:
            r0 >>= 1
            if not ((s0|t0)&1):
                s0 >>= 1
                t0 >>= 1
            else:
                s0 = (s0+b)>>1
                t0 = (t0-a)>>1
    return (s0, t0, r0*(1<<k))

def binary_poly_EEA(a,b):
    """
    Return (u,v,g) with u*a+v*b = g = gcd(a,b)
    Throughout, a*si+b*ti = ri for i = 0,1.
    Assumes a, b >= 0.

    Slightly modified from above.
    In this version, a, b are non-negative integers interpreted as polynomials over GF(2),
    e.g. DEC 14 = BIN 1101 = x^3 + x^2 + 1.

    "Binary" in that it works with the property of divisible by x or not ("even/odd")
    instead of using expensive polynomial division.
    """
    k = 0
    s0 = 1
    s1 = 0
    t0 = 0
    t1 = 1

    # trivial cases
    if a == 0:
        return a
    elif b == 0:
        return a

    # pull out common power of 2
    while not ((a|b)&1):
        a >>= 1
        b >>= 1
        k += 1

    # odd starting values
    r0 = a
    r1 = b

    while r0 != r1:
        if r0&1 == 1:
            degr0 = int(math.log2(r0))
            degr1 = int(math.log2(r1))
            if degr0 >= degr1:
                #r0 = r0 - r1
                #s0 = s0 - s1
                #t0 = t0 - t1
                r0 = r0 ^ (r1<<(degr0-degr1))
                s0 = s0 ^ (s1<<(degr0-degr1))
                t0 = t0 ^ (t1<<(degr0-degr1))
            else:
                #r1 = r1 - r0
                #s1 = s1 - s0
                #t1 = t1 - t0
                r1 = r1 ^ (r0<<(degr1-degr0))
                s1 = s1 ^ (s0<<(degr1-degr0))
                t1 = t1 ^ (t0<<(degr1-degr0))
        else:
            r0 >>= 1
            if not ((s0|t0)&1):
                s0 >>= 1
                t0 >>= 1
            else:
                #s0 = (s0+b)>>1
                #t0 = (t0-a)>>1
                s0 = (s0^b)>>1
                t0 = (t0^a)>>1
    return (s0, t0, r0*(1<<k))

def EEA(a,b):
    """
    Return (x, y, g) such that gcd(a,b) = g = a*x + b*y.
    Throughout, a*s_i+b*t_i = r_i for i = old, new.
    Assumes a, b >= 0.

    Standard (a,b) -> (b, a mod b) with integer division.
    """
    r_old, r_new = a, b
    s_old, s_new = 1, 0
    t_old, t_new = 0, 1
    while r_new != 0:
        q = r_old // r_new
        r_old, r_new = r_new, r_old - q*r_new
        s_old, s_new = s_new, s_old - q*s_new
        t_old, t_new = t_new, t_old - q*t_new
    return (s_old, t_old, r_old)
