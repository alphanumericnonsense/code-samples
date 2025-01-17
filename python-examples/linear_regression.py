def least_squares(X, Y):
    """
    Return (a,b), y = ax+b the least square line for (x,y) in X x Y
    """
    if len(X) != len(Y):
        print("X and Y not of equal lengths; returning None.")
        return None
    else:
        N = len(X)
    a11, a12, a21, a22, b1, b2 = (0,0,0,0,0,0)
    for i in range(N):
        a11 += X[i]
        a12 += X[i]*X[i]
        a21 += 1
        a22 += X[i]
        b1 += X[i]*Y[i]
        b2 += Y[i]
    A11 = a22
    A12 = -a12
    A21 = -a21
    A22 = a11
    det = a11*a22-a12*a21
    return ((A21*b1+A22*b2)/det, (A11*b1+A12*b2)/det)
