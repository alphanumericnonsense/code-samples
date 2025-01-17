import itertools
import math

def cycle_counts_alt1(g):
    """
    return a list [z1(g), ..., zn(g)]
    where zi(g) is the number of cycles of length i in g, g in S_n.
    """
    n = len(g)
    left = set(g)
    z = {i:0 for i in range(1,n+1)}
    while left != set():
        #print(left)
        left_min = min(left)
        left -= {left_min}
        cyc_next = g[left_min]
        cyc_len = 1
        while cyc_next != left_min:
            left -= {cyc_next}
            cyc_next = g[cyc_next]
            cyc_len += 1
        z[cyc_len] += 1
    return z

def cycle_counts_alt2(g):
    """
    return a dict {1:z1(g), ..., n:zn(g)}
    where zi(g) is the number of cycles of length i in g, g in S_n,
    !!! only for non-zero zi(g) !!!
    """
    n = len(g)
    left = set(g)
    z = {} # {i:0 for i in range(1,n+1)}
    while left != set():
        #print(left)
        left_min = min(left)
        left -= {left_min}
        cyc_next = g[left_min]
        cyc_len = 1
        while cyc_next != left_min:
            left -= {cyc_next}
            cyc_next = g[cyc_next]
            cyc_len += 1
        if cyc_len in z:
            z[cyc_len] += 1
        else:
            z[cyc_len] = 1
    return z

def gcd(a,b):
    while b!=0:
        a,b=b,a%b
    return abs(a)

# def factorial(n):
#     if n == 0:
#         return 1
#     else:
#         return n*factorial(n-1)

def factorial(n):
    prod = 1
    while n > 0:
        prod *= n
        n -=1
    return prod

def solution1(w, h, s):
    """
    counting the number of orbits of
    G = S_w\times S_h acting on X = Mat^{w\times h}(Z/sZ),
    (a,b)*m = a*m*b^{-1}.

    Burnside's lemma:
    |G\X| = \frac{1}{|G|}\sum_{g\in G}|Fix(g)|

    rephrased:
    In the Polya enumeration theory, we have S_w\times S_h acting on
    functions from the product {1,2,...,w}\times{1,2,...,h}
    to {1,2,...,s} by independently permuting elements of domain,
    (g,h)*(i,j) = (g(i), h(j)).

    Frank Harary, "ON THE NUMBER OF BI-COLORED GRAPHS",
    Pacific Journal of Mathematics, vol. 4 no. 8, June 1958,
    describes cycle index of product acting on product.

    Evaluate cycle index at s to get orbit count.

    Plenty of ways to be more efficient.
    Really only need cycle types with multiplicities, not all perms.

    Think I'm timing out on some test cases?
    Passes 1, 2, 4, and 9...
    """
    G1 = list(itertools.permutations(range(w)))
    G2 = list(itertools.permutations(range(h)))
    mysum = 0
    for g1 in G1:
        z_g1 = cycle_counts_alt1(g1)
        for g2 in G2:
            z_g2 = cycle_counts_alt1(g2)
            myprod = 1
            for r1 in range(1,w+1):
                for r2 in range(1,h+1):
                    #myprod *= s**(gcd(r1,r2)*num_cycles(r1,g1)*num_cycles(r2,g2))
                    myprod *= s**(gcd(r1,r2)*z_g1[r1]*z_g2[r2])
            #print(g1, g2)
            #print(myprod)
            mysum += myprod
    return str(mysum//factorial(w)//factorial(h))

def solution2(w, h, s):
    """
    Reduced loop sizes and replaced some mult with add.
    Still too slow (or the solution is wrong).

    Passes 1,2,4,9
    """
    G1 = list(itertools.permutations(range(w)))
    G2 = list(itertools.permutations(range(h)))
    mysum = 0
    for g1 in G1:
        z_g1 = cycle_counts_alt2(g1)
        for g2 in G2:
            z_g2 = cycle_counts_alt2(g2)
            myexp = 0
            for r1 in z_g1:
                for r2 in z_g2:
                    myexp += gcd(r1,r2)*z_g1[r1]*z_g2[r2]
            mysum += s**myexp
            #print(mysum)
    return str(mysum//factorial(w)//factorial(h))

def solution(w, h, s):
    """
    evaluate "cycle index" for the product of symmetric groups S_w\times S_h
    acting on the product [1,w]\times[1,h] diagonally, (p1, p2) * (i, j) = (p1(i), p2(j)).

    Sources:
        Van Lint and Wilson, "A Course in Combinatorics" (i.e a book I have laying around).
        Harary, "On the number of bi-colored graphs", for the cycle index of product acting on product.
        enumerating partitions, https://www.geeksforgeeks.org/generate-unique-partitions-of-an-integer/

    Definitely faster, disagrees with other solutions sometimes???
    Only passes 1,2,4
    """
    mysum = 0
    G1 = partition_generator(w)
    for g1 in G1:
        z1 = cycle_counts(g1)
        m1 = cycle_type_mult(g1,z1) # multiplicity of cycle type
        G2 = partition_generator(h)
        for g2 in G2:
            z2 = cycle_counts(g2)
            m2 = cycle_type_mult(g2,z2)# multiplicity of cycle type
            myexp = 0
            for r1 in z1:
                for r2 in z2:
                    myexp += gcd(r1,r2)*z1[r1]*z2[r2]
            mysum += m1*m2*s**myexp
            #print(g1,g2)
            #print(m1,m2)
    return str(mysum//factorial(w)//factorial(h))

def cycle_type_mult(g,z):
    myprod = factorial(sum(g)) # n!
    #myprod = myprod//math.prod(g) # math.prod not available until python 3.8???
    listprod = 1
    for part in g:
        listprod *= part
    myprod = myprod//listprod # divide out by cyclic shifts of each cycle in perm
    for i in z:
        myprod = myprod//factorial(z[i]) # divide out by permuting cycles of same size
    return myprod

def cycle_counts(p):
    """
    returns dict {k:z_k} with z_k the number of parts of size k in p
    """
    z = {}
    for part in p:
        if part in z:
            z[part] += 1
        else:
            z[part] = 1
    return z

def partition_generator_bad(n):
    """
    yields parition of n as list of non-increasing positive integers,
    in order, e.g. [4], [3,1], [2,2], [2,1,1], [1,1,1,1].

    Code modified from https://www.geeksforgeeks.org/generate-unique-partitions-of-an-integer/
    """
    p = [0] * n
    k = 0
    p[k] = n
    while True:
        yield [part for part in p if part != 0]
        # Generate next partition

        # Find the rightmost non-one value in p[].
        # Also, update the rem_val so that we know
        # how much value can be accommodated
        rem_val = 0
        while k >= 0 and p[k] == 1:
            rem_val += p[k]
            k -= 1
        # if k < 0, all the values are 1 so
        # there are no more partitions
        if k < 0:
            return
        # Decrease the p[k] found above
        # and adjust the rem_val
        p[k] -= 1
        rem_val += 1
        # If rem_val is more, then the sorted
        # order is violated. Divide rem_val in
        # different values of size p[k] and copy
        # these values at different positions after p[k]
        while rem_val > p[k]:
            p[k + 1] = p[k]
            rem_val = rem_val - p[k]
            k += 1
        # Copy rem_val to next position
        # and increment position
        p[k + 1] = rem_val
        k += 1

def partition_generator(n):
    """
    taken from https://jeromekelleher.net/generating-integer-partitions.html
    """
    a = [0 for i in range(n + 1)]
    k = 1
    a[1] = n
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while x <= y:
            a[k] = x
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]

def print_solutions(N):
    for i in range(1,N):
        for j in range(1,N):
            for k in range(1,N):
                #print(f"test ({i},{j},{k})") # formatted strings not available in python 2.7???
                print(solution1(i,j,k), solution2(i,j,k), solution(i,j,k))

def test_solutions(N):
    for i in range(1,N):
        for j in range(1,N):
            for k in range(1,N):
                s1 = solution1(i,j,k)
                s2 = solution2(i,j,k)
                s3 = solution(i,j,k)
                if s1 != s2 or s1 != s3 or s2 != s3:
                    print("whs: ",i,j,k)
                    print("sols: ", s1,s2,s3)
