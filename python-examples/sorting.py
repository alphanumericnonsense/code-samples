import math
from random import random
import time

def bubblesort(A):
    for j in range(len(A)):
        idle = True
        for i in range(len(A)-1):
            if A[i] > A[i+1]:
                A[i], A[i+1] = A[i+1], A[i]
                idle = False
        if idle:
            break
    return A

def merge(A, B):
    result = []
    while len(A) != 0 and len(B) != 0:
        if A[0] <= B[0]:
            result.append(A[0])
            A.remove(A[0])
        else:
            result.append(B[0])
            B.remove(B[0])
    return result+A+B

def mergesort(A):
    if len(A) <= 1:
        return A
    split = len(A)//2
    A1 = mergesort(A[:split])
    A2 = mergesort(A[split:])
    return merge(A1, A2)

def quicksort(A):
    if len(A) <= 1:
        return A
    p = (A[0]+A[-1])/2
    left = []
    right = []
    for a in A:
        if a <= p:
            left.append(a)
        else:
            right.append(a)

    if right == []:
        right.append(left[-1])
        left.remove(left[-1])
    elif left == []:
        left.append(right[-1])
        right.remove(right[-1])
    #print(left, right)
    return quicksort(left)+quicksort(right)

def selectionsort(A):
    if len(A) <= 1:
        return A
    m = min(A)
    A.remove(m)
    return [m] + selectionsort(A)

def insertionsort(A):
    result = []
    for a in A:
        idle = True
        for i in range(len(result)):
            if a < result[i]:
                result = result[:i+1] + [a] + result[i+1:]
                idle = False
                break
        if idle:
            result.append(a)
    return result

#--------------#
# heapsort
#--------------#
def heapify(A, i, heaplen):
    largest = i
    l = 2*i+1
    r = 2*i + 2
    if l < heaplen:
        if A[l] > A[i]:
            largest = l
    if r < heaplen:
        if A[r] > A[largest]:
            largest = r
    if largest != i:
        A[i], A[largest] = A[largest], A[i]
        heapify(A, largest, heaplen)

def build_heap(A, heaplen):
    for i in range(len(A)//2 - 1, -1, -1):
        heapify(A, i, heaplen)
    #print(A)

def heapsort(A):
    heaplen = len(A)
    build_heap(A, heaplen)
    for i in range(len(A)-1,0,-1):
        A[0], A[i] = A[i], A[0]
        heaplen = heaplen - 1
        heapify(A, 0, heaplen)

def bitslice_shiftrank(A,n,k):
    """
    Sorts array A of n k-bit unsigned integers
    """
    rankmat = [0 for i in range(n)]
    for i in range(k):
        updates = [0 for i in range(n)]
        slicei = [(a >> (k-1-i)) & 1 for a in A]
        print(slicei)
        for j in range(n):
            if slicei[j]:
                for l in range(n):
                    if rankmat[l] > rankmat[j]:
                        updates[l] += 1
                updates[j] += 1
        rankmat = [rankmat[i] + updates[i] for i in range(n)]
    return rankmat

def sort_test(size, bound, rounds):
    for sort in [bubblesort, mergesort, quicksort, selectionsort, insertionsort, heapsort]:
        t0 = time.time()
        for i in range(rounds):
            A = [math.floor(bound*random()) for i in range(size)]
            sort(A)
        t1 = time.time()
        print(sort.__name__, '\t', t1-t0)

if __name__ == "__main__":
    size = 100
    bound = 100
    rounds = 100
    print(f"timing {rounds} sorts of random lists of length {size} with random integers in range [0,{bound})")
    sort_test(100,100,100)
