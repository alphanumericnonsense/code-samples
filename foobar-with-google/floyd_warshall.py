import itertools

def has_timelike_loops(times):
    """
    check for a distinct vertex cycle with negative edge-sum.
    not efficient.
    e.g. [[0,N,-1],[N,0,N],[0,N,0]] requires looping Nish times to save bunny.
    """
    n = len(times)
    # enumerate all cyclic perms of length <= n
    # overcounting because I'm lazy
    loops = []
    for k in range(1,n):
        loops += list(itertools.permutations(range(n), k))
    for L in loops:
        s = loop_edge_sum(L, times)
        if s < 0:
            return True
    return False

def loop_edge_sum(loop, times):
    # e.g. loop = (a,b,c), want ab+bc+ca
    s = 0
    L = len(loop)
    for k in range(L):
        i = loop[k]
        j = loop[(k+1)%L]
        s += times[i][j]
    return s

def schedule_time(schedule, distances):
    # e.g. schedule = [a,b,c], want ab+bc+ca
    s = 0
    L = len(schedule)
    for k in range(L-1):
        i = schedule[k]
        j = schedule[k+1]
        s += distances[i][j]
    return s

def floyd_warshall(G):
    """
    stolen from https://www.programiz.com/dsa/floyd-warshall-algorithm
    """
    nV = len(G)
    distance = list(map(lambda i: list(map(lambda j: j, i)), G))
    # Adding vertices individually
    for k in range(nV):
        for i in range(nV):
            for j in range(nV):
                distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j])
    return distance

def solution(times, times_limit):
    """
    https://en.wikipedia.org/wiki/Shortest_path_problem#Directed_graph

    1)  Check for timelike loops (where all bunnies can be saved),
        i.e. reduce to "no negative cycles" in the reference terminology.
    2)  Find min time between all pairs (i,j), > -inf after 1).
        https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
    3)  Go through all visitation schedules start-b[0]-b[1]-...-b[k-1]-end in order,
        longest to shortest, bunny-lexicographic in each length.
    4)  Exit at first instance under times_limit
    5)  Probably some edge case I missed, no bunnies, etc.
    """
    N = len(times)
    N_b = N - 2 # number of bunnies
    #print("number of bunnies", N_b)
    # 1) if a loop of negative edge-sum exists, all bunnies can be saved
    if has_timelike_loops(times):
        return [i for i in range(N_b)]
    # 2) find all-pairs shortest paths
    all_pairs_dist = floyd_warshall(times) # stolen code
    #print("all pairs distances")
    #print(all_pairs_dist)
    # 3), 4) piece together all visitiation schedules and exit on the best
    for k in range(N_b,-1,-1):
        #print(k)
        for bunnies in itertools.permutations(range(1,N_b+1), k):
            #print(bunnies)
            schedule = [0] + list(bunnies) + [N-1] # add start/finish
            #print(schedule)
            sched_time = schedule_time(schedule, all_pairs_dist)
            #print(sched_time)
            if times_limit - sched_time >= 0:
                ret_list = schedule[1:-1] # remove start/finish
                ret_list = [r-1 for r in ret_list] # reindex
                ret_list.sort() # sort
                return ret_list # get the hell out
    return "FAIL"

    # 5) some boundary cases I missed???



