def solution_clear(map):
    """
    map column-major {0,1}-matrix

    Following wikipedia's Dijkstra's alg description,
    https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

    obviously need to modify for wall removal?
    """
    #import math???
    # math.inf not supported until python 3.5, so inf > 200^2
    inf = 1000000
    W = len(map)
    H = len(map[0])
    #print(W, H)
    # tentative distances to (0,0)
    tent_dists = [[inf for i in range(H)] for j in range(W)]
    tent_dists[0][0] = 1 # initial node at distance 1

    # unvisited nodes only including empty/zero nodes
    unvisited_set = {(i,j) for i in range(W) for j in range(H) if map[i][j] == 0}

    # return unvisited neighbors of node
    def unvisited_neighbors(node, unvisited_set):
        ret_list = []
        for i in [1, -1]:
            pot_neigh = (node[0]+i, node[1])
            if pot_neigh in unvisited_set:
                ret_list.append(pot_neigh)
            pot_neigh = (node[0], node[1]+i)
            if pot_neigh in unvisited_set:
                ret_list.append(pot_neigh)
        return ret_list

    #set current node
    current_node = (0,0)

    while True:
        # look at unvisited neighbors
        neighbors = unvisited_neighbors(current_node, unvisited_set)
        #print(neighbors)
        for N in neighbors:
            old = tent_dists[N[0]][N[1]]
            new = 1 + tent_dists[current_node[0]][current_node[1]]
            tent_dists[N[0]][N[1]] = min([old, new])

        # remove current node from unvisited
        # and set new smallest current node or exit
        unvisited_set.remove(current_node)
        if (W-1, H-1) not in unvisited_set:
            # already went through terminal node
            break
        else:
            # choose an unvisited node with minimal distance to source
            current_node = min(unvisited_set, key=lambda x : tent_dists[x[0]][x[1]])
        if tent_dists[current_node[0]][current_node[1]] == inf:
            #print("No path!")
            break
    return tent_dists[W-1][H-1]

import copy

def solution(map):
    """
    inefficiently consider all single-wall variations...
    """
    W = len(map)
    H = len(map[0])
    lengths = []
    lengths.append(solution_clear(map))
    for node in [(i,j) for i in range(W) for j in range(H) if map[i][j] == 1]:
        new_map = copy.deepcopy(map)
        new_map[node[0]][node[1]] = 0
        lengths.append(solution_clear(new_map))
    return min(lengths)
