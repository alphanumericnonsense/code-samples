def augment(entrances, exits, paths):
    """
    augment the capacity matrix with a supersource and supersink with capacities the sum of the capacities out of the sources/sinks (although any big enough number shoudl work).
    """
    N = len(paths)
    # sum of capacities out of sources and into sinks
    smax = sum([paths[i][j] for i in entrances for j in range(N)])
    tmax = sum([paths[i][j] for i in range(N) for j in exits])
    # augmented capacities
    C = []
    # from supersource to sources
    sourcerow = [0] + [smax if i in entrances else 0 for i in range(N)] + [0]
    C.append(sourcerow)
    # stuff inbetween, need stuff out of sinks to supersink
    for i in range(N):
        if i in exits:
            midrow = [0] + paths[i] + [tmax]
        else:
            midrow = [0] + paths[i] + [0]
        C.append(midrow)
    # nothing out of supersink
    sinkrow = [0 for i in range(N+2)]
    C.append(sinkrow)
    return C

import collections

class Graph:
    '''
    This class represents a directed graph using
    adjacency matrix representation.
    '''

    def __init__(self, graph):
        self.graph = graph  # residual graph
        self.row = len(graph)

    def bfs(self, s, t, parent):
        '''
        Returns true if there is a path from
        source 's' to sink 't' in residual graph.
        Also fills parent[] to store the path.
        '''

        # Mark all the vertices as not visited
        visited = [False] * self.row

        # Create a queue for BFS
        queue = collections.deque()

        # Mark the source node as visited and enqueue it
        queue.append(s)
        visited[s] = True

        # Standard BFS loop
        while queue:
            u = queue.popleft()

            # Get all adjacent vertices of the dequeued vertex u
            # If an adjacent has not been visited, then mark it
            # visited and enqueue it
            for ind, val in enumerate(self.graph[u]):
                if (visited[ind] == False) and (val > 0):
                    queue.append(ind)
                    visited[ind] = True
                    parent[ind] = u

        # If we reached sink in BFS starting from source, then return
        # true, else false
        return visited[t]

    # Returns the maximum flow from s to t in the given graph
    def edmonds_karp(self, source, sink):

        # This array is filled by BFS and to store path
        parent = [-1] * self.row

        max_flow = 0  # There is no flow initially

        # Augment the flow while there is path from source to sink
        while self.bfs(source, sink, parent):

            # Find minimum residual capacity of the edges along the
            # path filled by BFS. Or we can say find the maximum flow
            # through the path found.
            path_flow = float("Inf")
            s = sink
            while s != source:
                path_flow = min(path_flow, self.graph[parent[s]][s])
                s = parent[s]

            # Add path flow to overall flow
            max_flow += path_flow

            # update residual capacities of the edges and reverse edges
            # along the path
            v = sink
            while v != source:
                u = parent[v]
                self.graph[u][v] -= path_flow
                self.graph[v][u] += path_flow
                v = parent[v]

        return max_flow

def test_cases():
    paths1 = [[0, 0, 4, 6, 0, 0], [0, 0, 5, 2, 0, 0], [0, 0, 0, 0, 4, 4], [0, 0, 0, 0, 6, 6], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    entrances1 = [0, 1]
    exits1 = [4, 5]
    C1 = augment(entrances1, exits1, paths1)
    X1 = Graph(C1)
    print(X1.edmonds_karp(0, len(C1)-1))

    paths2 = [[0, 7, 0, 0], [0, 0, 6, 0], [0, 0, 0, 8], [9, 0, 0, 0]]
    entrances2 = [0]
    exits2 = [3]
    C2 = augment(entrances2, exits2, paths2)
    X2 = Graph(C2)
    print(X2.edmonds_karp(0, len(C2)-1))
