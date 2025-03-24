
import networkx as nx

def process_data(data):
    G = data['graph']
    origin = data['origin']
    endpoint = data['endpoint']
    count = {}
    try:
        path = nx.shortest_path(G, origin, endpoint, weight="weight")
        path.pop(0)
        path.pop()
        for node in path:
            if node not in count:
                count[node] = 1
            else:
                count[node] += 1
    except:
        pass
        
    return count
