
import pandas as pd
import os
from datetime import datetime
from tqdm import tqdm
import statistics
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing
from graph_tool import topology, Graph

# loading and cleaning ENSG converter
def createGeneConverter():
    test = []
    geneId_geneName = {}
    with open('Homo_sapiens.GRCh37.74.gtf', 'r') as file:
        for line in file:
            line = line.strip()
            data = line.split('\t')[-1]
            test.append(data)
            if 'gene_name' in data:
                attributes = data.split(';')
                geneId = attributes[0].split(' ')[1].strip('"')
                for attr in attributes:
                    if 'gene_name' in attr:
                        geneName = attr.split(' ')[2].strip('"')
                        if geneId not in geneId_geneName:
                            geneId_geneName[geneId] = geneName
    geneName_geneId = {v: k for k, v in geneId_geneName.items()}
    return geneName_geneId, geneId_geneName

def createRegulatory(regulatory_filepath):
    geneName_geneId, geneId_geneName = createGeneConverter()
    # loading, cleaning, and permutating regulatory dataset
    print(f'Loading: {regulatory_filepath}')
    regulatory = pd.read_csv(regulatory_filepath, index_col=0)

    exceptions = []
    for name in regulatory.index:
        try:
            regulatory = regulatory.rename(index={name: geneName_geneId[name]})
        except:
            exceptions.append(name)

    print(f'Row exception count: {len(exceptions)}')

    for exc in exceptions:
        regulatory = regulatory.drop(exc)

    # if columns are ENSG IDs, uncomment this part
    exceptions = []
    for ID in regulatory.columns.tolist():
        if ID not in geneId_geneName:
            exceptions.append(ID)

    # if columns are gene names, uncomment this part
    # exceptions = []
    # for name in regulatory.columns.tolist():
    #     try:
    #         regulatory = regulatory.rename(columns={name: geneName_geneId[name]})
    #     except:
    #         exceptions.append(name)


    print(f'Column exception count: {len(exceptions)}')

    for exc in exceptions:
        regulatory = regulatory.drop(exc)


    def inverse(x):
        return 1/x

    def absolute(x):
        return abs(x)

    regulatory = regulatory.map(inverse)
    regulatory = regulatory.map(absolute)

    return regulatory

def createEndpoints(filepath):
    geneName_geneId, geneId_geneName = createGeneConverter()
    # loading and cleaning dataset
    endpoints = pd.read_csv(filepath, sep='\t', index_col=0)

    exceptions = []
    for name in endpoints.index:
        try:
            endpoints = endpoints.rename(index={name: geneName_geneId[name.upper()]})
        except:
            exceptions.append(name)

    print(f'Endpoints exception count: {len(exceptions)}')

    for exc in exceptions:
        endpoints = endpoints.drop(exc)

    return list(endpoints.index)

def createRandomEndpoints(regulatory, num, seed):
    rand = regulatory.sample(n=num, random_state=seed)
    return list(rand.index)

def createGraph(regulatory, endpoints):
    regMatrix = regulatory.to_numpy().tolist()
    nodes = list(set(list(regulatory.index) + regulatory.columns.tolist() + endpoints))
    g = Graph()
    g.set_directed(False)
    name_prop = g.new_vertex_property('string')
    vertexNameDic = {}
    for node in nodes:
        v = g.add_vertex()
        name_prop[v] = node
        vertexNameDic[v] = node
    g.vertex_properties['name'] = name_prop
    nameVertexDic = {v: k for k, v in vertexNameDic.items()}

    edgeCount = 0
    g.ep.weight = g.new_edge_property("double")
    for rowName, row in zip(regulatory.index, regMatrix):
        for columnName, cell in zip(regulatory.columns.tolist(), row):

            rowVertex = nameVertexDic[rowName]
            columnVertex = nameVertexDic[columnName]

            edge = g.add_edge(rowVertex, columnVertex)
            g.ep.weight[edge] = cell
            edgeCount += 1

    print(f'EdgeCount: {edgeCount}')

    return g, nameVertexDic, vertexNameDic

def processData(data):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    dataset = data['dataset']
    origins = data['origins']
    endpoints = data['endpoints']
    regulatory = createRegulatory(f'data/{dataset}')
    g, nameVertexDic, vertexNameDic = createGraph(regulatory, endpoints)

    count = {}
    exceptionCount = 0
    for i, origin in enumerate(origins):
        for j, endpoint in enumerate(endpoints):
            path = topology.shortest_path(g, nameVertexDic[origin], nameVertexDic[endpoint], weights=g.ep.weight)
            path = path[0]
            # print(i, j, len(origins), len(endpoints))
            print(f'origin{i}:{nameVertexDic[origin]} endpoint{j}:{nameVertexDic[endpoint]} path:{path}')
            if len(path) == 0:
                continue
            path.pop(0)
            path.pop()
            for node in [vertexNameDic[node] for node in path]:
                if node not in count:
                    count[node] = 1
                else:
                    count[node] += 1
    return count

def connectionEnrichment(origins, endpoints):
    datasets = os.listdir('data')
    datasets.remove('.DS_Store')
    count = {}

    data = [{'dataset': dataset, 'origins': origins, 'endpoints': endpoints} for dataset in datasets]
    # with multiprocessing.Pool(processes=12) as pool:
    #     results = pool.map(processData, data)
    results = [processData(data[0])]

    for dic in results:
        for key, val in dic.items():
            if key not in count:
                count[key] = 1
            else:
                count[key] += 1
    return count

def sortDic(dic):
    return dict(reversed(sorted(dic.items(), key=lambda item: item[1])))

def printResultsWithStats(dic):
    print(f'Size: {len(dic)}')
    print(f'Average: {statistics.mean(dic.values())}')
    print(f'Median: {statistics.median(dic.values())}')
    print(dic)
