import networkx as nx
import pandas as pd
import os
from tqdm import tqdm
import multiprocessing

# Suppose you have already defined:
#   geneId_geneName
#   geneName_geneId
# above or in some shared module.
# Make sure any global dict is accessible from within the child process.
# loading and cleaning ENSG converter
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

def process_data_multiple_endpoints(task):
    """
    For a single 'origin', run single_source_dijkstra for all endpoints.
    Return a dict counting how many times each node appears
    on any shortest path from that origin to each endpoint.
    """
    G = task['graph']
    origin = task['origin']
    endpoints = task['endpoints']
    
    node_count = {}

    try:
        # single_source_dijkstra -> (distances, {node: path})
        dist_dict, path_dict = nx.single_source_dijkstra(G, source=origin, weight='weight')
        
        for endpoint in endpoints:
            path = path_dict.get(endpoint)
            if path:
                # Exclude origin (path[0]) and endpoint (path[-1]) from the count
                for node in path[1:-1]:
                    node_count[node] = node_count.get(node, 0) + 1

    except nx.NetworkXNoPath:
        # If no path is found for this origin, ignore
        pass
    except Exception as e:
        # Log other unexpected errors to help debug
        print(f"[Worker error] origin={origin}, error={e}")
    
    return node_count


def create_regulatory(regulatory_filepath):
    """
    Example: read regulatory file, map row/column labels to ENSG IDs,
    drop anything unmappable, then transform (1/x, abs).
    """
    print(f"Loading: {regulatory_filepath}")
    df = pd.read_csv(regulatory_filepath, index_col=0)

    # 1) Map row names from geneName -> ENSG if possible
    new_index = []
    drop_rows = []
    for old_name in df.index:
        upper_name = old_name.upper()
        if upper_name in geneName_geneId:
            new_index.append(geneName_geneId[upper_name])
        else:
            # no valid mapping -> mark for drop
            new_index.append(None)
            drop_rows.append(old_name)
    
    # Drop unmappable rows
    df['new_index'] = new_index
    df.drop(drop_rows, inplace=True)
    df.set_index('new_index', inplace=True)
    
    # 2) Map column names from geneName -> ENSG if needed
    #    (depends on your data; example below might differ from your real code)
    keep_cols = []
    new_cols = []
    for col in df.columns:
        upper_col = col.upper()
        if upper_col in geneName_geneId:
            keep_cols.append(col)
            new_cols.append(geneName_geneId[upper_col])
        else:
            # Not mapped
            pass
    
    # Keep only the mapped columns
    df = df[keep_cols]
    df.columns = new_cols

    # 3) Transform each cell: absolute(1/x)
    #    (assuming no zero values — watch for division by zero)
    df = df.applymap(lambda x: abs(1/x))

    return df


def create_graph(regulatory_df, endpoints):
    """
    Build an undirected Graph from the regulatory dataframe + endpoint nodes
    (ensuring no Nones end up in the graph).
    """
    G = nx.Graph()

    # Collect all unique nodes
    row_nodes = set(regulatory_df.index)
    col_nodes = set(regulatory_df.columns)
    endpoint_nodes = set(endpoints)  # Make sure 'endpoints' also has no None
    
    # Remove None from the endpoint list, in case it sneaks in
    endpoint_nodes.discard(None)

    # Add nodes in one go
    all_nodes = row_nodes.union(col_nodes).union(endpoint_nodes)
    G.add_nodes_from(all_nodes)

    # Add weighted edges from the regulatory matrix
    for row_name in regulatory_df.index:
        row_series = regulatory_df.loc[row_name]
        for col_name, val in row_series.items():
            G.add_edge(row_name, col_name, weight=val)

    return G


def create_endpoints(filepath):
    """
    Load endpoints from a TSV file, map gene names to ENSG IDs,
    drop unmappable ones, and return them as a list.
    """
    df = pd.read_csv(filepath, sep='\t', index_col=0)

    # Map index -> new_index if possible
    new_index = []
    drop_rows = []
    for old_name in df.index:
        upper_name = old_name.upper()
        if upper_name in geneName_geneId:
            new_index.append(geneName_geneId[upper_name])
        else:
            new_index.append(None)
            drop_rows.append(old_name)

    # Attach new_index, then drop rows that didn't map
    df['new_index'] = new_index
    df.drop(drop_rows, inplace=True)
    df.set_index('new_index', inplace=True)

    # The final index is ENSG IDs. Return as a list:
    endpoints = list(df.index)

    # Double-check no None in the list:
    endpoints = [e for e in endpoints if e is not None]
    return endpoints


def connection_enrichment(origins, endpoints):
    """
    Build a graph for each dataset, run single_source_dijkstra for each origin
    in parallel, combine the node counts.
    """
    data_dir = 'data'
    datasets = [f for f in os.listdir(data_dir) if f.endswith('.csv') or f.endswith('.tsv')]
    if '.DS_Store' in datasets:
        datasets.remove('.DS_Store')

    total_count = {}

    for i, dataset in enumerate(datasets):
        print(f"\nDataset {i}/{len(datasets)}: {dataset}")
        regulatory_df = create_regulatory(os.path.join(data_dir, dataset))
        G = create_graph(regulatory_df, endpoints)

        # Prepare tasks for each origin
        tasks = [{'graph': G, 'origin': o, 'endpoints': endpoints} for o in origins]

        # Parallel execution
        with multiprocessing.Pool(processes=16) as pool:
            results_iter = pool.imap(process_data_multiple_endpoints, tasks)

            # Consume results
            for node_count_dict in tqdm(results_iter, total=len(tasks)):
                for node, cnt in node_count_dict.items():
                    total_count[node] = total_count.get(node, 0) + cnt

    return total_count

if __name__ == "__main__":
    # Top-level code that calls the functions
    globalAgingGenes = create_endpoints('global_aging_genes.tsv')
    
    # This will now run the parallel code properly
    agingCount = connection_enrichment(globalAgingGenes, globalAgingGenes)
    
    # Print, save, or analyze `agingCount`
    print("Done. Number of nodes with counts:", len(agingCount))