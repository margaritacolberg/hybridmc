#import asyncio
import csv
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

def create_bit_flip_graph_with_edge_properties_async(n, transition_info):
    graph = {}
    for i in range(2**n):
        node = format(i, '0' + str(n) + 'b')
        graph[node] = []
        for j in range(n):
            neighbor = list(node)
            neighbor[j] = '1' if neighbor[j] == '0' else '0'
            neighbor = ''.join(neighbor)

            for transition in transition_info:
                if transition[0] == node and transition[1] == neighbor:
                    entropy = transition[2]
                    variance = transition[3] ** 2
                    graph[node].append((neighbor, {'entropy': entropy, 'variance': variance}))



    return graph

def dfs_async(graph, start, end, maxDepth=None, visited=None, path=None):
    if visited is None:
        visited = set()
    if path is None:
        path = []
    visited.add(start)
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    if maxDepth and len(path) > maxDepth:
        return []

    paths = []
    for neighbor, _ in graph[start]:
        if neighbor not in visited and (start, neighbor) not in path:
            new_paths = dfs_async(graph, neighbor, end, maxDepth, visited.copy(), path.copy())
            for new_path in new_paths:
                paths.append(new_path)
    return paths

def find_all_simple_paths_async(graph, start_node, end_node,maxDepth):
    return dfs_async(graph, start_node, end_node, maxDepth)

def sum_path_weights_async(graph, path):
    total_variance = 0
    for i in range(len(path) - 1):
        current_node = path[i]
        next_node = path[i + 1]
        neighbors = graph[current_node]
        for neighbor, properties in neighbors:
            if neighbor == next_node:
                total_variance += properties['variance']
                break
    return total_variance
def sum_path_properties(graph, path):
    total_variance = 0
    total_entropy = 0
    for i in range(len(path) - 1):
        current_node = path[i]
        next_node = path[i + 1]
        neighbors = graph[current_node]
        for neighbor, properties in neighbors:
            if neighbor == next_node:
                total_variance += properties['variance']
                total_entropy += properties['entropy']
                break
    return total_variance, total_entropy

# Example usage:
def process_paths(graph, simple_paths):
    path_entropy = []
    path_weight = []
    for path in simple_paths:
        summed_variance, summed_entropy  = sum_path_properties(graph, path)
        path_entropy.append(summed_entropy)
        path_weight.append(1.0/summed_variance)
        print("Path:", path, ' has variance ', summed_variance, ' and entropy ', summed_entropy, ' with weight',
              1.0 / summed_variance)

    return path_entropy, path_weight

def print_paths(graph, simple_paths):
    for path in simple_paths:
        summed_variance, summed_entropy  = sum_path_properties(graph, path)
        print("Path:", path, ' has variance ', summed_variance, ' and entropy ', summed_entropy, ' with weight', 1.0/summed_variance)
def get_transition_info(diff_s_bias_csv):
    transition_info = []  # Load transition info from CSV file
    # Simulated transition info, replace this with actual loading from CSV
    with open(diff_s_bias_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            transition_info.append((row[0], row[1], float(row[2]), float(row[3])))
            transition_info.append((row[1], row[0], -float(row[2]), float(row[3])))
    return transition_info



def main():
    n = 4
    transition_info = get_transition_info('diff_s_bias_with_error.csv')
    graph = create_bit_flip_graph_with_edge_properties_async(n, transition_info)

    print('Graph is:\n', graph)

    start_node = '1100'
    end_node = '1101'
    maxDepth = n + 6
    simple_paths = find_all_simple_paths_async(graph, start_node, end_node, maxDepth)
    sorted_paths = sorted(simple_paths, key=lambda t: (len(t)))

    path_entropy, path_weight = process_paths(graph, sorted_paths)
    print('There are a total of ', len(simple_paths), ' paths connecting ', start_node, ' to ', end_node)

    d1 = DescrStatsW(path_entropy, weights=path_weight)

    # find the standard error of the weighted average by taking 1 / sqrt(the sum of the weights for each path)
    std_error = 1 / np.sqrt(d1.sum_weights) if d1.sum_weights else 0
    print('  Sbias over path = ', d1.mean, ' with standard error ', std_error)


main()