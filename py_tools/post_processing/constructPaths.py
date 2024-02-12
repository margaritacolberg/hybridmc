# import asyncio
import csv
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW


def create_bit_flip_graph_with_edge_properties_async(n, transition_info):
    graph = {}
    for i in range(2 ** n):
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


def find_all_simple_paths_async(graph, start_node, end_node, maxDepth):
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
        summed_variance, summed_entropy = sum_path_properties(graph, path)
        path_entropy.append(summed_entropy)
        path_weight.append(1.0 / summed_variance)
        lower_confidence = summed_entropy - 1.96 * np.sqrt(summed_variance)
        upper_confidence = summed_entropy + 1.96 * np.sqrt(summed_variance)
        print("Path:", path, ' has entropy ', summed_entropy, ' with CI (', lower_confidence, ",",
              upper_confidence, ') with weight ', 1.0 / summed_variance)

    return path_entropy, path_weight


def print_paths(graph, simple_paths):
    for path in simple_paths:
        summed_variance, summed_entropy = sum_path_properties(graph, path)

        print("Path:", path, ' has confidence interval (', lower_confidence, ",",
              upper_confidence, ') and entropy ',
              summed_entropy, ' with weight', 1.0 / summed_variance)


def get_transition_info(diff_s_bias_csv):
    transition_info = []  # Load transition info from CSV file
    # Simulated transition info, replace this with actual loading from CSV
    with open(diff_s_bias_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            transition_info.append((row[0], row[1], float(row[2]), float(row[3])))
            transition_info.append((row[1], row[0], -float(row[2]), float(row[3])))
    return transition_info


def main(start_node='1101', end_node='1111', depth=2):
    n = len(start_node)

    if n != len(end_node):
        raise ValueError("lengths of start and end node have to be equal")

    transition_info = get_transition_info('diff_s_bias_with_error.csv')
    graph = create_bit_flip_graph_with_edge_properties_async(n, transition_info)

    print('Graph is:\n', graph)

    maxDepth = n + depth
    simple_paths = find_all_simple_paths_async(graph, start_node, end_node, maxDepth)
    sorted_paths = sorted(simple_paths, key=lambda t: (len(t)))

    path_entropy, path_weight = process_paths(graph, sorted_paths)
    print('There are a total of ', len(simple_paths), ' paths connecting ', start_node, ' to ', end_node)

    d1 = DescrStatsW(path_entropy, weights=path_weight)

    # find the standard error of the weighted average by taking 1 / sqrt(the sum of the weights for each path)
    std_error = 1 / np.sqrt(d1.sum_weights) if d1.sum_weights else 0

    lower_ci = d1.mean - 1.96 * std_error
    upper_ci = d1.mean + 1.96 * std_error

    print('  Sbias over path = ', d1.mean, ' with 95% confidence interval (', lower_ci,
          ',', upper_ci, ')')

    mean1 = d1.mean
    SE1 = std_error
    mean2 = path_entropy[0]
    SE2 = 1 / np.sqrt(path_weight[0])

    return diff_test(mean1, mean2, SE1, SE2)


def diff_test(mean1, mean2, SE1, SE2, Z=2.57):
    mean_diff = abs(mean1 - mean2)
    SE = np.sqrt(SE1 ** 2 + SE2 ** 2)

    conf_int = mean_diff - SE * Z, mean_diff + SE * Z

    print(conf_int)

    if conf_int[0] < 0 < conf_int[1]:
        print("There is no significant difference")
        return mean1, mean2, False, conf_int
    else:
        print("There is a significant difference")
        return mean1, mean2, True, conf_int


if __name__ == "__main__":
    transition_info = get_transition_info('diff_s_bias_with_error.csv')

    src_dst_nodes = list(zip(*list(zip(*transition_info))[:2]))

    field_names = ["src_node", "dst_node", "path mean", "direct mean","Sig_diff?", "diff_conf_int"]
    output = []
    n_trues = 0
    for src, dst in src_dst_nodes:

        mean1, mean2, significant_difference_test_result, diff_conf_int = main(src, dst, depth=4)
        output.append([src, dst, str(mean1), str(mean2), str(significant_difference_test_result), str(diff_conf_int)])
        if significant_difference_test_result:
            n_trues += 1

    print(f"The failure percentage was: {100 * n_trues / len(output)}")

    # Write data to CSV file
    with open("sig_diffs.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write the header
        writer.writerow(field_names)
        # Write the data rows
        writer.writerows(output)


