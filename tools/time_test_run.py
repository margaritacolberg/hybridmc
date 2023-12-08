from functools import wraps
import time
import run
import argparse
import json
import csv

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return total_time

    return timeit_wrapper

def get_json_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def update_json_data(json_data, new_values):
    for key, value in new_values.items():
        json_data[key] = value
    return json_data


def save_json_data(json_data, file_path):
    with open(file_path, 'w') as f:
        json.dump(json_data, f, indent=4)


def test_runtimes(file_path, new_values):
    data = get_json_data(file_path)
    data = update_json_data(data, new_values)
    save_json_data(data, 'runtime_test.json')

    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='runtime_test.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../release/hybridmc")
    parser.add_argument('-ov', '--old_version', help='set version code for old structure simulation run if needed',
                        default='old_2')

    args = parser.parse_args()

    return hmc(args)


@timeit
def hmc(args):
    """
    Simple function that runs main
    """
    return run.main(args)


if __name__ == '__main__':
    nbeads = [30, 40]
    del_t = [1, 2, 3, 4, 5]

    times = [[], [], []]

    for nbead in nbeads:
        for dt in del_t:
            time = test_runtimes('time_test.json', {'nbeads': nbead, 'del_t': dt})
            times[0].append(nbead), times[1].append(dt), times[2].append(time)

    with open('time_test.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['nbeads', 'del_t', 'time'])
        writer.writerows(zip(*times))




