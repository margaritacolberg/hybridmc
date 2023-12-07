from functools import wraps
import time
import run
import argparse


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


@timeit
def hmc(args):
    """
    Simple function that returns sum of all numbers up to the square of num.
    """
    return run.main(args)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='simple_stair_0.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../release/hybridmc")
    parser.add_argument('-ov', '--old_version', help='set version code for old structure simulation run if needed',
                        default='old_2')

    args = parser.parse_args()
    hmc(args)