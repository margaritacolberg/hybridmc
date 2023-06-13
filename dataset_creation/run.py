import os
import subprocess
import sys
import argparse

sys.path.append(os.getcwd())

parser = argparse.ArgumentParser()
parser.add_argument('nbonds', help='max number of bonds in structure')
parser.add_argument('json_in', help='Input json file containing the structural information')
parser.add_argument('tag', help="Version number or the trial number")
args = parser.parse_args()

# Function to get file paths
def get_path(file_name, search_path=os.getenv('HOME')):
    """
    Find the path of file_name named file.
    Parameters
    ----------
    file_name (str): The file that needs its path found
    search_path (str). The path in which to find file_name

    Returns
    -------
    THe string contianing the path of file_name
    """
    for root, _, files in os.walk(search_path, topdown=0):
        if file_name in files:
            return os.path.join(root, file_name)


def main():

    # Find the SLURM scheduler details for this job

    # SLURM_ARRAY_JOB_ID = os.getenv('SLURM_ARRAY_JOB_ID', default=1)
    SLURM_ARRAY_TASK_ID = os.getenv('SLURM_ARRAY_TASK_ID', default=1)
    nproc = os.getenv('SLURM_CPUS_PER_TASK', default=os.cpu_count())

    # Find the paths for the source for running the simulation
    file_name = get_path('init_json.py')
    exe = get_path('hybridmc')

    # Specify descriptive names for directory where the simualation will output results
    dir_name = f'train_{args.nbonds}_{args.tag}_{SLURM_ARRAY_TASK_ID}'
    tmp_dir_name = f'{dir_name}.tmp'

    # Create said directory if it does not already exist
    if not os.path.isdir(dir_name):
        if not os.path.isdir(tmp_dir_name):
            os.mkdir(tmp_dir_name)

        # Run the simulation for the specified structure (based on the input json file)
        subprocess.run(['python', file_name, args.json_in, exe, '1', '--nproc={}'.format(nproc)], cwd=tmp_dir_name, check=1)

        # Once done rename temporary directory
        os.rename(src=tmp_dir_name, dst=dir_name)


if __name__ == '__main__':
    main()
