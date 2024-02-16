import os
import sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from py_tools.generate_input import ConfigGenerator, ConfigGeneratorDriver, JobSubmitter
import configparser


def main(args):

    if args.gen:
        # Generate some configurations
        config_generator = ConfigGeneratorDriver(settings_config_file="settings.cfg")
        config_generator.generate_configs()

    # read config settings
    config = configparser.ConfigParser()
    config.read("settings.cfg")

    # submit job to slurm for these configs
    jobsubmitter = JobSubmitter(
        json_dir=config.get('master_settings', 'target_directory', fallback="generated_configs"),
        out_dir=config.get('slurm_settings', 'out_dir', fallback="config_run"),
        Nconfigs=config.get('slurm_settings', 'job_arrays', fallback="1-10"),
        cpus_per_task=2, mem_per_cpu=1000, time='0-5:55:00',
        exe="/scratch/vignesh9/hybridmc/py_bin/run.py",
        hmc_exe="/scratch/vignesh9/hybridmc/release/hybridmc",
    )

    #jobsubmitter.create_job_script()
    jobsubmitter.submit_job()
    #jobsubmitter.exec_script()
    print("DONE SUBMITTING")


# Example usage
if __name__ == "__main__":
    import argparse

    # Create the parser
    parser = argparse.ArgumentParser(description="Example script showing how to use argparse for a bool and a value.")

    # Add a boolean flag (default=False). If the flag is used, the value is set to True.
    parser.add_argument('-g', '--gen', action='store_true', help="Generate structures or no")

    # Parse the command-line arguments
    args = parser.parse_args()

    main(args)
