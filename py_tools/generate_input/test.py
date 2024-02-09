import os
import sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from py_tools.generate_input import ConfigGenerator, ConfigGeneratorDriver, JobSubmitter
import configparser


def main():
    # Generate some configurations
    # config_generator = ConfigGeneratorDriver(settings_config_file="settings.cfg")
    # config_generator.generate_configs()

    # read config settings
    config = configparser.ConfigParser()
    config.read("settings.cfg")

    Nconfigs = config.get('master_settings', 'job_arrays', fallback="1-10")

    # submit job to slurm for these configs
    jobsubmitter = JobSubmitter(
        target_dir=config.get('master_settings', 'target_directory', fallback="generated_configs"),
        Nconfigs=Nconfigs,
        cpus_per_task=4, mem_per_cpu=500, time='0-01:00:00',
        exe="/scratch/vignesh9/hybridmc/py_bin/run.py"
    )

    print(jobsubmitter.exe)
    jobsubmitter.submit_job()
    print("DONE SUBMITTING")


# Example usage
if __name__ == "__main__":
    main()
