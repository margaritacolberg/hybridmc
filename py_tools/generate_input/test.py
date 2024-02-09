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

    # submit job to slurm for these configs
    jobsubmitter = JobSubmitter(
        json_dir=config.get('master_settings', 'target_directory', fallback="generated_configs"),
        out_dir=config.get('slurm_settings', 'out_dir', fallback="config_run"),
        Nconfigs=config.get('slurm_settings', 'job_arrays', fallback="1-10"),
        cpus_per_task=2, mem_per_cpu=1000, time='0-1:55:00',
        exe="/scratch/vignesh9/hybridmc/py_bin/run.py",
        hmc_exe="/scratch/vignesh9/hybridmc/release/hybridmc",
        abspath=1
    )

    #jobsubmitter.create_job_script()
    jobsubmitter.submit_job()
    #jobsubmitter.exec_script()
    print("DONE SUBMITTING")


# Example usage
if __name__ == "__main__":
    main()
