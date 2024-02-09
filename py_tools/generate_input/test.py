from . import ConfigGenerator, ConfigGeneratorDriver, JobSubmitter


def main():
    #config_generator = ConfigGeneratorDriver(settings_config_file="settings.cfg")
    #config_generator.generate_configs()

    jobsubmitter = JobSubmitter()
    jobsubmitter.exe = '../../../py_bin/run.py'

    jobsubmitter.create_job_script()
    print(1)


# Example usage
if __name__ == "__main__":
    main()
