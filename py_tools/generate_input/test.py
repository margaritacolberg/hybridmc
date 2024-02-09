from . import ConfigGenerator, ConfigGeneratorDriver, JobSubmitter
import configparser


def main():

    # Generate some configurations
    config_generator = ConfigGeneratorDriver(settings_config_file="settings.cfg")
    config_generator.generate_configs()

    # read config settings
    config = configparser.ConfigParser()
    config.read("settings.cfg")

    # submit job to slurm for these configs
    jobsubmitter = JobSubmitter(target_dir=config.get('master_settings', 'target_directory', fallback="generated_configs"))
    jobsubmitter.submit_job()
    print("DONE SUBMITTING")


# Example usage
if __name__ == "__main__":
    main()
