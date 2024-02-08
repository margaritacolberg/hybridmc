from . import ConfigGenerator, ConfigGeneratorDriver


def main():
    config_generator = ConfigGeneratorDriver(settings_config_file="settings.cfg")
    config_generator.generate_configs()


# Example usage
if __name__ == "__main__":
    main()
