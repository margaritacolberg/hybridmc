import json
import random
import csv


class ConfigGenerator:
    def __init__(self, base_config=None, num_files=5, target_dir=None):
        self.base_config = base_config if base_config else {}
        self.param_settings = {}
        self.num_files = num_files

        self.target_dir = target_dir
        if self.target_dir:
            create_and_enter_dir(directory_path=target_dir)

    def set_base_config(self, base_config):
        self.base_config = base_config

    def set_num_files(self, num_files):
        self.num_files = num_files

    def add_parameter_setting(self, parameter, setting):
        self.param_settings[parameter] = setting

    def generate_configs(self) -> object:
        new_configs = []
        csv_data = [["Filename"] + list(self.param_settings.keys()) + ["Other specified parameters"]]

        for i in range(self.num_files):
            new_config = self.base_config.copy()
            row = [f"config_{i}.json"]
            for param, settings in self.param_settings.items():
                if 'range' in settings:
                    new_value = random.uniform(*settings['range'])
                elif 'values' in settings:
                    new_value = random.choice(settings['values'])
                new_config[param] = new_value
                row.append(new_value)
            new_configs.append(new_config)
            row.append("Other values")  # Add other specified parameters here
            csv_data.append(row)

            # Write the new JSON file
            with open(f"config_{i}.json", 'w') as json_file:
                json.dump(new_config, json_file, indent=4)

        self._write_csv(csv_data)

    def _write_csv(self, csv_data):
        # Write the CSV file
        with open("configurations.csv", 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows(csv_data)

# some helper functions

import os
def create_and_enter_dir(directory_path):
    # Check if the directory already exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create it
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' was created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

    # Change the current working directory to the specified directory
    os.chdir(directory_path)
    print(f"Current working directory: {os.getcwd()}")


# Example usage
directory_path = 'path/to/your/directory'
create_and_enter_dir(directory_path)


import json
import random
import csv
import os
from typing import Any, Dict

class ConfigGenerator:
    target_directory = 'generated_configs'  # Default target directory

    def __init__(self, base_config: Dict[str, Any] = None, num_files: int = 5, csv_file: str = "configurations.csv"):
        self.base_config = base_config if base_config else {}
        self.param_settings = {}
        self.num_files = num_files
        self.csv_file = os.path.join(ConfigGenerator.target_directory, csv_file)
        self.configs = []  # To store generated configurations for uniqueness check
        self.id_counter = 1  # Start from 1
        self._ensure_directory_exists(ConfigGenerator.target_directory)

    @staticmethod
    def _ensure_directory_exists(directory: str) -> None:
        """Ensure the target directory exists, creating it if necessary."""
        os.makedirs(directory, exist_ok=True)

    def add_parameter_setting(self, parameter: str, setting: Dict[str, Any]) -> None:
        self.param_settings[parameter] = setting

    def _is_unique_config(self, config: Dict[str, Any]) -> bool:
        """Check if a generated configuration is unique."""
        return config not in self.configs

    def _update_csv(self, csv_data: list) -> None:
        """Update the CSV file with new configurations."""
        file_exists = os.path.exists(self.csv_file)
        with open(self.csv_file, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if not file_exists:
                # Write headers only if the file did not exist
                writer.writerow(["ID"] + list(self.param_settings.keys()) + ["Filename"])
            writer.writerows(csv_data)

    def generate_configs(self) -> None:
        csv_data = []

        for _ in range(self.num_files):
            while True:
                new_config = self.base_config.copy()
                row = [self.id_counter]
                for param, settings in self.param_settings.items():
                    if 'range' in settings:
                        new_value = random.uniform(*settings['range'])
                    elif 'values' in settings:
                        new_value = random.choice(settings['values'])
                    new_config[param] = new_value
                    row.append(new_value)
                if self._is_unique_config(new_config):
                    break

            self.configs.append(new_config)
            filename = os.path.join(ConfigGenerator.target_directory, f"config_{self.id_counter}.json")
            row.append(filename)
            csv_data.append(row)

            with open(filename, 'w') as json_file:
                json.dump(new_config, json_file, indent=4)

            self.id_counter += 1

        self._update_csv(csv_data)


# Example usage
ConfigGenerator.set_target_directory('path/to/your/directory')  # Set custom directory
generator = ConfigGenerator()
# Set base_config, add_parameter_setting, etc.
# generator.generate_configs()
