import configparser
import json
import random
import csv
import os
import warnings
from typing import Any, Dict, List


class ConfigGenerator:

    def __init__(self, base_config: Dict[str, Any] = None, num_files: int = 5,
                 target_directory: str = "Generated_configs", csv_file_name: str = "configs.csv"):

        self.base_config = base_config
        self.num_files = num_files
        self._configs = []  # To store generated configurations for uniqueness check
        self.id_counter = None
        self.param_settings = {}
        self.csv_file_name = csv_file_name

        # set target directory, also ensure all paths like the csv name are changed accordingly
        self.target_directory = target_directory

    @property
    def id_counter(self):
        """Property for dynamically getting the next ID counter value."""
        if self._id_counter is None:
            self._initialize_id_counter()

        return self._id_counter

    @id_counter.setter
    def id_counter(self, value):
        """Allow manual setting of the ID counter."""
        self._id_counter = value

    def _initialize_id_counter(self):
        """Initialize ID counter by finding the last ID used in the CSV file."""
        try:
            with open(self.csv_file_name, 'r') as csv_file:
                reader = csv.reader(csv_file)
                last_row = list(reader)[-1]  # Go to the last row
                last_id = int(last_row[0])  # Assuming the ID is the first column
        except (FileNotFoundError, IndexError):
            last_id = 0  # Start from 1 if no file exists or file is empty

        self._id_counter = last_id + 1

    @property
    def base_config(self) -> Dict[str, Any]:
        return self._base_config

    @base_config.setter
    def base_config(self, value: Dict[str, Any]) -> None:

        if type(value) == str:

            with open(value, 'r') as file:
                base_config = json.load(file)
            self._base_config = base_config

        elif type(value) == dict:
            self._base_config = value

        elif value is None:
            self._base_config = {}

        else:
            raise ValueError("Need a Dict or json string. Defaults to empty dictionary")

    @property
    def target_directory(self) -> str:
        return self._target_directory

    @target_directory.setter
    def target_directory(self, value: str) -> None:
        self._target_directory = value
        os.makedirs(value, exist_ok=True)
        self.csv_file_name = os.path.join(self.target_directory, self.csv_file_name)

    def _update_csv(self, csv_rows: List[List[Any]]) -> None:
        """Update the CSV file with new configurations."""
        file_exists = os.path.exists(self.csv_file_name)

        # Step 1: Determine the existing field names in the CSV file
        if file_exists:
            with open(self.csv_file_name, mode='r', newline='') as file:
                reader = csv.DictReader(file)
                existing_field_names = reader.fieldnames
        else:
            existing_field_names = []

        # Ensure rows get written to correct column to these fields
        field_names = ["ID"] + list(self.param_settings.keys()) + ["Filename"]

        csv_data = zip(field_names, zip(*csv_rows))

        with open(self.csv_file_name, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)

            if not file_exists:
                writer.writerow(field_names)

            writer.writerows([k, *v] for k, v in csv_data)

    def generate_configs(self, N_tries_for_unique=1000):
        return self._generate_configs(N_tries_for_unique)

    def _generate_configs(self, N_tries_for_unique) -> None:
        csv_data = []

        # Generate new files num files number of times
        for _ in range(self.num_files):

            counter = 0
            # limit of searching 1000 times for a new configuration
            while counter < N_tries_for_unique:
                new_config, row = self.get_new_config()

                # stop looking for more configurations if the new config was found to be unique
                if new_config not in self._configs:
                    break

                else:
                    print("Trying again to find a new configuration")
                    counter += 1

            if counter > 1000:
                warnings.warn("Too many attempts to generate a new configuration. Stopping now", RuntimeWarning)
                break

            self._configs.append(new_config)
            filename = f"config_{self.id_counter}.json"
            row.append(filename)
            csv_data.append(row)

            with open(os.path.join(self.target_directory, filename), 'w') as json_file:
                json.dump(new_config, json_file, indent=4)

            self.id_counter += 1

        self._update_csv(csv_data)

    def get_new_config(self):
        # initialize new file based on template base config
        new_config = self._base_config.copy()
        # initialize the row for writing to csv with the id counter number
        row = [self.id_counter]
        # iterate through all the parameter keys and its values in param settings
        for param, settings in self.param_settings.items():

            # if param is rc then modify the rcs in the nonlocal bonds list
            if param == "rc":
                new_config["nonlocal_bonds"] = self._randomize_rcs(new_config["nonlocal_bonds"], settings)

                # append value for nonlocal bonds to the row for this file
                row.append(str(new_config["nonlocal_bonds"]))

            else:

                # check if param value is given by range or by values
                if 'range' in settings:
                    # pick values from a uniform distribution within the range given
                    new_config[param] = random.uniform(*settings['range'])

                elif 'values' in settings:
                    # pick values from the set of possible choices randomly
                    new_config[param] = random.choice(settings['values'])

                # append value for this param to the row for this file
                row.append(new_config[param])

        return new_config, row

    @staticmethod
    def _randomize_rcs(nonlocal_bonds, settings) -> list[list[Any]]:
        if 'range' in settings:
            # pick values from a uniform distribution within the range given for each bond in nonlocal bonds

            return [[bond[0], bond[1], round(random.uniform(*settings['range']), 2)] for bond in nonlocal_bonds]

        elif 'values' in settings:
            # pick values from the set of possible choices randomly
            return [[bond[0], bond[1], random.choice(settings['values'])] for bond in nonlocal_bonds]

        else:
            raise ValueError('settings does not contain correct rc information')


class ConfigGeneratorDriver(ConfigGenerator):
    def __init__(self, settings_config_file: str = "settings.cfg"):

        config = configparser.ConfigParser()
        config.read(settings_config_file)

        # Load file master settings
        if 'master_settings' in config:
            # initialize super class with the master settings
            super().__init__(
                base_config=config.get('master_settings', 'template_structure', fallback="template.json"),
                target_directory=config.get('master_settings', 'target_directory', fallback="generated_configs"),
                csv_file_name=config.get('master_settings', 'csv_name', fallback="configurations.csv"),
                num_files=config.getint('master_settings', 'num_files', fallback=5)
            )

            self.N_tries_for_unique = config.getint('master_settings', 'N_tries_for_unique', fallback=2000)

        # if not given then just use defaults
        else:
            super().__init__()

        # read in the param settings
        if 'param_settings' in config:
            for param, settings_str in config['param_settings'].items():
                settings_type, *values = settings_str.split(', ')
                if settings_type == 'range':
                    self.param_settings[param] = {'range': (float(values[0]), float(values[1]))}
                elif settings_type == 'values':
                    values = [json.loads(value) for value in values]
                    self.param_settings[param] = {'values': values}

    def generate_configs(self, N_tries_for_unique=None):

        # in case N_tries is given manually, use the manual value else use settings value
        if N_tries_for_unique is None:
            N_tries_for_unique = self.N_tries_for_unique

        return self._generate_configs(N_tries_for_unique)

    def submit_slurm(self):
        slurm_command = ('sbatch'
                         '')
