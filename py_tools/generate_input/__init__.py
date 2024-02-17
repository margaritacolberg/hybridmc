import configparser
import json
import random
import os
import warnings
from functools import cached_property
from typing import Any, Dict, List
import sqlite3
import subprocess
import tempfile


class ConfigGenerator:

    def __init__(self, base_config: Dict[str, Any] = None, num_files: int = 5, start_id=None,
                 target_directory: str = "Generated_configs", database_file_name: str = "configs.db"):

        self.base_config = base_config
        self.num_files = num_files
        self._configs = []  # To store generated configurations for uniqueness check
        self.param_settings = {}

        # set target directory
        self.target_directory = target_directory

        # Set the database name and find the id counter
        self.database_name = database_file_name
        self.id_counter = start_id

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

    @property
    def id_counter(self) -> int:
        return self._id_counter

    @cached_property
    def database_path(self):
        return os.path.join(self.target_directory, self.database_name)

    @id_counter.setter
    def id_counter(self, value):
        """Initialize id_counter with the last ID in the database plus one or given manual value."""
        if value:
            self._id_counter = value
        elif not os.path.exists(self.database_path):
            self._id_counter = 1
        else:
            db_manager = DatabaseManager(self.database_path)
            last_id = db_manager.get_max_id()
            self._id_counter = last_id + 1
            db_manager.close()

    def _update_database(self, rows: List[List[Any]]) -> None:
        """Update the CSV file with new configurations."""
        db_manager = DatabaseManager(self.database_path)

        # field_names =
        # Update the database with the new configurations
        db_manager.update_configurations(
            rows=rows,
            field_names=[el if el != "rc" else "nonlocal_bonds" for el in self.param_settings.keys()] + ["Filename"],
            default_values=self.base_config
        )

        # Close the database connection when done
        db_manager.close()

    def generate_configs(self, N_tries_for_unique=1000):
        return self._generate_configs(N_tries_for_unique)

    def _generate_configs(self, N_tries_for_unique) -> None:
        data = []

        # Generate new files num files number of times
        for _ in range(self.num_files):
            # initialize the row for writing to database with the id counter number
            row = []
            counter = 0
            # limit of searching 1000 times for a new configuration
            while counter < N_tries_for_unique:
                new_config = self.get_new_config(row)

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
            data.append(row)

            with open(os.path.join(self.target_directory, filename), 'w') as json_file:
                json.dump(new_config, json_file, indent=4)

            self.id_counter += 1

        self._update_database(data)

    def get_new_config(self, row):
        # initialize new file based on template base config
        new_config = self._base_config.copy()

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

        return new_config

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
                database_file_name=config.get('master_settings', 'database_name', fallback="configurations.db"),
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


class DatabaseManager:
    def __init__(self, db_file_name: str):
        self.db_file_name = db_file_name
        self.connection = sqlite3.connect(self.db_file_name)
        self.cursor = self.connection.cursor()

    def _ensure_table(self, field_names: List[str]):
        """Ensure the table exists and has the required columns."""
        # The ID column is defined as INTEGER PRIMARY KEY, which auto-increments.
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS configurations (ID INTEGER PRIMARY KEY)''')
        for field in field_names:
            try:
                self.cursor.execute(f'''ALTER TABLE configurations ADD COLUMN {field} TEXT''')
            except sqlite3.OperationalError:
                pass  # Column already exists

    def update_configurations(self, rows: List[List[Any]], field_names: List[str], default_values: Dict[str, Any]):
        """Update the database with new configurations, applying default values as necessary."""
        self._ensure_table(field_names)

        for row in rows:
            # Apply default values for any missing fields
            data = {field: default_values.get(field) for field in field_names}
            # Update data with actual row values
            data.update(dict(zip(field_names, row)))

            # Exclude 'ID' from the insert operation
            if 'ID' in data:
                del data['ID']

            # Prepare and execute the insert query
            columns = ', '.join(data.keys())
            placeholders = ', '.join(['?'] * len(data))
            values = tuple(data.values())
            query = f"INSERT INTO configurations ({columns}) VALUES ({placeholders})"
            self.cursor.execute(query, values)

        self.connection.commit()

    def get_max_id(self) -> int:
        """Fetch the maximum ID from the configurations table."""
        self.cursor.execute('SELECT MAX(ID) FROM configurations')
        max_id = self.cursor.fetchone()[0]
        return max_id if max_id is not None else 0

    def close(self):
        """Close the database connection."""
        self.connection.close()


class JobSubmitter:
    def __init__(self, account='def-jmschofi', job_name='get_training', cpus_per_task=4,
                 mem_per_cpu=500, time='0-01:00:00', Nconfigs="1-5", json_dir='simulation_configs',
                 out_dir="train_configs", exe="/scratch/vignesh9/hybridmc/py_bin/run.py",
                 hmc_exe="/scratch/vignesh9/hybridmc/release/hybridmc",
                 ):

        self.temp_script_path = None
        self.account = account
        self.job_name = job_name
        self.cpus_per_task = cpus_per_task
        self.mem_per_cpu = mem_per_cpu
        self.time = time
        self.Nconfigs = Nconfigs
        self.json_dir = json_dir
        self.out_dir = out_dir
        self.exe = exe
        self.hmc_exe = hmc_exe

    @property
    def json_dir(self) -> str:
        return self._json_dir

    @json_dir.setter
    def json_dir(self, value: str) -> None:
        self._json_dir = value
        os.makedirs(self._json_dir, exist_ok=True)

    @property
    def out_dir(self) -> str:
        return self._out_dir

    @out_dir.setter
    def out_dir(self, value: str) -> None:
        self._out_dir = value
        os.makedirs(self._out_dir, exist_ok=True)
        os.chdir(self._out_dir)
        self.json_dir = f"../{self.json_dir}"

    def create_job_script(self):
        # Create the SLURM script content
        slurm_script_content = f"""#!/bin/bash
#SBATCH --account={self.account}
#SBATCH --job-name={self.job_name}
#SBATCH --cpus-per-task={self.cpus_per_task}
#SBATCH --mem-per-cpu={self.mem_per_cpu}
#SBATCH --time={self.time}
#SBATCH --output=slurm_out/config_%A_%a.out
#SBATCH --error=slurm_err/config_%A_%a.err
#SBATCH --array={self.Nconfigs}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vignesh.rajesh@mail.utoronto.ca

# intitialize shell
source /home/vignesh9/.bashrc
module --force purge

# Capture start time
start_time=$(date +%s)

micromamba activate HMC

# python execute with time tracking
time python {self.exe} --json {os.path.join(self.json_dir, "config")}_"$SLURM_ARRAY_TASK_ID".json --exe {self.hmc_exe}

# Capture end time and calculate duration
end_time=$(date +%s)

duration=$((end_time - start_time))
echo "Job Duration: $duration seconds"

# Additional SLURM job information
echo "Detailed job information:"
scontrol show job $SLURM_JOB_ID --details

echo "Accounting information for the job:"
sacct -j $SLURM_JOB_ID

"""

        # Write the SLURM script to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_script:
            temp_script.write(slurm_script_content)
            self.temp_script_path = temp_script.name

    def submit_job(self):

        # Check if job script was made, if not create one
        if self.temp_script_path is None:
            self.create_job_script()
        # Submit the SLURM job using the script
        try:
            subprocess.run(['sbatch', self.temp_script_path], check=True)
        finally:
            # Ensure the temporary file is removed after submission
            os.remove(self.temp_script_path)
