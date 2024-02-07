from . import ConfigGenerator

from ..helpers.data_processing_helpers import read_json_file

generator = ConfigGenerator()

# Set the base configuration
generator.set_base_config(read_json_file('template.json'))

# Set the number of files to generate
generator.set_num_files(5)

# Add parameter settings
generator.add_parameter_setting("rc", {"range": (1.0, 2.0)})
generator.add_parameter_setting("nbeads", {"values": [15, 20, 25]})
generator.add_parameter_setting("del_t", {"range": (0.5, 1.5)})
generator.add_parameter_setting("total_iter", {"values": [1200, 1500, 1800]})

# Generate the configurations and the CSV file
generator.generate_configs()