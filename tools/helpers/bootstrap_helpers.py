import csv
import os
import h5py
import numpy as np
from scipy.stats import bootstrap
import csv
import matplotlib.pyplot as plt

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from data_processing_helpers import set_defaults

else:
    from .data_processing_helpers import set_defaults


class ConfigBoot:
    """
    Class used to obtain bootstrapping result for simulation run s_bias value (using a single trajectory result).

    Upon initialization, this class stores the bootstrapping result which you can get using bootstrap_result.

    Attributes:
        s_bias: float: The s_bias for the original simulation run for config_set
        config_set: List: The set of configurations from the original simulation run
        bootstrap_result: The bootstrapping result

    Methods:
        compute_bootstrapping(**kwargs): This value for bootstrap result can be updated with other parameters.

        boostrap_hist(**kwargs): View a histogram of the results of the bootstrapping.

        write_bootstrapping(csv_name): Write out the results of the bootstrapping to a csv_file.

    Returns:
            The bootstrapping result, configuration array, the s_bias value as a tuple.

    """

    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        """
        Initialize the ConfigBoot. This allows you to specify the parameters of the bootstrapping procedure.

        The s_bias and config_set parameters are ignored if the simulation name has been provided. The manual
        override for these two only works if you do not specify the simulation name.

        Parameters
        ----------
        simulation_name: str: Basename to the HDF5 file containing the simulation data
        s_bias: float: The s_bias for the simulation run given manually. ignored if simulation_name provided
        config_set: array: An array of integers representing the configurations explored. Ignored if simulation_name given
        bootstrap_kwargs: dict: Keyword arguments for the bootstrapping procedure
        """

        # check if s_bias and config_set provided
        if not simulation_name:
            if s_bias and config_set:
                self.s_bias = s_bias
                self.config_set = config_set
            else:
                raise ValueError('Simulation data path not specified, provide both s_bias and config_set parameters')
        # simulation name is provided
        else:
            # set hdf5 name to retrieve configuration and sbias info if given
            self.simulation_name = simulation_name
            # find the config set and corresponding s_bias value using simulation data file
            self._set_configs()

        # set the bootstrap kwargs
        self._bootstrap_kwargs = bootstrap_kwargs

        # add the minimum required kwargs parameters for bootstrap computation
        set_defaults(self._bootstrap_kwargs, dict(confidence_level=0.9,
                                                  random_state=np.random.default_rng()))
        # initialize bootstrap result
        self.__bootstrap_result = None
        # get the bootstrap value using the configs
        self.compute_bootstrap(**self._bootstrap_kwargs)

    def __repr__(self):
        return (f"Simulation had s_bias: {self.s_bias}\n"
                f"Standard Error: {self.__bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.__bootstrap_result.confidence_interval.low} to "
                f"{self.__bootstrap_result.confidence_interval.high}")

    def __str__(self):
        return (f"Simulation had s_bias: {self.s_bias}\n"
                f"Standard Error: {self.__bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.__bootstrap_result.confidence_interval.low} to "
                f"{self.__bootstrap_result.confidence_interval.high}")

    def __call__(self):
        return self.__bootstrap_result, self.config_set, self.s_bias

    def _set_configs(self):
        """
        Function to load the config counts for the iteration of the simulation and the sbias value for it
        """
        # load the distances for this simulation
        with h5py.File(f"{self.simulation_name}.h5", 'r') as f:
            total_configs = int(f['config_count'][-1].sum())
            self.config_set = f['config']['int'][-total_configs:]
            self.s_bias = f['s_bias'][0]

    def _s_bias_estimator(self, config_data):
        return self.s_bias + np.log(config_data.sum() / (config_data.size - config_data.sum()))

    @property
    def bootstrap_result(self):
        return self.__bootstrap_result

    def compute_bootstrap(self, **kwargs):
        # Add minimum required arguments with default values to the kwargs passed
        # (needed when run again with different kwargs)
        set_defaults(kwargs, dict(confidence_level=0.9,
                                  random_state=np.random.default_rng()))

        self.__bootstrap_result = bootstrap(data=(self.config_set,),
                                            statistic=self._s_bias_estimator,
                                            **kwargs)

    def bootstrap_hist(self, **kwargs):
        fig, ax = plt.subplots()
        ax.hist(self.__bootstrap_result.bootstrap_distribution, **kwargs)
        ax.set_title('Bootstrap Distribution')
        ax.set_xlabel('s_bias value')
        ax.set_ylabel('frequency')
        plt.show()

    def write_bootstrap(self, csv_name):
        with open(csv_name, 'r') as fread:
            reader = csv.reader(fread)
            diff_data = list(reader)

        with open(csv_name, 'w') as fwrite:
            if len(diff_data[0]) == 3:
                diff_data[0].append(self.bootstrap_result.confidence_interval)
            elif len(diff_data[0]) == 4:
                diff_data[0][3] = self.bootstrap_result.confidence_interval

            writer = csv.writer(fwrite)
            writer.writerow(diff_data[0])


def main(args):
    config_boot = ConfigBoot(args.simulation_name, s_bias=args.s_bias, config_set=args.config_set)
    print(config_boot)
    config_boot.bootstrap_hist()
    config_boot.write_bootstrap(f'{os.getcwd()}/diff_s_bias.csv')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--simulation_name', type=str, default='hybridmc_0_0000000000_0000000001')
    parser.add_argument('--s_bias', type=float, default=0)
    parser.add_argument('--config_set', type=float, default=0)
    parser.add_argument('--output-file', type=str, default='diff_s_bias.csv')
    args = parser.parse_args()

    main(args)
