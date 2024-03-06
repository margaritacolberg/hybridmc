import os
from dataclasses import make_dataclass

import h5py
from scipy.stats import bootstrap
import csv
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from data_processing_helpers import set_defaults, if_stair

else:
    from .data_processing_helpers import set_defaults, if_stair


class ConfigEntropyDiffBoot:
    """
    Class for providing configuration entropy difference and bootstrapping results to find associated error in it.

    -Allows you to specify bootstrapping procedure parameters.
        used to obtain bootstrapping result for simulation run s_bias value using a single trajectory result.

    -Works for a single transition result with full functionality for a single transition result without
        any staircase potential.

    -Also contains utility functions to help with configuration entropy calculations.

    """

    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        """
        Initialization routine for ConfigBoot:

        Upon initialization, this class stores the bootstrapping result which you can get using bootstrap_result.

        The s_bias and config_set parameters are ignored if the simulation name has been provided. The manual
            override for these two only works if you do not specify the simulation name.

        (Note: s_bias is the biased entropy for each configuration.)

        Parameters
        ----------
        simulation_name: str: Basename to the HDF5 file containing the simulation data

        s_bias: float: The s_bias for the simulation run given manually. ignored if simulation_name provided

        config_set: array: An array of integers representing the configurations explored. Ignored if simulation_name given

        bootstrap_kwargs: dict: Keyword arguments for the bootstrapping procedure

        Returns:
        _________
        The bootstrapping result, configuration array, the s_bias value in a dataclass
        """

        # check if s_bias and config_set provided
        if not simulation_name:
            if isinstance(s_bias, (int, float)) and config_set:
                self.s_bias = s_bias
                self.config_set = config_set
            else:
                raise ValueError('Simulation data path not specified, provide both s_bias and config_set parameters')
        # simulation name is provided
        else:
            # set hdf5 name to retrieve configuration and sbias info if given
            self.simulation_name = self.get_base_simulation_id(simulation_name)
            # find the config set and corresponding s_bias value using simulation data file
            self._set_configs()

        # set the bootstrap kwargs
        self.bootstrap_kwargs = bootstrap_kwargs

        # add the minimum required kwargs parameters for bootstrap computation
        set_defaults(self.bootstrap_kwargs, dict(confidence_level=0.95,n_resamples=2000, batch=1,
                                                 random_state=np.random.default_rng()))

        # initialize bootstrap result
        self.bootstrap_result = None

        # compute bootstrap result
        self.compute_bootstrap(**self.bootstrap_kwargs)

    def __repr__(self):
        return (f"Simulation had s_bias: {self.s_bias_mean}\n"
                f"Standard Error: {self.bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.bootstrap_result.confidence_interval.low} to "
                f"{self.bootstrap_result.confidence_interval.high}")

    def __str__(self):
        return (f"Simulation had s_bias_mean: {self.s_bias_mean}\n"
                f"Standard Error: {self.bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.bootstrap_result.confidence_interval.low} to "
                f"{self.bootstrap_result.confidence_interval.high}")

    def __call__(self):
        fields = ['bootstrap_result', 'config_set', 's_bias', 's_bias_mean']
        BootstrapInfo = make_dataclass('BootstrapInfo', fields)
        return BootstrapInfo(bootstrap_result=self.bootstrap_result,
                             config_set=self.config_set[0], s_bias=self.s_bias_mean)

    @staticmethod
    def _get_config_from_h5(h5_file_name):
        """
        Function to load the config counts for the iteration of the simulation and the sbias value for it
        """
        # load the distances for this simulation
        with h5py.File(h5_file_name, 'r') as f:
            total_configs = int(f['config_count'][-1].sum())
            return f['config']['int'][-total_configs:], f['s_bias'][0]


    @staticmethod
    def sbias_from_config_estimator(config_data, s_bias):
        """

        Parameters
        ----------
        config_data: List[int]: The list of configurations (0,1) visited during simulation.
        s_bias: float: The s_bias value simulation found to have based on the original configuration.

        Returns
        -------

        """
        # add sbias and log(number of 0 states in config_data / number of 1 states)
        return float(s_bias) + np.log((config_data.size - config_data.sum()) / config_data.sum())

    def _set_configs(self):
        """
        Sets the configuration bitstring list in a form compatible with the scipy bootstrap method,
        Also sets the s_bias value.
        """
        config_set, self.s_bias = self._get_config_from_h5(f"{self.simulation_name}.h5")
        self.s_bias_mean = self.sbias_from_config_estimator(config_set, self.s_bias)
        self.config_set = (config_set,)

    def _s_bias_estimator(self, config_data):
        return self.sbias_from_config_estimator(config_data, self.s_bias)

    def compute_bootstrap(self, **kwargs):
        set_defaults(kwargs, dict(confidence_level=0.95,
                                  random_state=np.random.default_rng()))

        self.bootstrap_result = bootstrap(data=self.config_set,
                                          statistic=self._s_bias_estimator,
                                          **kwargs)

    def bootstrap_hist(self, **kwargs):
        fig, ax = plt.subplots()
        ax.hist(self.bootstrap_result.bootstrap_distribution, **kwargs)
        ax.set_title('Bootstrap Distribution')
        ax.set_xlabel('s_bias value')
        ax.set_ylabel('frequency')
        plt.show()

    @staticmethod
    def get_base_simulation_id(simulation_name):
        simulation_name = simulation_name.split('.')[0]
        return '_'.join(os.path.basename(simulation_name).split('_')[:4])

    def get_diff_sbias_output(self):
        # Initialize diff_s_bias output
        diff_data = []
        # Obtain the input and output configs for the simulation as first two columns
        diff_data += self.simulation_name.split('_')[2:4]

        # Add the bootstrap result to the diff data output
        #diff_data += [self.s_bias, self.bootstrap_result.standard_error, self.bootstrap_result.confidence_interval.low,
        #              self.bootstrap_result.confidence_interval.high]

        diff_data += [self.s_bias_mean, self.bootstrap_result.standard_error,
                      self.bootstrap_result.confidence_interval.low,
                      self.bootstrap_result.confidence_interval.high]

        return diff_data

    def write_bootstrap(self, csv_name):
        """
        Only works for single transition writing.
        Parameters
        ----------
        csv_name: Name of the csv file to write output

        Returns
        -------
        None
        """

        if os.path.isfile(csv_name):
            with open(csv_name, 'r') as fread:
                reader = csv.reader(fread)
                diff_data = list(reader)
        else:
            diff_data = [[]]

        if not len(diff_data):
            diff_data = [[]]

        with (open(csv_name, 'w') as fwrite):
            if len(diff_data[0]) == 3:
                diff_data[0].append(self.bootstrap_result.confidence_interval)
            elif len(diff_data[0]) == 4:
                diff_data[0][3] = self.bootstrap_result.confidence_interval
            else:
                # Obtain the input and output configs for the simulation as first two columns and boostrap result as
                # the third if csv file found not to have this info already.
                diff_data[0] = self.get_diff_sbias_output()

            writer = csv.writer(fwrite)
            writer.writerow(diff_data[0])


class StairConfigEntropyDiffBoot(ConfigEntropyDiffBoot):
    """
    Subclass of ConfigEntropyDiffBoot to provide bootstrapping procedure for a stair configuration
    """

    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        super().__init__(simulation_name, s_bias, config_set, **bootstrap_kwargs)

    @staticmethod
    def stair_s_bias(stair_sbias_list):
        exp_s = []
        for i in range(1, len(stair_sbias_list) + 1):
            [exp_s.append(np.exp(sum(j))) for j in combinations(stair_sbias_list, i)]
        return np.log(np.sum(exp_s))

    def _set_configs(self):

        # get the staircase rc values and add 0 to indicate the final rc
        stair_rc_list = if_stair(self.simulation_name, files=os.listdir()) + [0]

        # Get list of the staircase s_bias values and associated configuration bitstring
        # if it is one of the intermediate staircase rc values use rc value to find simulation path
        # else use the basename

        sim_paths = [f"{self.simulation_name}_{rc}" if rc else self.simulation_name for rc in stair_rc_list]
        self.config_set, self.s_bias_list = list(zip(*(self._get_config_from_h5(f"{file_id}.h5") for file_id in sim_paths)))

        # Set the sbias as the actual sbias value for the simulation: This is the "mean" sbias here.
        self.s_bias_mean = self._stair_sbias_estimator(self.config_set)

    def _stair_sbias_estimator(self, config_data):
        stair_sbias_list = []
        for idx, config in enumerate(config_data):
            stair_sbias_list.append(self.sbias_from_config_estimator(config, self.s_bias_list[idx]))

        return self.stair_s_bias(stair_sbias_list)

    def _s_bias_estimator(self, *config_data):
        return self._stair_sbias_estimator(config_data)


def main(params):
    entropy_diff_bootstrap = StairConfigEntropyDiffBoot(simulation_name=params.simulation_name)
    print(entropy_diff_bootstrap)
    entropy_diff_bootstrap.bootstrap_hist()
    entropy_diff_bootstrap.write_bootstrap(params.output_file)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--simulation_name', type=str, default='hybridmc_0_00_10.h5')
    parser.add_argument('--s_bias', type=float, default=0)
    parser.add_argument('--config_set', type=float, default=0)
    parser.add_argument('--output-file', type=str, default='diff_s_bias_with_error.csv')
    args = parser.parse_args()

    main(args)
