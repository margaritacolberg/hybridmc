import os
from dataclasses import make_dataclass

import h5py
from scipy.stats import bootstrap
from scipy.stats._common import ConfidenceInterval
import csv
import matplotlib.pyplot as plt
from sklearn.utils import resample

import numpy as np

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from data_processing_helpers import set_defaults
    from py_tools.post_processing.diff_s_bias import stair_s_bias
    from py_tools.helpers.mfpt_helpers import if_stair

else:
    from .data_processing_helpers import set_defaults
    from ..post_processing.diff_s_bias import stair_s_bias
    from .mfpt_helpers import if_stair


class Boot:
    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        """
        Parent class for bootstrapping. This allows you to specify the parameters of the bootstrapping procedure.

        Child classes must specify the statistic for bootstrapping. While the configuration can be redefined,
        this class will automatically be using the h5 file provided to set the sample set.

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
            if isinstance(s_bias, (int, float)) and config_set:
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
        self.bootstrap_kwargs = bootstrap_kwargs

        # add the minimum required kwargs parameters for bootstrap computation
        set_defaults(self.bootstrap_kwargs, dict(confidence_level=0.9,
                                                 random_state=np.random.default_rng()))

        # initialize result for bootstrap_result
        self.bootstrap_result = None

        # compute bootstrap
        self.compute_bootstrap(**self.bootstrap_kwargs)

    def __repr__(self):
        return (f"Simulation had s_bias: {self.s_bias}\n"
                f"Standard Error: {self.bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.bootstrap_result.confidence_interval.low} to "
                f"{self.bootstrap_result.confidence_interval.high}")

    def __str__(self):
        return (f"Simulation had s_bias: {self.s_bias}\n"
                f"Standard Error: {self.bootstrap_result.standard_error}\n"
                f"Confidence interval: {self.bootstrap_result.confidence_interval.low} to "
                f"{self.bootstrap_result.confidence_interval.high}")

    def __call__(self):
        make_dataclass('')
        return self.bootstrap_result, self.config_set, self.s_bias

    def _get_config_from_h5(self, h5_file_name=None):
        """
        Function to load the config counts for the iteration of the simulation and the sbias value for it
        """

        # if a manual sim_name provided, use it
        if h5_file_name:
            sim_name = h5_file_name
        else:
            # initialize sim file name with the simulation name
            sim_name = self.simulation_name

        # load the distances for this simulation
        with h5py.File(f"{sim_name}.h5", 'r') as f:
            total_configs = int(f['config_count'][-1].sum())
            return f['config']['int'][-total_configs:], f['s_bias'][0]

    def _set_configs(self):
        self.config_set, self.s_bias = self._get_config_from_h5()
        self.config_set = (self.config_set,)

    def _s_bias_estimator(self, config_data):
        pass

    def compute_bootstrap(self, **kwargs):
        set_defaults(kwargs, dict(confidence_level=0.9,
                                  random_state=np.random.default_rng()))

        self.bootstrap_result = bootstrap(data=self.config_set,
                                          statistic=self._s_bias_estimator,
                                          **kwargs)

    def __bootstrap_hist(self, **kwargs):
        fig, ax = plt.subplots()
        ax.hist(self.bootstrap_result.bootstrap_distribution, **kwargs)
        ax.set_title('Bootstrap Distribution')
        ax.set_xlabel('s_bias value')
        ax.set_ylabel('frequency')
        plt.show()

    def bootstrap_hist(self, **kwargs):

        return self.__bootstrap_hist(**kwargs)


class ConfigBoot(Boot):
    """
    SubClass of Boot, used to obtain bootstrapping result for simulation run s_bias value (using a single trajectory result).

    Upon initialization, this class stores the bootstrapping result which you can get using bootstrap_result.

    This class works for a single transition result with full functionality for a single transition result.

    Attributes:
        simulation_name: str: Basename to the HDF5 file containing the simulation data. Used by default even if sbias and config set given manually
        s_bias: optional float: The s_bias for the original simulation run for config_set in case simulation name not given
        config_set: optional List: The set of configurations from the original simulation run in case simulation name not given
        bootstrap_result: The bootstrapping result

    Methods:
        compute_bootstrapping(**kwargs): This value for bootstrap result can be updated with other parameters.

        boostrap_hist(**kwargs): View a histogram of the results of the bootstrapping.

        write_bootstrapping(csv_name): Write out the results of the bootstrapping to a csv_file.

    Returns:
            The bootstrapping result, configuration array, the s_bias value as a tuple.

    """

    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        super().__init__(simulation_name, s_bias, config_set, **bootstrap_kwargs)

    @staticmethod
    def sbias_from_config_estimator(config_data, s_bias):
        return s_bias + np.log(config_data.sum() / (config_data.size - config_data.sum()))

    def _s_bias_estimator(self, config_data):
        return self.sbias_from_config_estimator(config_data, self.s_bias)

    def __write_bootstrap(self, csv_name):
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

    def write_bootstrap(self, csv_name):
        return self.__write_bootstrap(csv_name)


class StairConfigBoot(ConfigBoot):
    """
    Subclass of ConfigBoot to provide bootstrapping procedure for a stair configuration
    """

    def __init__(self, simulation_name=None, s_bias=None, config_set=None, **bootstrap_kwargs):
        super().__init__(simulation_name, s_bias, config_set, **bootstrap_kwargs)

    def _set_configs(self):

        sim_basename = os.path.basename(self.simulation_name)
        if len(sim_basename.split('_')) == 4:
            ref_sim_id = sim_basename
        else:
            ref_sim_id = '_'.join(sim_basename[:4])

        # get the staircase rc values
        stair_rc_list = if_stair(ref_sim_id, files=os.listdir())
        # add 0 to indicate the final rc
        stair_rc_list += [0]

        self.s_bias_list = []
        self.config_set = []
        for rc in stair_rc_list:
            # if it is one of the intermediate staircase rc values
            if rc:
                file_id = f"{ref_sim_id}_{rc}"
            # if it is the final step
            else:
                file_id = ref_sim_id

            config, s_bias = self._get_config_from_h5(file_id)

            self.config_set.append(config)
            self.s_bias_list.append(s_bias)

        self.s_bias = str(stair_s_bias(self.s_bias_list))

    @staticmethod
    def resample_MultiDim(dataset, **kwargs):
        resample_result = []
        for sample in dataset:
            # change to array in case sample in data is not an array
            if not isinstance(sample, np.ndarray):
                sample = np.array(sample)

            resample_result.append(resample(sample, **kwargs))

        return resample_result

    def __resample_configs(self, n_resamples, **kwargs):
        return [self.resample_MultiDim(self.config_set, **kwargs) for _ in range(n_resamples)]

    def _s_bias_estimator(self, config_data):
        stair_sbias_list = []
        for idx, config in enumerate(config_data):
            stair_sbias_list.append(super().sbias_from_config_estimator(config, self.s_bias_list[idx]))

        return stair_s_bias(stair_sbias_list)

    def stair_bootstrap(self, **kwargs):

        # initialize bootstrap result output class
        fields = ['mean', 'confidence_interval', 'standard_error']
        BootstrapResult = make_dataclass("BootstrapResult", fields)

        # prepare resample set
        data = self.__resample_configs(n_resamples=kwargs['n_boots'])

        # go through each resample set row and find estimated sbias
        bootstrap_sbias = np.zeros(kwargs['n_boots'])  # initialize sbias result array

        for idx, row in enumerate(data):
            # calculate the sbias for the row configs
            bootstrap_sbias[idx] = self._s_bias_estimator(row)

        # Calculate confidence interval
        ci_l = np.quantile(bootstrap_sbias, (1 - kwargs['confidence_level']) / 2)
        ci_u = np.quantile(bootstrap_sbias, (1 + kwargs['confidence_level']) / 2)

        # Return summary result
        return BootstrapResult(
            confidence_interval=ConfidenceInterval(ci_l, ci_u),
            mean=np.mean(bootstrap_sbias),
            standard_error=np.std(bootstrap_sbias, ddof=1, axis=-1)
        )

    def compute_bootstrap(self, **kwargs):

        set_defaults(kwargs, dict(confidence_level=0.9,
                                  random_state=np.random.default_rng(),
                                  n_boots=10
                                  ))

        self.bootstrap_result = self.stair_bootstrap(data=self.config_set,
                                                     statistic=self._s_bias_estimator,
                                                     **kwargs)


def main(params):
    bootstrap_result = StairConfigBoot(simulation_name=params.simulation_name, s_bias=params.s_bias,
                                       config_set=params.config_set, n_boots=20, confidence_level=0.95)
    print(bootstrap_result)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--simulation_name', type=str, default='hybridmc_3_1000000110_1010000110')
    parser.add_argument('--s_bias', type=float, default=0)
    parser.add_argument('--config_set', type=float, default=0)
    parser.add_argument('--output-file', type=str, default='diff_s_bias.csv')
    args = parser.parse_args()

    main(args)
