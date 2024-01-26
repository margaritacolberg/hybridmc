import h5py
import numpy as np
from scipy.stats import bootstrap


class SimConfigs:
    def __init__(self, simulation_name):
        self.simulation_name = simulation_name
        config_set, self.s_bias = self.get_configs()
        self.config_set = (config_set,)
        self.bootstrap_result = self.get_bootstrap()

    def get_configs(self):
        """
        Function to load the config counts for the iteration of the simulation
        """
        # load the distances for this simulation
        with h5py.File(f"{self.simulation_name}.h5", 'r') as f:
            total_configs = int(f['config_count'][-1].sum())
            config_set = f['config']['int'][-total_configs:]
            s_bias = f['s_bias'][0]

        return config_set, s_bias

    def sbias_estimator(self, config_data):
        return self.s_bias + np.log(config_data.sum() / (config_data.size - config_data.sum()))

    def get_bootstrap(self):
        return bootstrap(data=self.config_set,
                         statistic=self.sbias_estimator,
                         confidence_level=0.9,
                         random_state=np.random.default_rng())


if __name__ == '__main__':
    name = 'hybridmc_0_0000000000_0000000001'

    config = SimConfigs(name)
    print(config.bootstrap_result)
    import pandas as pd
    import os
    os.chdir('../')
    df = pd.read_csv('summary_data.csv')
    print(df.query('ID == "sbias"')['mean'])
