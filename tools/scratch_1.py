import h5py

def get_configs(name):
    """
    Function to load the config counts for the last
    Parameters
    ----------
    name: str: the name of the simulation

    Returns
    -------
    np.array: the distances for the transient bond active now
    """
    # load the distances for this simulation
    with h5py.File(f"{name}.h5", 'r') as f:
        total_configs = int(f['config_count'][-1].sum())
        config_set = f['config']['int'][-total_configs:]
        s_bias = f['s_bias'][0]

    return config_set, s_bias

def sbias_estimator(configs):
    return


if __name__ == '__main__':
    name = 'hybridmc_0_0000000000_0000000001'

    dataset = get_configs(name)

