from .data_processing_helpers import *
from .run_layer_helpers import run_sim, run_stairs


def run_layer(common_data, in_queue, out_queue, seed_increment, exe):
    """
    Run a single layer of simulations for all transitions in that layer in parallel. The assumption
    for a "layer" is that all transitions in that layer have the same permanent bonds turned on.

    Parameters
    ----------
    common_data: dict: JSON input for the HMC program as a python dictionary
    in_queue: multiprocessing.Queue: queue of transitions to run
    out_queue: multiprocessing.Queue: queue of transitions that have been run
    seed_increment: int: seed increment for each trajectory
    exe: str: path to the HMC executable

    Returns
    -------
    None
    """
    while True:
        item = in_queue.get()
        if item is None:
            break

        i, bits_in, count, input_hdf5 = item

        # copy of initial state for every nonbonded pair of beads
        bits_out = list(bits_in)

        # flip bit to form bond
        bits_out[i] = True
        bits_out = tuple(bits_out)
        bonds_in = bits_to_bonds(bits_in, common_data['nonlocal_bonds'])
        layer = len(bonds_in)

        # Define the output files identifier for this transition (used to name hdf5 files, json files and log files)
        output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}'

        # Set the transient bond to form this trajectory
        common_data['transient_bonds'] = [common_data['nonlocal_bonds'][i]]

        # Set the existing bonds as permanent bonds
        common_data['permanent_bonds'] = bonds_in

        # Obtain configuration information
        common_data['config_in'] = int(format_bits(bits_in), 2)
        common_data['config_out'] = int(format_bits(bits_out), 2)

        # Set the seed
        common_data['seeds'] = [count, seed_increment]

        # Check for stairs using wang_landau (WL) run
        stair_bp, stair_rc_list = stair_check(common_data, output_name, input_hdf5)

        # Check if staircase needed and do a staircase run for structure; if not do regular run
        if stair_bp and stair_rc_list:

            # ensure cell size can accommodate the rc
            if (common_data["length"] / common_data["ncell"]) < stair_rc_list[0]:
                common_data["length"] = stair_rc_list[0] * common_data["ncell"] + 0.1

            run_stairs(common_data, input_hdf5, output_name, exe, stair_rc_list)

        else:
            run_sim(common_data, input_hdf5, output_name, exe)

        # Put the output file name in the queue to indicate that this transition has been run
        out_queue.put((bits_out, f'{output_name}.h5'))
