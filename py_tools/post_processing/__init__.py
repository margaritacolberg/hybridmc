from . import diff_s_bias, avg_s_bias, mfpt


def post_processing():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias(out_csv='diff_s_bias.csv')
    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv")
    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt()
    # put together the mfpts in one file
    mfpt.compile_mfpts()
