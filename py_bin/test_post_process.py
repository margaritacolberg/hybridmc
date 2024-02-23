import argparse
from py_tools.post_processing import diff_s_bias, avg_s_bias, mfpt, constructPaths


def main():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias("diff_s_bias_new.csv")
    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias_new.csv")
    # Obtain the mfpt for each bonding state
    #mfpt.get_mfpt(rewrite=True)
    # put together the mfpts in one file
    #mfpt.compile_mfpts()
    constructPaths.diff_sbias_state_function_check("diff_s_bias_new.csv", csv_out="sig_diff_new.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", help="json_path", default='../test_st.json')
    args = parser.parse_args()
    main()