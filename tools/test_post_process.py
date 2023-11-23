import argparse
from post_processing import diff_s_bias, avg_s_bias, mfpt
import os


def main(args):
    os.chdir(args.path)
    # Obtain the differences in the sbias for each transition
    #diff_s_bias.get_diff_sbias()
    # Obtain the average sbias for each bonding state
    #avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv", structure_sim_json=args.json)
    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt(rewrite=True)
    # put together the mfpts in one file
    mfpt.compile_mfpts()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", help="results_path", default='test_st_morestair')
    #parser.add_argument("--json", help="json_path", default='../test_st.json')
    args = parser.parse_args()
    main(args)