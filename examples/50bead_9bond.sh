#!/bin/bash
#SBATCH --account=def-jmschofi
#SBATCH --time=9-0
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=256M

# abort script on error
set -e

[[ -n "$date" ]] || date=$(date -Is)

mkdir -p "run_$date".tmp
pushd "run_$date".tmp
python ../../tools/init_json.py ../../examples/50bead_9bond.json ../../release/hybridmc 1 ${SLURM_CPUS_PER_TASK:+--nproc="$SLURM_CPUS_PER_TASK"}
python ../../tools/diff_s_bias.py ../../examples/50bead_9bond.json
python ../../tools/avg_s_bias.py ../../examples/50bead_9bond.json diff_s_bias.csv
popd
mv "run_$date".tmp "run_$date"
