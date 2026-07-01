#!/usr/bin/env python3
"""
3_run_purity.py

Run CalicoST tumor purity estimation for a tissue. Expects parsed_inputs/
checkpoint files to already exist (built by 2_build_calicost_inputs.py).

CalicoST models each spot as a mixture of tumor and normal cells. With
n_clones=5 and no prior tumorprop_file, it infers tumor proportion jointly
with clone labels. Output is used as input to 4_run_cna.py.

Usage:
    python 3_run_purity.py <config_yaml>
    Example: python 3_run_purity.py tissue/lowseq_489/config_purity.yaml
"""

import sys
import os
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()

CALICOST_SRC = "/projectnb/paxlab/presh/software/CalicoST/src"
sys.path.insert(0, CALICOST_SRC)


def main():
    if len(sys.argv) < 2:
        print("Usage: python 3_run_purity.py <config_yaml>")
        sys.exit(1)

    config_file = sys.argv[1]
    if not os.path.exists(config_file):
        logger.error(f"Config file not found: {config_file}")
        sys.exit(1)

    # Verify parsed_inputs exist before launching CalicoST
    from calicost.arg_parse import read_configuration_file
    config = read_configuration_file(config_file)
    out_dir = config["output_dir"]

    parsed_dir = Path(out_dir) / "parsed_inputs"
    required = [
        "table_bininfo.csv.gz", "table_rdrbaf.csv.gz", "table_meta.csv.gz",
        "adjacency_mat.npz", "smooth_mat.npz", "exp_counts.pkl"
    ]
    missing = [f for f in required if not (parsed_dir / f).exists()]
    if missing:
        logger.error(f"Missing parsed_inputs files: {missing}")
        logger.error(f"Run 2_build_calicost_inputs.py first, then symlink "
                     f"parsed_inputs/ into {out_dir}/parsed_inputs")
        sys.exit(1)

    logger.info(f"All parsed_inputs present in {parsed_dir}")
    logger.info(f"Running CalicoST purity estimation with config: {config_file}")
    logger.info(f"Output dir: {out_dir}")

    from calicost.calicost_main import main as calicost_main
    calicost_main(config_file)

    # Report expected output location (clone_labels.tsv, cnv_*.tsv)
    n_clones = config.get("n_clones", 5)
    sw = config.get("spatial_weight", 1.0)
    clone_labels = Path(out_dir) / f"clone{n_clones}_rectangle0_w{sw:.1f}" / "clone_labels.tsv"
    if clone_labels.exists():
        logger.info(f"Purity estimation complete. Clone labels: {clone_labels}")
    else:
        logger.warning(f"Expected clone_labels.tsv not found at: {clone_labels}")
        logger.warning("Check CalicoST logs for errors.")


if __name__ == "__main__":
    main()
