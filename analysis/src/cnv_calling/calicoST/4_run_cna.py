#!/usr/bin/env python3
"""
4_run_cna.py

Run CalicoST CNA + clone calling for a tissue. Expects:
  - parsed_inputs/ checkpoint files (built by 2_build_calicost_inputs.py)
  - tumorprop_file from 3_run_purity.py output

Infers allele-specific integer copy numbers, clone labels per spot, and
phylogeography of cancer clones.

Usage:
    python 4_run_cna.py <config_yaml>
    Example: python 4_run_cna.py tissue/lowseq_489/config_cna.yaml
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
        print("Usage: python 4_run_cna.py <config_yaml>")
        sys.exit(1)

    config_file = sys.argv[1]
    if not os.path.exists(config_file):
        logger.error(f"Config file not found: {config_file}")
        sys.exit(1)

    from calicost.arg_parse import read_configuration_file
    config = read_configuration_file(config_file)
    out_dir = config["output_dir"]

    # Verify parsed_inputs
    parsed_dir = Path(out_dir) / "parsed_inputs"
    required = [
        "table_bininfo.csv.gz", "table_rdrbaf.csv.gz", "table_meta.csv.gz",
        "adjacency_mat.npz", "smooth_mat.npz", "exp_counts.pkl"
    ]
    missing = [f for f in required if not (parsed_dir / f).exists()]
    if missing:
        logger.error(f"Missing parsed_inputs files: {missing}")
        logger.error("Ensure parsed_inputs/ is symlinked or present in the CNA output_dir.")
        sys.exit(1)

    # Verify tumorprop_file
    tumorprop_file = config.get("tumorprop_file")
    if tumorprop_file and not os.path.exists(tumorprop_file):
        logger.error(f"tumorprop_file not found: {tumorprop_file}")
        logger.error("Run 3_run_purity.py first.")
        sys.exit(1)

    logger.info(f"Running CalicoST CNA calling with config: {config_file}")
    logger.info(f"Output dir: {out_dir}")
    if tumorprop_file:
        logger.info(f"Using tumor proportion file: {tumorprop_file}")
    else:
        logger.warning("No tumorprop_file in config — running without purity prior.")

    from calicost.calicost_main import main as calicost_main
    calicost_main(config_file)

    # Report expected outputs
    n_clones = config.get("n_clones", 3)
    sw = config.get("spatial_weight", 1.0)
    result_dir = Path(out_dir) / f"clone{n_clones}_rectangle0_w{sw:.1f}"
    if result_dir.exists():
        outputs = list(result_dir.glob("*.tsv"))
        logger.info(f"CNA calling complete. Results in {result_dir}:")
        for o in sorted(outputs):
            logger.info(f"  {o.name}")
    else:
        logger.warning(f"Expected result dir not found: {result_dir}")


if __name__ == "__main__":
    main()
