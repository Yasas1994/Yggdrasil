#!/usr/bin/env python3
"""Convert a tab-separated file to parquet (runs inside the vContact3 env, which has pyarrow)."""
import argparse
import sys

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    a = ap.parse_args()
    pd.read_csv(a.input, sep="\t").to_parquet(a.output, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
