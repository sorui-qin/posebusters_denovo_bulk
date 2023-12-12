"""
Author: Rui Qin
Date: 2023-12-12 11:01:29
LastEditTime: 2023-12-12 11:01:55
Description:Run Posebusters for the generated molecular conformations.
Codes inspired from https://posebusters.readthedocs.io/en/latest/api_notebook.html
Run the script as: python conf_check.py [input_path] [output_path]
"""
from posebusters import PoseBusters
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys


def df_metrics_to_freq(df: pd.DataFrame) -> pd.DataFrame:
    """
    Counting the frequency of the conformations that passed the check.
    """
    count = df.apply(pd.Series.value_counts)["True":]
    freq = count.apply(lambda x: x / df.shape[0]).reset_index(drop=True)
    return freq


def run_buster(path) -> pd.DataFrame:
    """
    Execute Posebuster.
    """
    overall_df = pd.DataFrame()
    for sdf in tqdm(Path(path).glob("*.sdf")):
        df = PoseBusters(config="mol").bust(sdf)
        freq = df_metrics_to_freq(df)
        overall_df = pd.concat([overall_df, freq])
    return overall_df


if __name__ == "__main__":
    input, out = sys.argv[1], sys.argv[2]
    Path(out).mkdir(parents=True, exist_ok=True)

    sdf_paths = Path(input).glob("*")
    for path in sdf_paths:
        output_csv = f"{out}/{Path(path).name}.csv"
        if Path(output_csv).exists():
            continue
        df = run_buster(path)
        df.to_csv(output_csv)
