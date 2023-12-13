"""
Author: Rui Qin
Date: 2023-12-13 10:45:33
LastEditTime: 2023-12-13 10:54:50
Description:Run Posebusters for the generated molecular conformations.
Codes inspired from https://posebusters.readthedocs.io/en/latest/api_notebook.html
Run the script as: python conf_check.py [input_path] [output_path]
"""
from posebusters import PoseBusters
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool
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


def run_buster(sdf) -> pd.DataFrame:
    """
    Execute Posebuster.
    """
    df = PoseBusters(config="mol").bust(sdf)
    freq = df_metrics_to_freq(df)
    return freq


if __name__ == "__main__":
    input, out = sys.argv[1], sys.argv[2]
    Path(out).mkdir(parents=True, exist_ok=True)

    sdf_paths = Path(input).glob("*")
    for path in sdf_paths:
        output_csv = f"{out}/{Path(path).name}.csv"
        if Path(output_csv).exists():
            continue

        overall_df = pd.DataFrame()
        with Pool() as pool:
            df_list = list(
                tqdm(pool.imap(run_buster, Path(path).glob("*.sdf")), total=len(list(Path(path).glob("*.sdf"))))
            )
        for freq in df_list:
            overall_df = pd.concat([overall_df, freq])
        overall_df.to_csv(output_csv)
