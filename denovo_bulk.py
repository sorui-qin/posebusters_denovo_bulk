"""
Author: Rui Qin
Date: 2023-12-11 14:03:58
LastEditTime: 2023-12-11 16:57:12
Description: Run Posebusters with bulk of de-novo generated molecules.
Codes inspired from https://posebusters.readthedocs.io/en/latest/api_notebook.html
Run the script as: python denovo_bulk.py [input_path] [output_path]
"""
from posebusters import PoseBusters
from pathlib import Path
import pandas as pd
import glob
import sys
from tqdm import tqdm


def run_buster(path):
    sdf_ls = glob.glob(f"{path}/*.sdf")
    if sdf_ls == []:
        raise ValueError("Wrong Path")
    sdf_ls = sorted(sdf_ls, key=len)
    pred_file, true_file = sdf_ls[:-1], sdf_ls[-1]
    cond_file = glob.glob(f"{path}/*.pdb")[0]

    buster = PoseBusters(config="dock")
    df = buster.bust(pred_file, true_file, cond_file, full_report=True)
    return df


if __name__ == "__main__":
    input, out = sys.argv[1], sys.argv[2]
    Path(out).mkdir(parents=True, exist_ok=True)

    pocket_ls = glob.glob(f"{input}/*")
    for pocket in tqdm(pocket_ls):
        df = run_buster(pocket)
        df.to_csv(f"{out}/{Path(pocket).name}.csv")
