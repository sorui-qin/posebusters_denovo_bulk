"""
Author: Rui Qin
Date: 2023-12-14 17:12:55
LastEditTime: 2023-12-14 17:15:13
Description:
"""
from posebusters import PoseBusters
from pathlib import Path
import pandas as pd
import glob
import sys
import subprocess
from tqdm import tqdm
from multiprocessing import Pool


def run_buster(path) -> pd.DataFrame:
    """
    Execute Posebuster.
    """
    sdf_ls = glob.glob(f"{path}/*out.sdf")
    if sdf_ls == []:
        raise ValueError("Wrong Path")
    pred_file = sdf_ls[:-1]
    buster = PoseBusters(config="mol")
    df = buster.bust(pred_file)
    return df


if __name__ == "__main__":
    input, out = sys.argv[1], sys.argv[2]
    Path(out).mkdir(parents=True, exist_ok=True)

    # pocket_ls = glob.glob(f"{input}/*")
    pocket_ls = glob.glob(f"{input}/*/SDF")
    for pocket in tqdm(pocket_ls):
        p = Path(pocket)

        output_csv = f"{out}/{p.name}.csv"
        if Path(output_csv).exists():
            continue
        df = run_buster(str(p))
        df.to_csv(output_csv)
