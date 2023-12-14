"""
Author: Rui Qin
Date: 2023-12-11 14:03:58
LastEditTime: 2023-12-14 15:52:10
Description: Run Posebusters with the bulk of de-novo generated molecules.
Codes inspired from https://posebusters.readthedocs.io/en/latest/api_notebook.html
Run the script as: python pocket_denovo.py [input_path] [output_path]
"""
from posebusters import PoseBusters
from pathlib import Path
import pandas as pd
import glob
import sys
import subprocess
from tqdm import tqdm
from multiprocessing import Pool


def pdbqt_to_sdf(pdbqt_path):
    """
    Transform the top1 conformations in pdbqt file (docking results)
    to the sdf file.
    """
    pdbqt_path = str(pdbqt_path)
    top_pdbqt = pdbqt_path.replace(".pdbqt", "_top1.pdbqt")
    top_sdf = pdbqt_path.replace(".pdbqt", "_top1.sdf")
    with open(pdbqt_path, "r") as i, open(top_pdbqt, "w") as o:
        for line in i:
            o.write(line)
            if "ENDMDL" in line:
                break

    command = f"obabel {top_pdbqt} \
                    -o sdf \
                    -O {top_sdf}"
    subprocess.run(args=command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    Path(top_pdbqt).unlink()


def run_buster(path) -> pd.DataFrame:
    """
    Execute Posebuster.
    """
    sdf_ls = glob.glob(f"{path}/*top1.sdf")
    if sdf_ls == []:
        raise ValueError("Wrong Path")
    pred_file, true_file = sdf_ls[:-1], sdf_ls[-1]
    cond_file = glob.glob(f"{path}/*.pdb")[0]
    buster = PoseBusters(config="dock")
    df = buster.bust(pred_file, true_file, cond_file, full_report=True)
    # if need a concise report, set full_report=False
    return df


if __name__ == "__main__":
    input, out = sys.argv[1], sys.argv[2]
    Path(out).mkdir(parents=True, exist_ok=True)

    pocket_ls = glob.glob(f"{input}/*")
    for pocket in tqdm(pocket_ls):
        p = Path(pocket)
        pdbqt_ls = p.glob("*_out.pdbqt")
        with Pool() as pool:
            pool.map(pdbqt_to_sdf, pdbqt_ls)

        output_csv = f"{out}/{p.name}.csv"
        if Path(output_csv).exists():
            continue
        df = run_buster(str(p))
        df.to_csv(output_csv)
