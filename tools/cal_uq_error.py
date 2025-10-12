#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from glob import glob

# Paths
# studies = ["UQ_analysis"]
studies = ["UQ_OH_analysis_noise_5", "UQ_OH_analysis_noise_10", "UQ_NH_analysis_noise_0","UQ_NH_analysis_noise_5","UQ_NH_analysis_noise_10"]

# cases = ["N_0","N_08","N_16","N_24","N_30","AH64_K80","AH64_K120","AH64_K140", "AH46", "AH64"]
cases = ["N_30"]

for study in studies:
    for case in cases: 
        truth_path = f"SIM_results/NH3_NP_CF_reduce/{case}.csv"
        recon_dir  = f"figs/NH3_NP_CF/{study}/{case}/outputs"
        out_dir    = f"figs/NH3_NP_CF/{study}/{case}/error_outputs"

        os.makedirs(out_dir, exist_ok=True)

        # 1) Load truth DataFrame, index on 'grid'
        truth = pd.read_csv(truth_path)
        if 'grid' not in truth.columns:
            raise KeyError("Truth file must contain 'grid' column")
        truth = truth.set_index('grid')

        # Columns to exclude from error computation
        exclude = {
            'velocity','spread_rate','lambda',
            'density','MF','Y_AR','normalized_grid'
        }

        # 2) Loop over each reconstructed‐output CSV
        for recon_path in sorted(glob(os.path.join(recon_dir, "*.csv"))):
            try:
                recon = pd.read_csv(recon_path)
            except pd.errors.EmptyDataError:
                print(f"Warning: {recon_path} is empty. Skipping...")
                recon = pd.DataFrame()  # 创建一个空的 DataFrame 作为替代
        
            # Must have a 'grid' column to align with truth
            if 'grid' not in recon.columns:
                print(f"Skipping {recon_path}: no 'grid' column")
                continue

            # 3) Columns to compute errors on:
            #    keep 'grid' separately, then any overlapping columns minus exclude
            common_cols = [
                c for c in recon.columns 
                if c in truth.columns and c not in exclude
            ]

            # Build recon_sub (indexed by grid) and truth_sub (same index, columns=common_cols)
            recon_sub = recon.set_index('grid')[common_cols]
            truth_sub = truth[common_cols].reindex(recon_sub.index)

            # 4) Compute relative error DataFrame
            error_df = (recon_sub - truth_sub) / truth_sub

            # 5) Reset index to get 'grid' back as a column and write out
            out_path = os.path.join(out_dir, f"error_{os.path.basename(recon_path)}")
            error_df.reset_index().to_csv(out_path, index=False)
            print("Wrote", out_path)