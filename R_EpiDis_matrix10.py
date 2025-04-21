#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[Complex Release Version] R_EpiDis_matrix10.py

Version: 1.1.2
Author: Wenxuan Xu (riffleriver@163.com)
Created on Fri Apr  5 02:56:33 2024

Usage:
    python3 R_EpiDis_matrix10.py --snp <dir_snp> --prefix <output_prefix> [--chunksize <int>] [--format <matrix|csv|parquet>]

Description:
    This script performs a transformation on a SNP matrix by converting each row
    into a binary vector, where the entry corresponding to the most frequent value
    is set to 1 and others to 0. The computation is carried out in parallel on data
    chunks.

License:
    Proprietary – internal use only.
"""

import sys
import gc
import os
import logging
import argparse
import multiprocessing
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# -----------------------------------------------------------------------------
def _init_logger():
    logger = logging.getLogger("R_EpiDis_matrix10")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    handler.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(handler)
    return logger

_logger = _init_logger()

# -----------------------------------------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(
        description="SNP Matrix Transformation with parallel processing."
    )
    parser.add_argument("--snp", required=True, help="Input SNP matrix file (tab-delimited or Parquet).")
    parser.add_argument("--prefix", required=True, help="Output file prefix.")
    parser.add_argument("--chunksize", type=int, default=2000,
                        help="Chunk size for parallel processing (default: 2000).")
    parser.add_argument("--format", choices=["matrix", "csv", "parquet"], default="matrix",
                        help="Output file format: matrix (default, uncompressed text), csv (gzip-compressed CSV), or parquet.")
    args = parser.parse_args()
    return args

# -----------------------------------------------------------------------------
def _load_data(snp_path):
    _logger.info("Loading SNP data from: %s", snp_path)
    ext = os.path.splitext(snp_path)[1].lower()
    if ext == ".parquet":
        _logger.info("Input file detected as Parquet. Reading with pd.read_parquet: %s", snp_path)
        df = pd.read_parquet(snp_path)
    else:
        _logger.info("Input file detected as text (tab-delimited). Reading with pd.read_table: %s", snp_path)
        df = pd.read_table(snp_path, index_col=0)
    gc.collect()
    return df

# -----------------------------------------------------------------------------
def _transform_row(row):
    value_counts = row.value_counts()
    max_value = value_counts.idxmax()
    return (row == max_value).astype(int)

# -----------------------------------------------------------------------------
def _apply_transform(chunk):
    return chunk.apply(_transform_row, axis=1)

# -----------------------------------------------------------------------------
def _process_matrix(df, chunksize):
    n_chunks = (len(df) - 1) // chunksize + 1
    _logger.info("Splitting data into %d chunks (chunksize = %d).", n_chunks, chunksize)
    chunks = [df.iloc[i * chunksize : (i + 1) * chunksize] for i in range(n_chunks)]
    
    num_cores = multiprocessing.cpu_count()
    _logger.info("Processing with %d CPU cores.", num_cores)
    results = Parallel(n_jobs=num_cores, backend="loky")(
        delayed(_apply_transform)(chunk) for chunk in tqdm(chunks, desc="Transforming chunks")
    )
    transformed_df = pd.concat(results)
    return transformed_df

# -----------------------------------------------------------------------------
def _save_data(df, output_prefix, out_format):
    if out_format == "parquet":
        output_file = output_prefix + ".10.matrix.parquet"
        _logger.info("Saving transformed data as Parquet to: %s", output_file)
        df.to_parquet(output_file)
    elif out_format == "csv":
        output_file = output_prefix + ".10.matrix.csv.gz"
        _logger.info("Saving transformed data as CSV (gzip compressed) to: %s", output_file)
        df.to_csv(output_file, sep="\t", compression="gzip")
    elif out_format == "matrix":
        output_file = output_prefix + ".10.matrix"
        _logger.info("Saving transformed data as plain text matrix to: %s", output_file)
        df.to_csv(output_file, sep="\t")
    else:
        raise ValueError("Unsupported output format: " + out_format)
    _logger.info("Data saved successfully.")

# -----------------------------------------------------------------------------
def main():
    args = _parse_args()
    df_snp = _load_data(args.snp)
    _logger.info("Starting transformation on %d rows.", df_snp.shape[0])
    
    transformed_matrix = _process_matrix(df_snp, args.chunksize)
    _save_data(transformed_matrix, args.prefix, args.format)
    
    _logger.info("Processing complete.")

if __name__ == "__main__":
    main()
