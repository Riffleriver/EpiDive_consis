#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[Complex Release Version] R_classify_snps_near_faiss.py

Version: 1.0.0
Author: Wenxuan Xu (riffleriver@163.com)
Created on April 20, 2025

Usage:
    python R_classify_snps_near_faiss.py \
      --input <input_snp_matrix> \
      --allowed_diff <allowed_hamming_distance> \
      --query_snp <snp_id> \
      --output <output_file>

Description:
    This script performs SNP grouping and similarity analysis using Faiss and 
    SciPy connected_components. It supports querying a single SNP for its 
    neighbors (within a given Hamming distance threshold), or performing 
    a genome-wide analysis to group SNPs with similar patterns.

    Functionality:
    1. Preprocess SNP data using packbits to pad bit-length to a multiple of 8.
    2. Construct Faiss IndexBinaryFlat for fast retrieval of similar SNP patterns.
    3. If --query_snp is specified, the script performs a range search on the 
       single SNP and outputs the neighbors list before exiting.
    4. Otherwise, the script performs a global range search, constructs an 
       adjacency matrix, and computes connected components.
    5. The connected components are re-labeled in descending order of group size, 
       and the SNP group labels are output.

License:
    Proprietary – internal use only.

Example Usage:
    # Query a single SNP's neighbors
    python R_classify_snps_near_faiss.py \
      --input matrix.tsv \
      --allowed_diff 5 \
      --query_snp snp12345 \
      --output neighbors.txt

    # Perform genome-wide grouping
    python R_classify_snps_near_faiss.py \
      --input matrix.tsv \
      --allowed_diff 7 \
      --output groups.tsv
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd
import faiss
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from collections import Counter

def build_index(X, bit_len):
    Xb = np.packbits(X, axis=1) 
    index = faiss.IndexBinaryFlat(bit_len)
    index.add(Xb)
    return index, Xb

def query_single(df, Xb, index, allowed_diff, query_snp, output_path=None):
    snp_name = str(query_snp)
    if snp_name not in df.index:
        logging.error(f"SNP '{snp_name}' not found in input matrix")
        sys.exit(1)
    q_idx = df.index.get_loc(snp_name)
    logging.info(f"Querying neighbors of '{snp_name}' (index {q_idx}) with diff ≤ {allowed_diff}")
    lims, D, I = index.range_search(Xb[q_idx:q_idx+1], allowed_diff)
    start, end = int(lims[0]), int(lims[1])
    neighbors = [int(j) for j in I[start:end] if int(j) != q_idx]
    neighbor_names = df.index[neighbors].tolist()
    logging.info(f"Found {len(neighbor_names)} neighbors for '{snp_name}'")
    if output_path:
        logging.info(f"Writing neighbors to {output_path}")
        with open(output_path, 'w') as fout:
            for name in neighbor_names:
                fout.write(name + '\n')
    else:
        for name in neighbor_names:
            print(name)
    sys.exit(0)

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S"
    )
    parser = argparse.ArgumentParser(description="Faiss SNP grouping & single-SNP neighbor query")
    parser.add_argument("-i", "--input",      required=True,
                        help="Input SNP matrix TSV: The first column is SNP names, followed by columns with 0/1 encoding")
    parser.add_argument("-o", "--output",     required=True,
                        help="Output path: In query mode, output the neighbor list; in global mode, output the grouping results")
    parser.add_argument("-d", "--allowed_diff", type=int, default=2,
                        help="Allowed mismatch count (Hamming distance threshold), default is 2")
    parser.add_argument("-q", "--query_snp",
                        help="Optional: Query all neighbors of the specified SNP, skipping global clustering")
    args = parser.parse_args()
    logging.info(f"Loading matrix from {args.input}")
    df = pd.read_csv(args.input, sep="\t", index_col=0, dtype=str)
    df.index = df.index.astype(str)
    n_snps, n_bits = df.shape
    logging.info(f"Matrix shape: {n_snps} SNPs × {n_bits} bits")
    X = df.values.astype(np.uint8)
    pad = (-n_bits) % 8
    if pad:
        logging.info(f"Padding {pad} zero-columns to reach {n_bits+pad} bits")
        X = np.hstack([X, np.zeros((n_snps, pad), dtype=np.uint8)])
        n_bits += pad
    index, Xb = build_index(X, n_bits)
    if args.query_snp:
        query_single(df, Xb, index, args.allowed_diff, args.query_snp, args.output)

    logging.info(f"Performing global range_search with diff ≤ {args.allowed_diff}")
    lims, D, I = index.range_search(Xb, args.allowed_diff)
    logging.info("Global range_search complete")
    lengths = np.diff(lims).astype(np.intp)
    rows = np.repeat(np.arange(n_snps, dtype=np.intp), lengths)
    cols = I.astype(np.intp)
    mask = rows != cols
    rows, cols = rows[mask], cols[mask]
    logging.info(f"Collected {rows.size} edges (excluding self-loops)")
    data = np.ones(rows.size, dtype=bool)
    adj = coo_matrix((data, (rows, cols)), shape=(n_snps, n_snps))
    adj = adj + adj.T
    logging.info("Computing connected components via SciPy")
    n_comp, labels0 = connected_components(adj, directed=False, return_labels=True)
    logging.info(f"Found {n_comp} components")
    cnt = Counter(labels0)
    order = [lab for lab, _ in cnt.most_common()]
    label_map = {lab: idx+1 for idx, lab in enumerate(order)}
    labels = np.array([label_map[lab] for lab in labels0], dtype=np.int64)
    logging.info(f"Writing group labels to {args.output}")
    out_df = pd.DataFrame({"snp": df.index, "group": labels})
    out_df.sort_values(by=["group", "snp"]).to_csv(args.output,
                                                  sep="\t", index=False, header=False)

    logging.info("Done.")

if __name__ == "__main__":
    main()
