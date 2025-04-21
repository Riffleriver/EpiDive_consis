#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[Complex Release Version] R_consis_snp.py

Version: 1.1.2
Author: Wenxuan Xu (riffleriver@163.com)
Created on Apri 20 2025

Usage
    python R_consis_snp.py \
      --input your_snp_matrix.txt \
      --output snp_classification.txt \
      --min_size 5

Description:
    Classify the SNP matrix by grouping SNPs with identical sample 0/1 
    encoding patterns. Then, sort the groups based on the number of SNPs 
    in each group, from largest to smallest, and assign category 
    numbers accordingly. Finally, only output the results for groups 
    where the number of SNPs exceeds a threshold (which can be adjusted 
    using the --min_size parameter).

License:
    Proprietary – internal use only.
"""

import pandas as pd
import numpy as np
import hashlib
from collections import defaultdict
import argparse

def classify_snps_by_pattern(snp_df):
    """
    对 SNP 矩阵进行分类，将具有相同样本 0/1 编码模式的 SNP 归为一类。
    返回一个字典 groups，键为 MD5 模式哈希值，值为 SNP 名称列表。
    """
    groups = defaultdict(list)
    for snp_name, row in snp_df.iterrows():
        row_bytes = row.values.astype(np.uint8).tobytes()
        hash_val = hashlib.md5(row_bytes).hexdigest()
        groups[hash_val].append(snp_name)
    return groups

def main():
    parser = argparse.ArgumentParser(
        description="Classify SNPs by binary pattern and assign ranked group IDs."
    )
    parser.add_argument(
        "--input", required=True,
        help="输入 SNP 矩阵文件（制表符分隔），第一行为表头，第一列为 SNP 名称。"
    )
    parser.add_argument(
        "--output", required=True,
        help="输出文件，每行：SNP_name<tab>group_ID。"
    )
    parser.add_argument(
        "--min_size", type=int, default=0,
        help="只输出那些类别中 SNP 数目大于该阈值的记录（默认 0，输出所有）。"
    )
    args = parser.parse_args()

    # 1. 读取数据
    snp_df = pd.read_csv(args.input, sep="\t", index_col=0)
    print(f"Loaded SNP matrix with shape: {snp_df.shape}")

    # 2. 分类
    groups = classify_snps_by_pattern(snp_df)
    print(f"Total distinct SNP patterns found: {len(groups)}")

    # 3. 按组大小排序并分配类别 ID
    sorted_groups = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)
    snp_to_category = {}
    for category_id, (hash_val, snp_list) in enumerate(sorted_groups, start=1):
        for snp in snp_list:
            snp_to_category[snp] = category_id

    # 4. 构造结果 DataFrame
    result_df = pd.DataFrame({
        "snp": list(snp_to_category.keys()),
        "category": list(snp_to_category.values())
    })
    # 按 category 和 snp 排序
    result_df = result_df.sort_values(["category", "snp"]).reset_index(drop=True)

    # 5. 根据 min_size 参数过滤
    if args.min_size > 0:
        # 计算每个类别的大小
        cat_counts = result_df["category"].value_counts()
        # 保留那些大小 > min_size 的类别
        valid_cats = cat_counts[cat_counts > args.min_size].index
        before = len(result_df)
        result_df = result_df[result_df["category"].isin(valid_cats)]
        after = len(result_df)
        print(f"Filtered by min_size={args.min_size}: {before} → {after} records remain.")

    # 6. 写出结果
    # 无表头，制表符分隔
    result_df.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Classification result saved to: {args.output}")

if __name__ == "__main__":
    main()
