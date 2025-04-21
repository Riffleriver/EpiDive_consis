
# **Bacterial Genome Consistency Pattern Analysis Tutorial**

## **Overview**
Bacterial genome consistency pattern analysis involves extracting and encoding SNP (Single Nucleotide Polymorphism) data from whole genome datasets into binary patterns. These patterns are used to identify consistency and similarity across samples. This tutorial demonstrates how to process large-scale SNP data, perform SNP consistency analysis, and group SNPs based on their similarity patterns using parallel processing, binary encoding, and the Faiss indexing library for fast retrieval.

## **Prerequisites**
Ensure the following are available:
- **Python 3.x**
- **Required libraries**: `numpy`, `pandas`, `faiss`, `scipy`, `joblib`, `tqdm`
- **SNP data** in tab-delimited, CSV, or Parquet format.

## **Step 1: Setup Environment**
Before proceeding, install the necessary libraries:

```bash
conda install numpy pandas faiss-cpu scipy joblib tqdm
```

## **Step 2: Parallel Binary Pattern Encoding**
The input data should be a tab-delimited SNP matrix where:
- **Rows** represent SNPs.
- **Columns** represent samples.
- Each entry is a 0 or 1 representing the presence or absence of an allele.

You can load SNP data from various file formats (tab-delimited, CSV, or Parquet) using the `R_EpiDis_matrix10.py` script. This script applies binary pattern encoding, where the most frequent allele is encoded as `1`, and all others as `0`.

#### Command for loading data:
```bash
python3 R_EpiDis_matrix10.py --snp <input_snp_file> --prefix <output_prefix>
```
- `--snp`: Path to your input SNP matrix file.
- `--prefix`: Prefix for the output files.
- `--format`: Specify output format as `matrix`, `csv`, or `parquet`.
- `--chunksize`: Optional. Define chunk size for parallel processing (default: 2000).

#### Example:
```bash
python3 R_EpiDis_matrix10.py --snp matrix.tsv --prefix output --format matrix --chunksize 1000
```

This command processes the SNP matrix in chunks of 1000 rows, accelerating the transformation by utilizing multiple CPU cores.

## **Step 3: SNP Grouping by Consistency Pattern**
Once the SNP matrix is transformed into binary vectors, SNPs are grouped based on their consistency patterns using MD5 hashing. Identical binary patterns result in the same hash value.

The input data should be a tab-delimited SNP matrix where:
- **Rows** represent SNPs.
- **Columns** represent samples.
- Each entry in the matrix is a `0` or `1`.

#### Example:
The `R_consis_snp.py` script classifies SNPs based on their binary patterns. You can run the script using the following command:

#### Command:
```bash
python3 R_consis_snp.py --input <input_snp_matrix> --output <output_file> --min_size <min_group_size>
```

- `--input`: Path to the input SNP matrix file (in tab-delimited format).
- `--output`: Prefix for the output file containing the SNP classification.
- `--min_size`: Optional. Minimum size for the SNP groups to be included in the output (default is 1).

#### Example:
```bash
python3 R_consis_snp.py --input your_snp_matrix.txt --output snp_classification.txt --min_size 5
```

This command will:
1. Classify SNPs based on their binary patterns.
2. Output the classified SNP groups to `snp_classification.txt`.
3. Only include groups with more than 5 SNPs.

## **Step 4: Querying Specific SNPs (Optional)**
The `R_classify_snps_near_faiss.py` script allows you to query a specific SNP and retrieve its neighboring SNPs that share similar binary patterns within a defined Hamming distance threshold.

#### Command for querying a single SNP:
```bash
python3 R_classify_snps_near_faiss.py --input matrix.tsv --query_snp snp12345 --allowed_diff 5 --output neighbors.txt
```

- `--query_snp`: SNP ID for which you want to find similar SNPs.
- `--allowed_diff`: Hamming distance threshold for similarity.
- `--output`: File where the neighbor SNPs will be saved.

This command retrieves all SNPs whose binary patterns are within a Hamming distance of 5 (not inclusive of 5) from `snp12345`.

## **Step 5: Genome-Wide Similarity Analysis**
For a genome-wide consistency pattern analysis, you can use the `R_classify_snps_near_faiss.py` script to perform a range search and identify all SNPs with similar patterns within a predefined threshold.

#### Command for genome-wide analysis:
```bash
python3 R_classify_snps_near_faiss.py --input matrix.tsv --allowed_diff 7 --output groups.tsv
```

This command performs a range search to identify similar SNPs across the entire dataset and groups them based on similarity.

## **Step 6: Output and Results**
The results are output in the specified format (e.g., matrix, CSV, or Parquet). If querying a single SNP, the output will contain the list of neighboring SNPs. For genome-wide analysis, the output contains clusters of SNPs with similar patterns.

#### Example output formats:
- **Neighbors List**: For the queried SNP, all SNPs with similar patterns within the specified Hamming distance are listed.
- **Cluster Labels**: Each SNP is assigned a cluster label indicating which group it belongs to.

## **Conclusion**
This tutorial demonstrates the process of bacterial genome consistency pattern analysis, utilizing efficient binary encoding, parallel processing, and Faiss-based similarity search to group and analyze SNPs based on their consistency patterns. This approach is scalable, allowing for both localized SNP queries and large-scale genome-wide analysis.
