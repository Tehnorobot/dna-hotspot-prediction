import sys
import os
import argparse
import pandas as pd

DEFAULT_OUTPUT_TSV = "~/pipeline/results/matrix_normal_seq_ori_gc.tsv"
GENOMIC_SEQ_COL = "gene_genomic_seq"

def calculate_gc_content(sequence: str) -> float:
    if not isinstance(sequence, str) or sequence.upper() == "N/A" or not sequence:
        return float('nan') 
    seq_upper = sequence.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total_bases = len(seq_upper)
    if total_bases == 0:
    	return 0.0

    return gc_count / total_bases

def main():
    parser = argparse.ArgumentParser(
		description ="Adds a gc_content column to a .tsv matrix using genomic sequences."
	)
    parser.add_argument(
        "input_matrix_tsv",
        help="Path to the input TSV matrix (must contain a column named 'gene_genomic_seq')."
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_TSV,
        help=f"Name of the output TSV file. Default: {DEFAULT_OUTPUT_TSV}"
    )
    parser.add_argument(
        "--seq_col",
        default=GENOMIC_SEQ_COL,
        help=f"Name of the column containing genomic sequences. Default: '{GENOMIC_SEQ_COL}'"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input_matrix_tsv):
        sys.exit(f"Error: Input TSV matrix file not found: {args.input_matrix_tsv}")

    try:
        df = pd.read_csv(args.input_matrix_tsv, sep='\t')
    except Exception as e:
        sys.exit(f"Error reading input TSV file '{args.input_matrix_tsv}': {e}")

    if args.seq_col not in df.columns:
        sys.exit(f"Error: Column '{args.seq_col}' not found in the input matrix. "
                 f"Please ensure the column name is correct or specify it with --seq_col.")

    print(f"Calculating GC content for sequences in column '{args.seq_col}'...")
    df['gc_content'] = df[args.seq_col].apply(calculate_gc_content)
    
    try:
        df.to_csv(args.output, sep='\t', index=False)
        print(f"Successfully added GC content. Results saved to '{args.output}'")
    except Exception as e:
        sys.exit(f"Error saving results to '{args.output}': {e}")

if __name__ == "__main__":
    main()