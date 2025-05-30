#!/usr/bin/env python3
import sys
import os
import argparse
import pandas as pd
import re

DEFAULT_GENOME_GFF = "~/pipeline/data/data/genomic.gff"
DEFAULT_OUTPUT_TSV = "~/pipeline/results/matrix_normal_seq_ori.tsv"

def get_attribute_value(attributes_str, key):
    for part in attributes_str.split(';'):
        if part.startswith(f'{key}='):
            return part[len(key)+1:]
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Adds a 'distance_to_ori' column to a TSV matrix using gene coordinates from GFF."
    )
    parser.add_argument(
        "input_matrix_tsv",
        help="Path to the input TSV matrix with a 'transcript_id' column (containing gene IDs)."
    )
    parser.add_argument(
        "--gff",
        default=DEFAULT_GENOME_GFF,
        help=f"Path to the GFF annotation file. Default: {DEFAULT_GENOME_GFF}"
    )
    parser.add_argument(
        "--ori_pos",
        type=int,
        default=0,
        help="The numerical position of the origin (oriC). Default: 0"
    )
    parser.add_argument(
        "--feature_type",
        default="gene",
        help="The type of feature to extract from GFF (e.g., 'gene', 'CDS'). Default: 'gene'"
    )
    parser.add_argument(
        "--position_type",
        default="start",
        choices=['start', 'end', 'midpoint'],
        help="The position within the feature to use for distance calculation ('start', 'end', or 'midpoint'). Default: 'start'"
    )
    parser.add_argument(
        "--id_attribute_key",
        default="ID",
        help="The key in the GFF 'attributes' column that holds the feature ID (e.g., 'ID', 'gene_id'). Default: 'ID'"
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_TSV,
        help=f"Name of the output TSV file. Default: {DEFAULT_OUTPUT_TSV}"
    )
    args = parser.parse_args()

    if not os.path.exists(args.gff):
        sys.exit(f"Error: GFF file not found: {args.gff}")
    if not os.path.exists(args.input_matrix_tsv):
        sys.exit(f"Error: Input TSV matrix file not found: {args.input_matrix_tsv}")

    gene_coordinates = {} 

    with open(args.gff, 'r') as gff_in:
        for line in gff_in:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, _, feature_type, start_str, end_str, _, strand, _, attributes = parts

            if feature_type == args.feature_type:
                gene_id = get_attribute_value(attributes, args.id_attribute_key)
                if gene_id:
                    try:
                        start = int(start_str)
                        end = int(end_str)
                        gene_coordinates[gene_id] = (seqid, start, end, strand)
                    except ValueError:
                        pass

    processed_rows_count = 0
    with open(args.input_matrix_tsv, 'r') as f_in, open(args.output, 'w') as f_out:
        header_line = f_in.readline().strip()
        header = header_line.split('\t')
        
        try:
            gene_id_col_idx = header.index('transcript_id')
        except ValueError:
            sys.exit("Error: Column 'transcript_id' not found in the input TSV file.")
        
        f_out.write('\t'.join(header + ['distance_to_ori']) + '\n')

        for line in f_in:
            parts = line.strip().split('\t')
            current_gene_id = None
            if len(parts) > gene_id_col_idx:
                current_gene_id = parts[gene_id_col_idx]
            
            distance = "N/A"
            if current_gene_id in gene_coordinates:
                seqid, start, end, strand = gene_coordinates[current_gene_id]
                
                gene_position = 0
                if args.position_type == 'start':
                    gene_position = start
                elif args.position_type == 'end':
                    gene_position = end
                elif args.position_type == 'midpoint':
                    gene_position = (start + end) / 2
                
                distance = gene_position - args.ori_pos
            
            f_out.write('\t'.join(parts + [str(distance)]) + '\n')
            processed_rows_count += 1
    
    print(f"Processed {processed_rows_count} rows. Results saved to '{args.output}'")

if __name__ == "__main__":
    main()

# Please, launch with explicit input path python3 orilen.py my_input_matrix.tsv