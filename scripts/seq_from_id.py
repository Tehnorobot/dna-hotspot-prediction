#!/usr/bin/env python3
import sys
import os
import argparse
from pyfaidx import Fasta

DEFAULT_GENOME_FASTA = "~/pipeline/data/data/genomic.fna"
DEFAULT_GENOME_GFF = "~/pipeline/data/data/genomic.gff"
DEFAULT_OUTPUT_TSV = "~/pipeline/results/normal_matrix_with_seqs.tsv"

def get_attribute_value(attributes_str, key):
    """Exctracting attribute value from GFF."""
    for part in attributes_str.split(';'):
        if part.startswith(f'{key}='):
            return part[len(key)+1:]
    return None

def reverse_complement(seq):
    """Returning reverse complement."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                      'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                      'N': 'N', 'n': 'n'}
    return "".join(complement_map.get(base, base) for base in reversed(seq))

def main():
    parser = argparse.ArgumentParser(
        description="Adding genome sequences to input TSV matrix using GFF and FASTA."
    )
    parser.add_argument(
        "input_matrix_tsv",
        help="Path to input TSV матрице with 'transcript_id' column."
    )
    parser.add_argument(
        "--fasta",
        default=DEFAULT_GENOME_FASTA,
        help=f".fasta path (default: {DEFAULT_GENOME_FASTA})"
    )
    parser.add_argument(
        "--gff",
        default=DEFAULT_GENOME_GFF,
        help=f".gff annotation path (default: {DEFAULT_GENOME_GFF})"
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_TSV,
        help=f"output matrix path (default: {DEFAULT_OUTPUT_TSV})"
    )
    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        sys.exit(f"Error: FASTA not found: {args.fasta}")
    if not os.path.exists(args.gff):
        sys.exit(f"Error: GFF not found: {args.gff}")
    if not os.path.exists(args.input_matrix_tsv):
        sys.exit(f"Error: TSV not found: {args.input_matrix_tsv}")

    try:
        genome = Fasta(args.fasta)
    except Exception as e:
        sys.exit(f"Error loading FASTA: {e}")

   
    requested_gene_ids = set()
    print(f"Reading genes ID from {args.input_matrix_tsv}...")
    try:
        with open(args.input_matrix_tsv, 'r') as f_in:
            header_line = f_in.readline().strip() # Читаем заголовок
            header = header_line.split('\t')
            try:
                gene_id_col_idx = header.index('transcript_id')
            except ValueError:
                sys.exit("Error: 'transcript_id' column not found in input .tsv.")
            
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) > gene_id_col_idx:
                    requested_gene_ids.add(parts[gene_id_col_idx])
        print(f"Found {len(requested_gene_ids)} unique ID's.")
    except Exception as e:
        sys.exit(f"Error reading input TSV: {e}")

   
    gene_coordinates = {} 

    print(f"Parcing GFF: {args.gff} for gene coordinate gathering...")
    with open(args.gff, 'r') as gff_in:
        for line in gff_in:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, _, feature_type, start_str, end_str, _, strand, _, attributes = parts

            
            if feature_type == 'gene':
                gene_id = get_attribute_value(attributes, 'ID')
                if gene_id:
                    try:
                        start = int(start_str)
                        end = int(end_str)
                        gene_coordinates[gene_id] = (seqid, start, end, strand)
                    except ValueError:
                        print(f"Warning: incorrect coordinates for gene {gene_id}: {line.strip()}")
                        pass

    print(f"Processing input matrix and adding genome sequences...")
    processed_rows_count = 0
    with open(args.input_matrix_tsv, 'r') as f_in, open(args.output, 'w') as f_out:
        f_in.seek(0) 
        header_line = f_in.readline().strip()
        header = header_line.split('\t')
        
        f_out.write('\t'.join(header + ['gene_genomic_seq']) + '\n')

        for line in f_in: 
            parts = line.strip().split('\t')
            current_gene_id = None
            if len(parts) > gene_id_col_idx:
                current_gene_id = parts[gene_id_col_idx]
            
            gene_seq = "N/A" 
            if current_gene_id in gene_coordinates:
                seqid, start, end, strand = gene_coordinates[current_gene_id]
                try:
                   
                    extracted_seq = str(genome[seqid][start-1:end]).upper()
                    if strand == '-':
                        extracted_seq = reverse_complement(extracted_seq)
                    gene_seq = extracted_seq
                except KeyError:
                    print(f"Warning: Contig '{seqid}' not found in FASTA for gene {current_gene_id}. Skipping.")
                    pass
                except Exception as e:
                    print(f"Errror exctracting sequence for gene: {current_gene_id}: {e}. Skipping.")
                    pass
            
            f_out.write('\t'.join(parts + [gene_seq]) + '\n')
            processed_rows_count += 1
    
    print(f"Job is done! Processed {processed_rows_count} strings. Result is saved in output file: {args.output}")

if __name__ == "__main__":
    main()

# Please, launch with explicit input path 'python3 seq_from_id.py my_input_matrix.tsv'