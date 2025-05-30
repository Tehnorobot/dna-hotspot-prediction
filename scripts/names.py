import pandas as pd
import re
import os

def create_mapping_file(gtf_path, output_mapping_path):
    """
    Parses a GFF3/GTF file to create a mapping of technical IDs (ID/gene_id/transcript_id)
    to gene names (Name/gene), focusing only on 'gene', 'mRNA', or 'transcript' features.
    
    Args:
        gtf_path (str): Path to the input GFF3/GTF file.
        output_mapping_path (str): Path to save the output TSV mapping file.
    """
    
    if not os.path.exists(gtf_path):
        print(f"Ошибка: Входной файл не найден по пути: {gtf_path}")
        return

    output_dir = os.path.dirname(output_mapping_path)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Создана выходная директория: {output_dir}")
        except OSError as e:
            print(f"Ошибка при создании директории {output_dir}: {e}")
            return

    mapping_data = {}
    print(f"Парсинг файла: {gtf_path}...")
    
    line_count = 0
    records_found = 0
    # debug_limit = 20

    with open(gtf_path, 'r') as f:
        for line in f:
            line_count += 1
            if line.startswith('#'):
                continue 
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                if line_count <= debug_limit:
                    print(f"DEBUG: String {line_count} has less than 9 columns: '{line.strip()}'")
                continue
            
            feature_type = parts[2] # Тype (gene, CDS, exon)
            attributes_str = parts[8] # Attributes

            technical_id = None
            gene_name = None

           
            if feature_type in ['gene', 'mRNA', 'transcript']:
                # GFF3-style regex (key=value;)
                id_match_gff3 = re.search(r'ID=([^;]+)', attributes_str)
                name_match_gff3 = re.search(r'Name=([^;]+)', attributes_str)
                gene_match_gff3 = re.search(r'gene=([^;]+)', attributes_str) # 'gene' attribute in GFF3-style

                # GTF-style regex (key "value";) - kept as fallback
                gene_id_match_gtf = re.search(r'gene_id "([^"]+)"', attributes_str)
                transcript_id_match_gtf = re.search(r'transcript_id "([^"]+)"', attributes_str)
                gene_name_match_gtf = re.search(r'gene "([^"]+)"', attributes_str) # 'gene' attribute in GTF-style

                # Debug messages
                if line_count <= debug_limit:
                    print(f"\n--- DEBUG: Str {line_count} (Type: {feature_type}) ---")
                    print(f"Attribute str: '{attributes_str}'")
                    print(f"  GFF3 ID: {id_match_gff3.group(1) if id_match_gff3 else 'NOT FOUND'}")
                    print(f"  GFF3 Name: {name_match_gff3.group(1) if name_match_gff3 else 'NOT FOUND'}")
                    print(f"  GFF3 gene: {gene_match_gff3.group(1) if gene_match_gff3 else 'NOT FOUND'}")
                    print(f"  GTF gene_id: {gene_id_match_gtf.group(1) if gene_id_match_gtf else 'NOT FOUND'}")
                    print(f"  GTF transcript_id: {transcript_id_match_gtf.group(1) if transcript_id_match_gtf else 'NOT FOUND'}")
                    print(f"  GTF gene name (from 'gene'): {gene_name_match_gtf.group(1) if gene_name_match_gtf else 'NOT FOUND'}")

                # Determine Technical_ID: Prioritize GFF3 'ID', then GTF 'gene_id', then GTF 'transcript_id'
                if id_match_gff3:
                    technical_id = id_match_gff3.group(1)
                elif gene_id_match_gtf:
                    technical_id = gene_id_match_gtf.group(1)
                elif transcript_id_match_gtf and transcript_id_match_gtf.group(1) != "":
                    technical_id = transcript_id_match_gtf.group(1)
                
                # Determine Gene_Name: Prioritize GFF3 'Name', then GFF3 'gene', then GTF 'gene'
                if name_match_gff3:
                    gene_name = name_match_gff3.group(1)
                elif gene_match_gff3:
                    gene_name = gene_match_gff3.group(1)
                elif gene_name_match_gtf:
                    gene_name = gene_name_match_gtf.group(1)
                
                # Fallback: if Technical_ID found but no Gene_Name, use Technical_ID as Gene_Name
                if technical_id and not gene_name:
                    gene_name = technical_id
                    if line_count <= debug_limit:
                        print(f"  -> Gene name not found, using technical gene name: '{technical_id}'")

                if technical_id and gene_name:
                    mapping_data[technical_id] = gene_name
                    records_found += 1
                    if line_count <= debug_limit:
                        print(f"  -> Added to mapping: '{technical_id}' -> '{gene_name}'")
                elif line_count <= debug_limit:
                    print(f"  -> NOT added to mapping (technical_id: {technical_id}, gene_name: {gene_name})")
            else: # If feature_type is not 'gene', 'mRNA', or 'transcript'
                if line_count <= debug_limit:
                    print(f"\n--- DEBUG: String {line_count} (Type: {feature_type}) ---")
                    print(f"  -> SKIPPED (is not 'gene', 'mRNA', 'transcript' type)")

    print(f"Parsed {line_count} strings.")
    print(f"Found {records_found} unique mapping records.")

    mapping_df = pd.DataFrame(mapping_data.items(), columns=['Technical_ID', 'Gene_Name'])
    
    mapping_df.to_csv(output_mapping_path, sep='\t', index=False)
    print(f"Mapping file saved at: {output_mapping_path} with {len(mapping_df)} records.")
    
    if not mapping_df.empty:
        print("Example:")
        print(mapping_df.head())
    else:
        print("Mapping file is empty. Please, look into debug messages.")
        print(f"Check input file path: {gtf_path}")



if __name__ == "__main__":
    gtf_path = "/home/mark/Documents/itmo/ecolibl21/genomic.gff" 
    output_mapping_path = '/home/mark/Documents/itmo/ecolibl21/id_to_gene_name_mapping.tsv'

    create_mapping_file(gtf_path, output_mapping_path)