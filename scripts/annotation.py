import pandas as pd
import os


def annotate_matrix(input_matrix_path, mapping_file_path, output_matrix_path):
    print(f"Uploading comparsion file: {mapping_file_path}...")
     
    try:
        gene_mapping = pd.read_csv(mapping_file_path, sep='\t')
        gene_mapping.set_index('Technical_ID', inplace=True)
        id_to_name_series = gene_mapping['Gene_Name']
    except FileNotFoundError:
        print(f"Comparsion file '{mapping_file_path}' not found")
        return
    print("Uploaded successfully")

    print(f"Loading aggregated .tsv: {input_matrix_path}...")
    try:
        df_matrix = pd.read_csv(input_matrix_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: '{input_matrix_path}' not found")
        return

    print("Uploaded successfully")

    target_id_column_name = df_matrix.columns[0]

    if target_id_column_name not in df_matrix.columns:
        print(f"Error: col '{target_id_column_name}' not found")
        return

    print(f"Comparin ID '{target_id_column_name}' with name...")
    df_matrix['Gene_Name'] = df_matrix[target_id_column_name].map(id_to_name_series)

    if df_matrix['Gene_Name'].isnull().any():
        null_count = df_matrix['Gene_Name'].isnull().sum()
        print(f"Warning: {null_count} ID were not found in comparsion file.")
        print(f"Saving original technical ID.")

    df_matrix['Gene_Name'].fillna(df_matrix[target_id_column_name], inplace=True)
    cols = ['Gene_Name'] + [col for col in df_matrix.columns if col not in ['Gene_Name']]
    df_matrix = df_matrix[cols]

    print(f"Saving annotated matrix: {output_matrix_path}...")
    df_matrix.to_csv(output_matrix_path, sep='\t', index=False)
    print("Annotated!")
    print(f"Example:")
    print(df_matrix.head())

if __name__ == "__main__":
    INPUT_MATRIX_FILE = '~/pipeline/results/normalized_kallisto_tpm_matrix.tsv' 
    MAPPING_FILE = '~/pipeline/results/id_to_gene_name_mapping.tsv'
    OUTPUT_MATRIX_FILE = '~/pipeline/results/aggregated_tpm_matrix_annotated.tsv' 

    annotate_matrix(INPUT_MATRIX_FILE, MAPPING_FILE, OUTPUT_MATRIX_FILE)
