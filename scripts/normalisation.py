import pandas as pd
import numpy as np
import os

INPUT_MATRIX = "~/pipeline/results/kallisto_tpm_matrix.tsv"

OUTPUT_MATRIX = "~/pipeline/results/normalized_kallisto_tpm_matrix.tsv"


def main():

    if not os.path.exists(INPUT_MATRIX):
        raise FileNotFoundError(f"Matrix file not found: {INPUT_MATRIX}")

    df = pd.read_csv(INPUT_MATRIX, sep="\t", index_col=0)
    df = df.astype(float)
    print(f"[1] Loaded matrix {df.shape[0]} × {df.shape[1]}")

    # 1) лог2(TPM + 1)
    print("[2] Applying log2(TPM + 1)…")
    df = np.log2(df + 1)

    # 2) Z-нормировка по генам (каждый ген центрируется и масштабируется)
    print("[3] Applying gene Z-normalisation…")
    gene_means = df.mean(axis=1)
    gene_stds  = df.std(axis=1).replace(0, np.nan)
    df = df.sub(gene_means, axis=0).div(gene_stds, axis=0)

    # Сохраняем результат
    out_dir = os.path.dirname(OUTPUT_MATRIX)
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(OUTPUT_MATRIX, sep="\t")
    print(f"[4] Saved: {OUTPUT_MATRIX}")

if __name__=="__main__":
    main()
