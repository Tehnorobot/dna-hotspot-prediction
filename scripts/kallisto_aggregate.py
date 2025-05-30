import pandas as pd

# snakemake — это уже готовый объект, его не надо импортировать
input_files = snakemake.input.kallisto_outputs
output_file = snakemake.output[0]

all_tpm = {}
for f in input_files:
    sample = f.split("/")[-2]
    df = pd.read_csv(f, sep="\t", index_col="target_id")
    all_tpm[sample] = df["tpm"]

combined = pd.DataFrame(all_tpm)
combined.index.name = "transcript_id"
combined.to_csv(output_file, sep="\t")
print(f"Aggregated {len(input_files)} samples into {output_file}")