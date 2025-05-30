# main snakefile pipeline
import pandas as pd
import os

os.makedirs("results", exist_ok=True)
os.makedirs("scripts", exist_ok=True)

configfile: "cfg.yaml"

SAMPLES = pd.read_csv(config["samples_file"], header=None)[0].str.strip().tolist()


rule all:
    input:
        expand("results/kallisto_tpm_matrix.tsv")

# transcriptome indexing
# rule build_kallisto_index:
#     input:
#         transcriptome=config["transcriptome_fasta"]
#     output:
#         index=os.path.join(config["kallisto_index_dir"], "kallisto.idx")
#     log:
#         "logs/kallisto_index.log"
#     conda:
#         "envs/kallisto.yaml"
#     shell: 
#         """
#         mkdir -p {config[kallisto_index_dir]} logs
#         kallisto index -i {output.index} {input.transcriptome} \
#                 2> {log}
#         """
        
# # sra fetch + convert
# rule download_sra:
#     output:
#         r1="data/raw/{sample}_1.fastq.gz",
#         r2="data/raw/{sample}_2.fastq.gz"
#     params:
#         sra_id="{sample}"
#     log:
#         "logs/download_sra/{sample}.log"
#     conda:
#         "envs/sra_toolkit.yaml"
#     shell:
#         """
#         mkdir -p data/raw logs/download_sra
#         fastq-dump --gzip --split-files {params.sra_id} -O data/raw 2> {log}
#         """

# fastp qc
rule trim_reads:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
        r1="data/trimmed/{sample}_1.trimmed.fastq.gz",
        r2="data/trimmed/{sample}_2.trimmed.fastq.gz",
        html="qc/fastp/{sample}.fastp.html", # Отчет fastp
        json="qc/fastp/{sample}.fastp.json"  # JSON для MultiQC
    params:
        options=config["fastp_options"]
    log:
        "logs/trim_reads/{sample}.log"
    conda:
        "envs/fastp.yaml"
    threads: 4
    shell:
        """
        mkdir -p data/trimmed qc/fastp logs/trim_reads
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --html {output.html} --json {output.json} \
              {params.options} \
              --thread {threads} 2> {log}
        """

# quantification via kallisto
rule kallisto_quant:
    input:
        r1="data/trimmed/{sample}_1.trimmed.fastq.gz",
        r2="data/trimmed/{sample}_2.trimmed.fastq.gz",
        idx=os.path.join(config["kallisto_index_dir"], "kallisto.idx")
    output:
        abundance="kallisto_results_2/{sample}/abundance.tsv"
    params:
        options=config["kallisto_options"],
        bootstrap=config["kallisto_bootstrap"],
        outdir="kallisto_results_2/{sample}" 
    log:
        "logs/kallisto_quant/{sample}.log"
    conda:
        "envs/kallisto.yaml"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir} logs/kallisto_quant
        kallisto quant -i {input.idx} \
                       -o {params.outdir} \
                       {params.options} \
                       -b {params.bootstrap} \
                       {input.r1} {input.r2} \
                       --threads {threads} \
                       2> {log}
        """

# aggregation
rule aggregate_kallisto_tpm:
    input:
        kallisto_outputs=expand("kallisto_results_2/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        "results/kallisto_tpm_matrix.tsv"
    log:
        "logs/aggregate_kallisto_tpm.log"
    script:
        "scripts/kallisto_aggregate.py"

    #input:
    #    kallisto_outputs=expand("kallisto_results/{sample}/abundance.tsv", sample=SAMPLES)
    #output:
    #    "results/kallisto_tpm_matrix.tsv"
    #log:
    #    "logs/aggregate_kallisto_tpm.log"
    #script:
    #    "scripts/aggregate_kallisto_tpm.py {output} {input.kallisto_outputs}"
    
