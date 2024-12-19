# DESCRIPTION: Snakemake workflow for comparative genome-centric metagenomics
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

import os
import sys
import math
shell.executable("/bin/bash")

def get_input(sample,subset_gb,datatype):
    if datatype == "reads": return ''.join(["reads/", str(sample), "_", str(subset_gb), "gb.fastq.gz"])
    else: return ''.join(["contigs/", str(sample), "_", str(subset_gb), ".fasta"])

def get_input_fa(sample,subset_gb,datatype):
    if datatype == "reads": return ''.join(["reads/", str(sample), "_", str(subset_gb), "gb.fasta"])
    else: return ''.join(["contigs/", str(sample), "_", str(subset_gb), ".fasta"])

def get_mem(yield_gb,attempt,multi,base_mb):
    return (base_mb+multi*yield_gb*1000)*((3+attempt)/4)

def get_partition(mem_mb):
    if mem_mb < 500000: return "general"
    else: return "high-mem"

onstart:
    from snakemake.utils import min_version
    min_version("8.0.0")
    if not os.path.exists("results"): os.makedirs("results")
    if not os.path.exists("reads"): os.makedirs("reads")
    if not os.path.exists("contigs"): os.makedirs("contigs")
    if not os.path.exists("assembly"): os.makedirs("assembly")
    if not os.path.exists("tiara"): os.makedirs("tiara")
    if not os.path.exists("kaiju"): os.makedirs("kaiju")
    if not os.path.exists("singlem"): os.makedirs("singlem")
    if not os.path.exists("melon"): os.makedirs("melon")
    if not os.path.exists("barrnap"): os.makedirs("barrnap")
    if not os.path.exists("jellyfish"): os.makedirs("jellyfish")
    if not os.path.exists("longshot"): os.makedirs("longshot")
    ruleorder: barrnap_pred > longshot

rule agr_all:
    input:
        mmlong2="results/mmlong2_mags.tsv",
        assembly="results/assembly.tsv",
        jellyfish="results/jellyfish_sum.tsv",
        cramino="results/cramino.tsv",
        tiara_r="results/tiara_reads.tsv",
        tiara_c="results/tiara_contigs.tsv",
        kaiju_r="results/kaiju_reads.tsv",
        kaiju_c="results/kaiju_contigs.tsv",
        singlem_r="results/singlem_reads.tsv",
        singlem_c="results/singlem_contigs.tsv",
        melon_r="results/melon_reads.tsv",
        melon_c="results/melon_contigs.tsv",
        barrnap_r="results/barrnap_reads.tsv",
        barrnap_c="results/barrnap_contigs.tsv",
        longshot="results/longshot_contigs.tsv",
    output: "done.txt"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        touch {output}
        """

rule rasusa:
    input: expand("{read_dir}/{{sample}}_R1041_trim_filt.fastq.gz",read_dir=config["read_dir"])
    output: "reads/{sample}_{subset_gb}gb.fastq.gz"
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_rasusa"],
        subset="{subset_gb}",
    threads: 1
    resources:
        runtime="2d",
        mem_mb=50*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        rasusa reads -O g -b {params.subset}GB -o {output} {input}
        conda deactivate
        """

rule mmlong2:
    input: "reads/{sample}_{subset_gb}gb.fastq.gz"
    output: "{sample}_{subset_gb}/results/{sample}_{subset_gb}_bins.tsv"
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_mmlong2"],
        sample="{sample}",
        subset="{subset_gb}",
        tmp=config["tmp_dir"],
    threads: 60
    resources:
        runtime="14d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),8,75000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),8,75000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        if [ ! -d "{params.sample}_{params.subset}" ]; then mkdir {params.sample}_{params.subset}; fi
        if [ ! -d "{params.sample}_{params.subset}/results" ]; then mkdir {params.sample}_{params.subset}/results; fi
        if [ -d "{params.sample}_{params.subset}/tmp/polishing" ]; then if [ ! -f "{params.sample}_{params.subset}/tmp/polishing/asm_pol.fasta" ]; then rm -r {params.sample}_{params.subset}/tmp/polishing && rm {params.sample}_{params.subset}/tmp/flye/assembly.fasta.map-ont.mmi; fi; fi
        mmlong2-lite -np {input} -o {params.sample}_{params.subset} -p {threads} -tmp {params.tmp} -sem soil
        conda deactivate
        """

rule mmlong2_agr:
    input: expand("{sample}_{subset_gb}/results/{sample}_{subset_gb}_bins.tsv",sample=config["sample"],subset_gb=config["subset_gb"])
    output: "results/mmlong2_mags.tsv"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        head -n 1 {input[0]} > {output}
        for file in {input}; do tail -n +2 $file >> {output}; done
        """

rule reads_prep:
    input: "reads/{sample}_{subset_gb}gb.fastq.gz",
    output: "reads/{sample}_{subset_gb}gb.fasta",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_barrnap"],
    threads: 10
    resources:
        runtime="6d",
        mem_mb=50*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        seqkit fq2fa -j {threads} {input} > {output}
        conda deactivate
        """

rule assembly_prep:
    input: "{sample}_{subset_gb}/results/{sample}_{subset_gb}_bins.tsv",
    output: "contigs/{sample}_{subset_gb}.fasta",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=1000,
        slurm_partition="general",
    shell:
        """
        cp {wildcards.sample}_{wildcards.subset_gb}/tmp/polishing/asm_pol_lenfilt.fasta {output}
        """
	
rule assembly_stats:
    input: "contigs/{sample}_{subset_gb}.fasta",   
    output:
        df1="assembly/nanoq_contigs_{sample}_{subset_gb}gb.tsv",
        df2="assembly/stats_contigs_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_nanoq"],
    threads: 1
    resources:
        runtime="1d",
        mem_mb=10*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        cat {input} | nanoq -s -f --report {output.df1}
        name="{wildcards.sample}_{wildcards.subset_gb}"
        asm_n=$(awk '{{print $1}}' {output.df1})
        asm_bp=$(awk '{{print $2}}' {output.df1})
        asm_n50=$(awk '{{print $3}}' {output.df1})
        echo -e "${{name}}\t${{asm_n}}\t${{asm_bp}}\t${{asm_n50}}" > {output.df2}
        conda deactivate
        """

rule assembly_agr:
    input: expand("assembly/stats_contigs_{sample}_{subset_gb}gb.tsv",sample=config["sample"],subset_gb=config["subset_gb"])
    output: "results/assembly.tsv"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        echo -e "sample\tcontig_n\tcontig_bp\tcontig_n50" > {output}
        cat {input} >> {output}
        """

rule cramino:
    input: "contigs/{sample}_{subset_gb}.fasta",   
    output:
        df1="assembly/cramino_{sample}_{subset_gb}gb.tsv",
        df2="assembly/stats_map_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_cramino"],
    threads: 20
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,20000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,20000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        cramino --threads {threads} {wildcards.sample}_{wildcards.subset_gb}/tmp/binning/mapping_tmp/1_cov.bam > {output.df1}
        awk 'BEGIN {{FS = "\t"}} {{col1[NR] = $1; col2[NR] = $2}} END {{print "{wildcards.sample}_{wildcards.subset_gb}\t" col2[4] "\t" col2[5] "\t" col2[11] "\t" col2[12] "\t" col2[13];}}' {output.df1} > {output.df2}
        conda deactivate
        """

rule cramino_agr:
    input: expand("assembly/stats_map_{sample}_{subset_gb}gb.tsv",sample=config["sample"],subset_gb=config["subset_gb"])
    output: "results/cramino.tsv"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        echo -e "sample\tmap_yield_gb\tmap_cov_mean\tmap_ident_med\tmap_ident_mean\tmap_ident_modal" > {output}
        cat {input} >> {output}
        """

rule kaiju_classify:
    input: lambda wildcards: get_input(wildcards.sample,wildcards.subset_gb,wildcards.datatype)
    output: 
        kj1="kaiju/kaiju_{datatype}_{sample}_{subset_gb}gb_innit.tsv",
        kj2="kaiju/kaiju_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_kaiju"],
        db=config["db_kaiju"],
    threads: 60
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,250000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,250000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        kaiju -z {threads} -t {params.db}/nodes.dmp -f {params.db}/*.fmi -i {input} -o {output.kj1} 
        kaiju-addTaxonNames -r superkingdom -t {params.db}/nodes.dmp -n {params.db}/names.dmp -i {output.kj1} -o {output.kj2}
        conda deactivate
        """

rule kaiju_nonprok:
    input:
        data=lambda wildcards: get_input(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
        kj="kaiju/kaiju_{datatype}_{sample}_{subset_gb}gb.tsv",
    output:
        id="kaiju/kaiju_{datatype}_{sample}_{subset_gb}gb_id.tsv",
        df="kaiju/nanoq_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_nanoq"],
    threads: 10
    resources:
        runtime="6d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,50000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        grep -v -e 'Bacteria' -e 'Archaea' {input.kj} | cut -f2 > {output.id}
        seqkit grep -f {output.id} {input.data} | nanoq -s -f --report {output.df}
        conda deactivate
        """

rule kaiju_sum:
    input: "kaiju/nanoq_{datatype}_{sample}_{subset_gb}gb.tsv",   
    output: "kaiju/stats_{datatype}_{sample}_{subset_gb}gb.tsv"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        name="{wildcards.sample}_{wildcards.subset_gb}"
        nonprok_n=$(awk '{{print $1}}' {input})
        nonprok_bp=$(awk '{{print $2}}' {input})
        nonprok_n50=$(awk '{{print $3}}' {input})
        echo -e "${{name}}\t${{nonprok_n}}\t${{nonprok_bp}}\t${{nonprok_n50}}" > {output}
        """

rule kaiju_agr:
    input: expand("kaiju/stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"])
    output:
        reads="results/kaiju_reads.tsv",
        contigs="results/kaiju_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        samples=({input})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tkaiju_read_n\tkaiju_read_bp\tkaiju_read_n50" > {output.reads}
        cat $reads >> {output.reads}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tkaiju_contig_n\tkaiju_contig_bp\tkaiju_contig_n50" > {output.contigs}
        cat $contigs >> {output.contigs}
        """

rule barrnap_pred:
    input: lambda wildcards: get_input_fa(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
    output: "barrnap/rrna_all_{datatype}_{sample}_{subset_gb}gb.fasta",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_barrnap"],
    threads: 25
    resources:
        runtime="6d",
        mem_mb=200*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        barrnap {input} --threads {threads} --outseq {output}
        conda deactivate
        """

rule barrnap_sub:
    input: "barrnap/rrna_all_{datatype}_{sample}_{subset_gb}gb.fasta",
    output: "barrnap/rrna_16s_{datatype}_{sample}_{subset_gb}gb.fasta",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_barrnap"],
    threads: 10
    resources:
        runtime="6d",
        mem_mb=50*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        grep -A1 ">16S" {input} | seqkit seq -m 1000 > {output}
        conda deactivate
        """

rule barrnap_clust:
    input: "barrnap/rrna_16s_{datatype}_{sample}_{subset_gb}gb.fasta",
    output:
        seq="barrnap/rrna_clust_{datatype}_{sample}_{subset_gb}gb.fasta",
        df="barrnap/rrna_clust_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_barrnap"],
        id = lambda wildcards: "0.97" if wildcards.datatype == "reads" else "0.987" if wildcards.datatype == "contigs" else None,
    threads: 50
    resources:
        runtime="6d",
        mem_mb=100*1000,
        slurm_partition="general",
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        vsearch --cluster_fast {input} --id {params.id} --iddef 0 --strand both --minseqlength 1000 --threads {threads} --centroids {output.seq} --otutabout {output.df}
        conda deactivate
        """

rule barrnap_sum:
    input: "barrnap/rrna_clust_{datatype}_{sample}_{subset_gb}gb.tsv"
    output:
        sum="barrnap/rrna_stats_{datatype}_{sample}_{subset_gb}gb.tsv",
        otu="barrnap/barrnap_otu_{datatype}_{sample}_{subset_gb}gb.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        name="{wildcards.sample}_{wildcards.subset_gb}"
        clust_uniq=$(wc -l {input} | awk '{{print $1-1}}')
        clust_total=$(awk 'NR > 1 {{sum += $2}} END {{print sum}}' {input})
        echo -e "${{name}}\t${{clust_uniq}}\t${{clust_total}}" > {output.sum}
        
        sed 1d {input} | sed 's/^/{wildcards.sample}_{wildcards.subset_gb}\t/' > {output.otu}
        """

rule barrnap_agr:
    input:
        sum=expand("barrnap/rrna_stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"]),
        otu=expand("barrnap/barrnap_otu_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"]),
    output:
        reads_sum="results/barrnap_reads.tsv",
        contigs_sum="results/barrnap_contigs.tsv",
        reads_otu="results/barrnap_otu_reads.tsv",
        contigs_otu="results/barrnap_otu_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        samples=({input.sum})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tbarrnap_read_otu_uniq\tbarrnap_read_otu_all" > {output.reads_sum}
        cat $reads >> {output.reads_sum}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tbarrnap_contig_otu_uniq\tbarrnap_contig_otu_all" > {output.contigs_sum}
        cat $contigs >> {output.contigs_sum}
        
        
        samples=({input.otu})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tbarrnap_read_otu\tbarrnap_read_abund" > {output.reads_otu}
        cat $reads >> {output.reads_otu}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tbarrnap_contig_otu\tbarrnap_contig_abund" > {output.contigs_otu}
        cat $contigs >> {output.contigs_otu}
        """

rule melon:
    input: lambda wildcards: get_input(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
    output: directory("melon/{datatype}_{sample}_{subset_gb}"),
    wildcard_constraints: datatype="reads|contigs",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_melon"],
        db=config["db_melon"],
    threads: 30
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,50000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        melon {input} -d {params.db} -o {output} -t {threads} 
        conda deactivate
        """

rule melon_sum:
    input: "melon/{datatype}_{sample}_{subset_gb}"
    output:
        sum="melon/stats_{datatype}_{sample}_{subset_gb}gb.tsv",
        otu="melon/melon_otu_{datatype}_{sample}_{subset_gb}gb.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        name="{wildcards.sample}_{wildcards.subset_gb}"
        clust_uniq=$(sed 1d {input}/*.tsv | grep -v "unclassified" | wc -l)
        clust_copy=$(sed 1d {input}/*.tsv | grep -v "unclassified" | awk '{{sum += $9}} END {{print sum}}')
        clust_abund=$(sed 1d {input}/*.tsv | grep -v "unclassified" | awk '{{sum += $10}} END {{print sum}}')
        echo -e "${{name}}\t${{clust_uniq}}\t${{clust_copy}}\t${{clust_abund}}" > {output.sum}
        
        sed 1d {input}/*.tsv | sed 's/^/{wildcards.sample}_{wildcards.subset_gb}\t/' > {output.otu}
        """

rule melon_agr:
    input:
        sum=expand("melon/stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"]),
        otu=expand("melon/melon_otu_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"]),
    output:
        reads_sum="results/melon_reads.tsv",
        contigs_sum="results/melon_contigs.tsv",
        reads_otu="results/melon_otu_reads.tsv",
        contigs_otu="results/melon_otu_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        samples=({input.sum})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tmelon_read_otu_uniq\tmelon_read_otu_all\tmelon_read_otu_abund" > {output.reads_sum}
        cat $reads >> {output.reads_sum}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tmelon_contig_otu_uniq\tmelon_contig_otu_all\tmelon_contig_otu_abund" > {output.contigs_sum}
        cat $contigs >> {output.contigs_sum}
        
        
        samples=({input.otu})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tmelon_read_superkingdom\tmelon_read_phylum\tmelon_read_class\tmelon_read_order\tmelon_read_family\tmelon_read_genus\tmelon_read_species\tmelon_read_copy\tmelon_read_abundance\tmelon_read_identity" > {output.reads_otu}
        cat $reads >> {output.reads_otu}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tmelon_contig_superkingdom\tmelon_contig_phylum\tmelon_contig_class\tmelon_contig_order\tmelon_contig_family\tmelon_contig_genus\tmelon_contig_species\tmelon_contig_copy\tmelon_contig_abundance\tmelon_contig_identity" > {output.contigs_otu}
        cat $contigs >> {output.contigs_otu}
        """

rule singlem:
    input: lambda wildcards: get_input(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
    output:
        prof="singlem/singlem_prof_{datatype}_{sample}_{subset_gb}gb.tsv",
        otu="singlem/singlem_otu_{datatype}_{sample}_{subset_gb}gb.tsv",
        mf="singlem/singlem_mf_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_singlem"],
        db=config["db_singlem"],
    threads: 30
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,50000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        export SINGLEM_METAPACKAGE_PATH={params.db}
        singlem pipe --reads {input} -p {output.prof} --otu-table {output.otu} --threads {threads} --hmmsearch-package-assignment
        singlem microbial_fraction --reads {input} -p {output.prof} > {output.mf} 
        conda deactivate
        """

rule singlem_sum:
    input:
        otu="singlem/singlem_otu_{datatype}_{sample}_{subset_gb}gb.tsv",
        mf="singlem/singlem_mf_{datatype}_{sample}_{subset_gb}gb.tsv",     
    output:
        stats="singlem/singlem_stats_{datatype}_{sample}_{subset_gb}gb.tsv",
        otu="singlem/singlem_otu_{datatype}_{sample}_{subset_gb}gb_id.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        name="{wildcards.sample}_{wildcards.subset_gb}"
        clust_uniq=$(grep "S3.5.ribosomal_protein_S2_rpsB" {input.otu} | wc -l)
        clust_perc=$(awk 'NR > 1 {{print $4}}' {input.mf})
        clust_abund=$(awk 'NR > 1 {{print $2}}' {input.mf})
        echo -e "${{name}}\t${{clust_uniq}}\t${{clust_perc}}\t${{clust_abund}}" > {output.stats}
        
        grep "S3.5.ribosomal_protein_S2_rpsB" {input.otu} | sed 's/^/{wildcards.sample}_{wildcards.subset_gb}\t/' > {output.otu}
        """

rule singlem_agr:
    input:
        sum=expand("singlem/singlem_stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"]),
        otu=expand("singlem/singlem_otu_{datatype}_{sample}_{subset_gb}gb_id.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"])
    output:
        reads_sum="results/singlem_reads.tsv",
        contigs_sum="results/singlem_contigs.tsv",
        reads_otu="results/singlem_otu_reads.tsv",
        contigs_otu="results/singlem_otu_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        samples=({input.sum})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tsinglem_read_otu\tsinglem_read_mf\tsinglem_read_mf_bp" > {output.reads_sum}
        cat $reads >> {output.reads_sum}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tsinglem_contig_otu\tsinglem_contig_mf\tsinglem_contig_mf_bp" > {output.contigs_sum}
        cat $contigs >> {output.contigs_sum}
        
        
        samples=({input.otu})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\tsinglem_read_gene\tsinglem_read_sample\tsinglem_read_sequence\tsinglem_read_num_hits\tsinglem_read_coverage\tsinglem_read_taxonomy" > {output.reads_otu}
        cat $reads >> {output.reads_otu}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\tsinglem_contig_gene\tsinglem_contig_sample\tsinglem_contig_sequence\tsinglem_contig_num_hits\tsinglem_contig_coverage\tsinglem_contig_taxonomy" > {output.contigs_otu}
        cat $contigs >> {output.contigs_otu}
        """

rule tiara:
    input: lambda wildcards: get_input_fa(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
    output: "tiara/tiara_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_tiara"],
    threads: 50
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),6,100000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),6,100000))
    shell:
        """
        set +eu
        source {params.conda_exe}
        conda activate {params.conda_env}
        tiara -i {input} -o {output} -t {threads} -m 200 
        conda deactivate
        """

rule tiara_nonprok:
    input:
        data=lambda wildcards: get_input_fa(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
        tiara="tiara/tiara_{datatype}_{sample}_{subset_gb}gb.tsv",
    output:
        id="tiara/tiara_{datatype}_{sample}_{subset_gb}gb_id.tsv",
        df="tiara/nanoq_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_nanoq"],
    threads: 10
    resources:
        runtime="6d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,50000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        cut -f1,2 {input.tiara} | grep -v -e "prokarya" -e "bacteria" -e "archaea" - | cut -f1 | cut -f1 -d' ' | sort > {output.id}
        seqkit grep -f {output.id} {input.data} | nanoq -s -f --report {output.df}
        conda deactivate
        """

rule tiara_sum:
    input: "tiara/nanoq_{datatype}_{sample}_{subset_gb}gb.tsv",   
    output: "tiara/stats_{datatype}_{sample}_{subset_gb}gb.tsv"
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        name="{wildcards.sample}_{wildcards.subset_gb}"
        nonprok_n=$(awk '{{print $1}}' {input})
        nonprok_bp=$(awk '{{print $2}}' {input})
        nonprok_n50=$(awk '{{print $3}}' {input})
        echo -e "${{name}}\t${{nonprok_n}}\t${{nonprok_bp}}\t${{nonprok_n50}}" > {output}
        """

rule tiara_agr:
    input: expand("tiara/stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype=config["datatypes"],sample=config["sample"],subset_gb=config["subset_gb"])
    output:
        reads="results/tiara_reads.tsv",
        contigs="results/tiara_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        samples=({input})
        
        reads=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "reads")
        echo -e "sample\ttiara_read_n\ttiara_read_bp\ttiara_read_n50" > {output.reads}
        cat $reads >> {output.reads}
        
        contigs=$(echo "${{samples[@]}}" | tr ' ' '\n' | grep "contigs")
        echo -e "sample\ttiara_contig_n\ttiara_contig_bp\ttiara_contig_n50" > {output.contigs}
        cat $contigs >> {output.contigs}
        """

rule jellyfish_count:
    input: lambda wildcards: get_input_fa(wildcards.sample,wildcards.subset_gb,wildcards.datatype),
    output: "jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb.jf",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_jellyfish"],
    threads: 60
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),8,100000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),8,100000))
    shell:
        """
        set +eu
        source {params.conda_exe}
        conda activate {params.conda_env}
        jellyfish count -m 24 -s 5000M -t {threads} -C -o {output} {input}
        conda deactivate
        """

rule jellyfish_hist:
    input: "jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb.jf",
    output: "jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb_hist.txt",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_jellyfish"],
    threads: 30
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,50000))
    shell:
        """
        set +eu
        source {params.conda_exe}
        conda activate {params.conda_env}
        jellyfish histo -t {threads} -o {output} {input}
        sed -i 's/$/ {wildcards.sample}_{wildcards.subset_gb}/' {output}
        conda deactivate
        """

rule jellyfish_stats:
    input: "jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb.jf",
    output:
        df1="jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb_stats.tsv",
        df2="jellyfish/stats_{datatype}_{sample}_{subset_gb}gb.tsv",
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_jellyfish"],
    threads: 10
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,30000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,30000)),
    shell:
        """
        set +eu
        source {params.conda_exe}
        conda activate {params.conda_env}
        jellyfish stats -o {output.df1} {input}
        awk -F: '{{gsub(/^[ \t]+|[ \t]+$/, "", $2); printf "%s\t", $2}}' {output.df1} | sed 's/\t$//' | sed 's/^/{wildcards.sample}_{wildcards.subset_gb}\t/' > {output.df2}
        conda deactivate
        """

rule jellyfish_agr:
    input:
        hist=expand("jellyfish/jellyfish_{datatype}_{sample}_{subset_gb}gb_hist.txt",datatype="reads",sample=config["sample"],subset_gb=config["subset_gb"]),
        sum=expand("jellyfish/stats_{datatype}_{sample}_{subset_gb}gb.tsv",datatype="reads",sample=config["sample"],subset_gb=config["subset_gb"]),
    output:
        hist="results/jellyfish_hist.txt",
        sum="results/jellyfish_sum.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        echo -e "kmer_freq kmer_count sample" > {output.hist}
        cat {input.hist} >> {output.hist}
        echo -e "sample\tkmer_uniq\tkmer_dist\tkmer_total\tkmer_max" > {output.sum}
        paste -d '\n' {input.sum} >> {output.sum}
        """

rule longshot_prep:
    input: "contigs/{sample}_{subset_gb}.fasta",   
    output: "{sample}_{subset_gb}/tmp/binning/mapping_tmp/1_cov.bam.bai",  
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_longshot"],
    threads: 5
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),1,20000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),1,20000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        samtools index {wildcards.sample}_{wildcards.subset_gb}/tmp/binning/mapping_tmp/1_cov.bam
        conda deactivate
        """

rule longshot:
    input: "{sample}_{subset_gb}/tmp/binning/mapping_tmp/1_cov.bam.bai",   
    output: "longshot/longshot_{sample}_{subset_gb}.vcf",  
    params:
        conda_exe=config["conda_exe"],
        conda_env=config["env_longshot"],
    threads: 60
    resources:
        runtime="12d",
        mem_mb=lambda wildcards, attempt: get_mem(int(wildcards.subset_gb),int(attempt),8,50000),
        slurm_partition=lambda wildcards, attempt: get_partition(get_mem(int(wildcards.subset_gb),int(attempt),8,50000))
    shell:
        """
        source {params.conda_exe}
        conda activate {params.conda_env}
        longshot --bam {wildcards.sample}_{wildcards.subset_gb}/tmp/binning/mapping_tmp/1_cov.bam --ref contigs/{wildcards.sample}_{wildcards.subset_gb}.fasta --out {output} -n --min_cov 6 --min_mapq 18 --min_alt_frac 0.1 --min_alt_count 3
        conda deactivate
        """

rule longshot_sum:
    input: "longshot/longshot_{sample}_{subset_gb}.vcf"
    output:
        asm="longshot/con_{sample}_{subset_gb}.tsv",
        var_sum="longshot/sum_{sample}_{subset_gb}.tsv",
        var_asm="longshot/asm_{sample}_{subset_gb}.tsv",
        var_bin="longshot/bin_{sample}_{subset_gb}.tsv",
    threads: 5
    resources:
        runtime="1d",
        mem_mb=30*1000,
        slurm_partition="general",
    shell:
        """
        awk -F "\t" '{{ if ($7 == "PASS") {{print $1}} }}' {input} | uniq -c - | awk '{{$1=$1;print}}' | awk -F" " '{{print $2,$1}}' OFS='\t' | sort -k1,1 > {output.var_sum}
        sed 1d {wildcards.sample}_{wildcards.subset_gb}/tmp/flye/assembly_info.txt | cut -f1,2,3,4 | sort -k1,1 > {output.asm}
        join -a1 -a2 -e 0 -o 'auto' -t $'\t' <(cat {output.asm}) <(cat {output.var_sum}) | sort -k1,1 > {output.var_asm}
        join -a1 -a2 -e "Unbinned" -o 'auto' -t $'\t' <(cat {output.var_asm}) <(sort -k1,1 {wildcards.sample}_{wildcards.subset_gb}/tmp/binning/contig_bin.tsv) | sed 's/^/{wildcards.sample}_{wildcards.subset_gb}\t/' > {output.var_bin}
        """

rule longshot_agr:
    input: expand("longshot/bin_{sample}_{subset_gb}.tsv",sample=config["sample"],subset_gb=config["subset_gb"])
    output: "results/longshot_contigs.tsv",
    threads: 1
    resources:
        runtime="1d",
        mem_mb=5*1000,
        slurm_partition="general",
    shell:
        """
        echo -e "sample\tcontig_id\tcontig_len\tcontig_cov\tcontig_circ\tlongshot_vars\tbin" > {output}
        cat {input} >> {output}
        """
