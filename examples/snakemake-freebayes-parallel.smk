## A simple example snakemake .smk file for parallelising freebayes
## Uses a simple R script to split the genome into regions of equal size based on the .fai index
## As snakemake automatically moves each cpu core to the next genome chunk, this works out faster
## than the freebayes-parallel wrapper and allows pooled sample calling. 
## This .smk file assumes we have a list of the bam files called bam.list
## This .smk file splits the genome by chromosome, which of course, is not necessary.
## One will want to edit the paths (for example, the path to bam files)
import numpy as np 

# these parameters should usually be stored in the snakemake configuration file (config.yaml) and accessed e.g. config['ref']
samples = ['SampleA', 'SampleB', 'SampleC']
reference = "path/to/reference"
chroms = [1,2,3]
nchunks = 9

bamlist = "path/to/bam.list"
chunks = np.arange(1, nchunks)

rule GenomeIndex:
    input:
        ref = reference
    output:
        idx = reference + ".fai"
    log: 
        "logs/GenomeIndex.log"
    wrapper: 
        "v0.69.0/bio/samtools/faidx"


rule GenerateFreebayesRegions:
    input:
        ref_idx = reference,
        index = reference + ".fai",
        bams = expand("resources/alignments/{sample}.bam", sample=samples)
    output:
        regions = expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=chroms, i = chunks)
    log:
        "logs/GenerateFreebayesRegions.log"
    params:
        chroms = chroms
        chunks = chunks
    conda:
        "../envs/freebayes-env.yaml"
    script:
        "../scripts/GenerateFreebayesRegions.R" # This is located in the scripts/ directory of freebayes


rule VariantCallingFreebayes:
	input:
		bams = expand("resources/alignments/{sample}.bam", sample=samples),
		index = expand("resources/alignments/{sample}.bam.bai", sample=samples),
		ref = reference,
		samples = bamlist,
		regions = "resources/regions/genome.{chrom}.region.{i}.bed"
	output:
		temp("results/variants/vcfs/{chrom}/variants.{i}.vcf")
	log:
		"logs/VariantCallingFreebayes/{chrom}.{i}.log"
	conda:
        "../envs/freebayes-env.yaml"
	threads:1
	shell:	"freebayes -f {input.ref} -t {input.regions} -L {input.samples} > {output} 2> {log}"


rule ConcatVCFs:
    input:
        calls = expand("results/variants/vcfs/{{chrom}}/variants.{i}.vcf", i=chunks)
    output:
        "results/variants/vcfs/variants.{chrom}.vcf"
    log:
        "logs/ConcatVCFs/{chrom}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads:4
    shell:  
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"