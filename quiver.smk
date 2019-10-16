import os
import tempfile
import numpy as np
import pandas as pd
import json
from pprint import pprint
import itertools 
import pysam 
from Bio import SeqIO
import random


#
# setup the env for each execution
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile) + "/"
shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

PILON_JAR = "/net/eichler/vol26/7200/software/modules-sw/pilon/1.22/Linux/RHEL6/x86_64/pilon-1.22.jar"

#
# Get tmp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

MAX_THREADS=4

configfile: "cor.yaml"
SMS = list(config.keys())
refs = {}
fofns = {}
fais = {}
ids = {}
algs = {}
rgns = {}
for SM in SMS:
	refs[SM] = config[SM]["ref"]
	# paried end reads must be given in a paired fofn, where paired fastqs are on the same line sperated by a comma	
	fofns[SM] = config[SM]["fofn"]
	if("alg" in config[SM]):
		alg = config[SM]["alg"]
	algs[SM] = alg	

	# read in fai make regions
	fais[SM] = config[SM]["ref"] + ".fai"
	assert os.path.exists(fais[SM])
	rgns[SM] = []
	for line in open(fais[SM]).readlines():
		tokens=line.strip().split()
		rgns[SM].append( "'{}:{}-{}'".format(tokens[0], 0, tokens[1] ) )
	ids[SM] = list(range(len(rgns[SM])))
	
		
def get_ref(wildcards):
	SM = str(wildcards.SM)
	return(refs[SM])

def get_fai(wildcards):
	SM = str(wildcards.SM)
	return(refs[SM] + ".fai")

def get_reads(wildcards):
	SM = str(wildcards.SM)
	idx = int(str(wildcards.READID))
	# split them in case they are paired
	reads = [ line.strip() for line in open(config[SM]["fofn"]).readlines() ]
	rtn = reads[idx]
	rtn = [ read.strip() for read in rtn.split(",") ]
	return(rtn)	

def get_fofn(wildcards):
	SM = str(wildcards.SM)
	return(fofns[SM])

def get_rgn(wildcards):
	SM = str(wildcards.SM)
	ID = int(str(wildcards.ID))
	return( rgns[SM][ID] )	
	
def get_alg(wildcards):
	SM = str(wildcards.SM)
	return( algs[SM] )	

wildcard_constraints:
	ID="\d+",
	READID="\d+",
	SM="|".join(SMS),

localrules: all,

def get_cor_fasta(wildcards):
	rtn = []
	for sm in algs:
		alg = algs[sm]
		if(alg == "racon"):
			rtn.append( "{SM}_racon.fasta".format(SM=sm) )
		elif(alg == "quiver" or alg=="arrow"):
			rtn.append( "{SM}_pb_cor.fasta".format(SM=sm) )
	return(rtn)

rule all:
	input:
		final = get_cor_fasta, 


######################################
### Arrow/quvier specifc part      ###
######################################



rule make_contig_bam:
	input:
		ref = get_ref,
		bam = get_fofn,  
	output:
		bam = temp("temp/{SM}_arrow_out/{ID}.bam"),
		bai = temp("temp/{SM}_arrow_out/{ID}.bam.bai"),
		pbi = temp("temp/{SM}_arrow_out/{ID}.bam.pbi"),
	benchmark:
		"logs/make_contig_bam_{SM}_{ID}.b"
	params:
		rgn = get_rgn,
	resources:
		mem=2
	threads: MAX_THREADS
	shell:"""
samtools merge -R {params.rgn} -@ {threads} {output.bam} -b {input.bam} -O BAM && \
	samtools index {output.bam} && \
	pbindex {output.bam}
"""

rule arrow:
	input:
		ref = get_ref,
		bam = "temp/{SM}_arrow_out/{ID}.bam",
		bai = "temp/{SM}_arrow_out/{ID}.bam.bai",
		pbi = "temp/{SM}_arrow_out/{ID}.bam.pbi",
	output:
		fasta = "temp/{SM}_arrow_out/{ID}.fasta",
		fastq = "temp/{SM}_arrow_out/{ID}.fastq",
		vcf = "temp/{SM}_arrow_out/{ID}.vcf",
	benchmark:
		"logs/arrow_{SM}_{ID}.b"
	params:
		rgn = get_rgn,
		alg = get_alg,
	resources:
		mem=lambda wildcards, attempt: 4*attempt, # tested and 4 is not enough for CCS input 
	threads: MAX_THREADS
	shell:"""
variantCaller --log-level INFO -j {threads} \
	--algorithm {params.alg} --noEvidenceConsensusCall lowercasereference \
	--referenceWindow {params.rgn} \
	-r {input.ref} \
	-o {output.fasta} -o {output.fastq} -o {output.vcf} \
	 {input.bam} 
"""

def get_out_fasta(wildcards):
	SM = str(wildcards.SM)
	IDS = ids[SM]
	fastas = []
	for ID in IDS:
		fastas.append( "temp/{SM}_arrow_out/{ID}.fasta".format(SM=SM, ID=ID) )
	return( fastas )	

rule merge_out:	
	input:
		fasta = get_out_fasta,
	output:
		final = "{SM}_pb_cor.fasta",
	benchmark:
		"logs/merge_out_{SM}.b"
	resources:
		mem=16
	threads:1
	shell:"""
cat {input.fasta}  | seqtk seq -l 80 > {output.final} && samtools faidx {output.final}
"""












