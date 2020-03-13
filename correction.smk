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

ALN_THREADS = 16

configfile: "cor.yaml"
SMS = list(config.keys())
refs = {}
fofns = {}
algs = {}
for SM in SMS:
	refs[SM] = config[SM]["ref"]
	# paried end reads must be given in a paired fofn, where paired fastqs are on the same line sperated by a comma	
	fofns[SM] = config[SM]["fofn"]
	alg = "racon"
	if("alg" in config[SM]):
		alg = config[SM]["alg"]
	if(alg not in algs):
		algs[alg] = []
	algs[alg].append(SM)
		
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

# TODO update pilon and qrrow/quiver to work with 
ALG="NONE"

wildcard_constraints:
	ID="\d+",
	READID="\d+",
	SM="|".join(SMS),

localrules: all,

rule all:
	input:
		final = expand("{SM}_racon.fasta", SM=algs["racon"]),


###############################
### Racon specifc part      ###
###############################

def chunk_by_size(df, num):
	out = {}
	counter = 0
	tigs = df["contig"]
	for tig in tigs:
		if(counter not in out):
			out[counter] = []	
		out[counter].append(tig)
			
		counter += 1
		if(counter >= num):
			counter = 0
	return(out)




checkpoint split_ref:
	input:
		ref=get_ref,
		fai=get_fai,
	output: 
		split = directory("temp/{SM}_split_ref/"),
		#splitref=temp(expand("temp/{{SM}}_split_ref/{ID}.fasta", ID=IDs)),
		#splitrgn=temp(expand("temp/{{SM}}_split_ref/{ID}.rgn", ID=IDs)),
	benchmark:
		"logs/split_ref_{SM}.b"
	resources:
		mem=6
	threads:1
	run:
		tigs = pd.read_table(input["fai"], header=None, names=["contig", "len", "x", "y", "z"])
		tigs.sort_values(by=["len"], inplace=True)
		chunks = min(100, len(tigs["contig"]))
		batches = chunk_by_size(tigs, chunks)

		d_ref = SeqIO.to_dict(SeqIO.parse(input["ref"], "fasta"))
		for ID in range(chunks):
			rgn = open("temp/{}_split_ref/{}.rgn".format(wildcards["SM"], ID)  , "w+")
			tigs = batches[ID]
			recs = []
			for tig in tigs:
				recs.append(d_ref[tig])
				rgn.write( tig + "\n")
			
			SeqIO.write(recs, "temp/{}_split_ref/{}.fasta".format(wildcards["SM"], ID), "fasta")
			rgn.close()

#
# CREATE ALIGNMENTS
#
rule minimap:
	input:
		ref = get_ref,
		reads = get_reads,
	output:
		bam = temp("temp/{SM}_{READID}_alignments.bam"),
	benchmark:
		"logs/minimap_{SM}_{READID}.b"
	resources:
		mem=8,
		smem=8,
	threads:ALN_THREADS
	shell:"""
minimap2  -a -x map-pb -s 2500 -t {threads} --secondary=no {input.ref}  {input.reads} \
| samtools view -u -F 1796 - | \
samtools sort -m {resources.smem}G -@ {threads} - > {output.bam} 
"""


def get_all_bams(wildcards):
	SM = str(wildcards.SM)
	reads = [ line.strip() for line in open(config[SM]["fofn"]).readlines() ]
	READIDS = list(range(len(reads)))
	return( expand("temp/{SM}_{READID}_alignments.bam", SM=SM, READID=READIDS) )
	
rule merge_align:
	input:
		ref = get_ref,
		bams = get_all_bams,
	output:
		bam = temp("temp/{SM}_alignments.bam"),
	benchmark:
		"logs/merge_align_{SM}.b"
	resources:
		mem=2
	threads:ALN_THREADS
	run:
		if(len(input["bams"]) == 1):
			shell("ln {input.bams} {output.bam}")
		else: 
			shell("samtools merge -@ {threads} {output.bam} {input.bams}")

rule index_merge_align:
	input:
		bam = "temp/{SM}_alignments.bam",
	output:
		bai = temp("temp/{SM}_alignments.bam.bai"),
	benchmark:
		"logs/index_merge_align_{SM}.b"
	resources:
		mem=8
	threads:1
	shell:"""
samtools index {input.bam}
"""


#
# RUN CORRECTION 
#
rule split_bam_to_sam:
	input:
		bam = "temp/{SM}_alignments.bam",
		bai = "temp/{SM}_alignments.bam.bai",
		splitrgn="temp/{SM}_split_ref/{ID}.rgn",
	output:
		splitsam = temp("temp/{SM}_alignments/{ID}.sam"),
	benchmark:
		"logs/splot_bam_to_sam_{SM}_{ID}.b"
	resources:
		mem=4
	threads:1
	shell:"""
samtools view -h {input.bam} $( cat {input.splitrgn} ) > {output.splitsam}
"""

rule split_fastq:
	input:
		splitsam = "temp/{SM}_alignments/{ID}.sam",
	output:
		splitfastq = temp("temp/{SM}_split_fasta/{ID}.fastq"),
	benchmark:
		"logs/split_fastq_{SM}_{ID}.b"
	resources:
		mem=4
	threads:1
	shell:"""
samtools fastq {input.splitsam} > {output.splitfastq}
"""	

rule run_racon:
	input:
		splitsam = "temp/{SM}_alignments/{ID}.sam",
		splitfastq = "temp/{SM}_split_fasta/{ID}.fastq",
		splitref= "temp/{SM}_split_ref/{ID}.fasta",
	output:
		splitracon= temp("temp/{SM}_split_cor/{ID}.fasta"),
	benchmark:
		"logs/run_racon_{SM}_{ID}.b"
	resources:
		mem=8
	threads:4
	shell:"""
racon {input.splitfastq} {input.splitsam} {input.splitref}  -u -t {threads} > {output.splitracon}
"""



def get_split_cor(wildcards):
	checkpoint_output = checkpoints.split_ref.get(**wildcards).output[0]
	IDS = glob_wildcards(os.path.join(checkpoint_output, "{ID}.fasta")).ID
	SM = wildcards.SM
	rtn = expand("temp/{SM}_split_cor/{ID}.fasta", SM=SM, ID=IDS) 
	return(rtn)

rule combine_cor:
	input:
		get_split_cor
		#splitcor= expand("temp/{{SM}}_split_cor/{ID}.fasta", ID=IDs),
		#splitref= expand("temp/{{SM}}_split_ref/{ID}.fasta", ID=IDs),
		#bam = "temp/{SM}_alignments.bam",
		#bai = "temp/{SM}_alignments.bam.bai",
	output:
		final = "{SM}_racon.fasta",
		fai = "{SM}_racon.fasta.fai",
	benchmark:
		"logs/combine_cor_{SM}.b"
	resources:
		mem=8
	threads:1
	shell:"""
cat {input} | seqtk seq -l 80 > {output.final}
samtools faidx {output.final}
"""









if(ALG == "pilon"):
	rule bwa_index:
		input:
			ref = REF,
		output:
			sa  = "{SM}_bwa_index/ref.sa",
			amb = "{SM}_bwa_index/ref.amb",
			ann = "{SM}_bwa_index/ref.ann",
			pac = "{SM}_bwa_index/ref.pac",
			bwt = "{SM}_bwa_index/ref.bwt",
		benchmark:
			"logs/bwa_index_{SM}.b"
		resources:
			mem=16
		threads:1
		shell:"""
mkdir -p {wildcards.SM}_bwa_index
bwa index -b 500000000 {input.ref} -p {wildcards.SM}_bwa_index/ref
"""

	rule bwa:
		input:
			sa  = "{SM}_bwa_index/ref.sa",
			amb = "{SM}_bwa_index/ref.amb",
			ann = "{SM}_bwa_index/ref.ann",
			pac = "{SM}_bwa_index/ref.pac",
			bwt = "{SM}_bwa_index/ref.bwt",
			reads = readByID,
		output:
			bam = temp("temp/{SM}_{READID}_alignments.bam"),
		benchmark:
			"logs/bwa_{SM}_{READID}.b"
		resources:
			mem=8,
			sortmem=4,
		threads:8
		shell:"""
samtools --version 
bwa mem -t {threads} {wildcards.SM}_bwa_index/ref {input.reads} \
| samtools view -u -F 1796 - | \
samtools sort -T tmp_{wildcards.SM}_{wildcards.READID} -m {resources.sortmem}G -@ {threads} - > {output.bam} 
"""


	rule split_bam:
		input:
			bam = "{SM}_alignments.bam",
			bai = "{SM}_alignments.bam.bai",
		output:
			splitbam = temp("temp/{SM}_alignments/{ID}.bam"),
			splitbai = temp("temp/{SM}_alignments/{ID}.bam.bai"),
		benchmark:
			"logs/split_bam_{SM}_{ID}.b"
		resources:
			mem=4
		threads:1
		run:
			ID = int(str(wildcards["ID"]))
			tigs = sorted(batches[ID])
			shell("samtools view -b {input.bam} {tigs} > {output.bam}")
			shell("samtools index {output.bam} ") 

	rule run_pilon:
		input:
			splitref=  "temp/{SM}_split_ref/{ID}.fasta",
			splitbam = "temp/{SM}_alignments/{ID}.bam",
			splitbai = "temp/{SM}_alignments/{ID}.bam.bai",
		output:
			splitracon= temp("temp/{SM}_split_cor/{ID}.fasta"),
		log:
			o="logs/{rule.name}_{SM}_{ID}.o",
			e="logs/{rule.name}_{SM}_{ID}.e",
		benchmark:
			"logs/run_pilon_{SM}_{ID}.b"
		resources:
			mem=16
		threads:1
		shell:"""
java -Xmx16G -jar {PILON_JAR} --genome {input.splitref} --bam {input.splitbam} --output temp/{wildcards.SM}_split_cor/{wildcards.ID}
"""


	rule combine_cor:
		input:
			splitcor= expand("temp/{{SM}}_split_cor/{ID}.fasta", ID=IDs),
			splitref= expand("temp/{{SM}}_split_ref/{ID}.fasta", ID=IDs),
			bam = "{SM}_alignments.bam",
			bai = "{SM}_alignments.bam.bai",
		output:
			final = "{SM}_pilon.fasta",
			fai = "{SM}_pilon.fasta.fai",
		benchmark:
			"logs/combine_cor_{SM}.b"
		resources:
			mem=8
		threads:1
		shell:"""
	cat {input.splitcor} | seqtk seq -l 80 > {output.final}
	samtools faidx {output.final}
	"""










######################################
### Arrow/quvier specifc part      ###
######################################
if(ALG in ["arrow", "quiver"] ):

	#### read in alignment workflow ####
	#subworkflow pbalign:
	#    workdir:
	#        SNAKEMAKE_DIR + "pb_align.smk"
	#    snakefile:
	#        "../path/to/otherworkflow/Snakefile"
	#    configfile:
	#        "path/to/custom_configfile.yaml"

	chunks = 200
	if("chunks" in config):
		chunks = config["chunks"]
	chunks = min(chunks, len(tigs["contig"]))

	IDs = list(range(chunks))
	XML = config["xml"]
	xmlpre =  "temp/{SM}_arrow_in/" + ".".join( os.path.basename(XML).split(".")[:-2] )
	xmlsuf = ".".join( os.path.basename(XML).split(".")[-2:] )
	print(xmlpre, xmlsuf)	
	assert xmlsuf == "alignmentset.xml", xmlsuf

	rule xml_split:
		input:
			xml = XML,
		output:
			xmls = expand( xmlpre + ".chunk{ID}." + xmlsuf, ID=IDs),
		resources:
			mem=32
		benchmark:
			"logs/xml_split_{SM}_{ID}.b"
		threads:1
		shell:"""
mkdir -p temp/{wildcards.SM}_arrow_in
dataset split --contigs --chunks {chunks} --outdir temp/{wildcards.SM}_arrow_in {input.xml}
"""

	rule xml_ref:
		input:
			ref = REF,
		output:
			xml = "{SM}_ref.xml",
		resources:
			mem=8
		benchmark:
			"logs/xml_ref_{SM}_{ID}.b"
		threads:1
		shell:"""
dataset create --type ReferenceSet {output.xml} {input.ref}
"""


	rule arrow_by_xml:
		input:
			ref = "{SM}_ref.xml",
			xml = xmlpre + ".chunk{ID}." + xmlsuf,  
		output:
			fasta = "temp/{SM}_arrow_out/{ID}.fasta",
			fastq = "temp/{SM}_arrow_out/{ID}.fastq",
			vcf = "temp/{SM}_arrow_out/{ID}.vcf",
		benchmark:
			"logs/arrow_by_xml_{SM}_{ID}.b"
		resources:
			mem=lambda wildcards, attempt: 8 + 4*attempt, # tested and 4 is not enough for CCS input 
		threads:6
		shell:"""
variantCaller --log-level INFO -j {threads} \
	--algorithm {ALG} --noEvidenceConsensusCall lowercasereference \
	--alignmentSetRefWindows -r {input.ref} \
	-o {output.fasta} -o {output.fastq} -o {output.vcf}\
	 {input.xml}
"""

	rule merge_out:	
		input:
			fasta = expand("temp/{SM}_arrow_out/{ID}.fasta", ID=IDs),
		output:
			final = "{SM}_corrected.fasta",
		benchmark:
			"logs/merge_out_{SM}_{ID}.b"
		resources:
			mem=16
		threads:1
		shell:"""
cat {input.fasta}  | seqtk seq -l 80 > {output.final} && samtools faidx {output.final}
"""












