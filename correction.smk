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




configfile: "cor.yaml"
REF = config["ref"]
FAI = REF + ".fai"


PRE=""
if("prefix" in config):
	PRE=config["prefix"] 

ALG="arrow"
if("alg" in config):
	ALG=config["alg"].lower()

localrules: all,

wildcard_constraints:
	ID="\d+",
	READID="\d+",

#
#
#
rule all:
	input:
		final = PRE+"corrected.fasta",


#
# try to make input index 
#
if(not os.path.exists(FAI)):
	shell( "samtools faidx {}".format(REF) )	
tigs = pd.read_table(FAI, header=None, names=["contig", "len", "x", "y", "z"]); tigs.sort_values(by=["len"], inplace=True)

rule index:
	input: REF
	output: FAI
	log:
		"logs/{rule}_{ID}.o",
		"logs/{rule}_{ID}.e",
	benchmark:
		"logs/{rule}_{ID}.b"
	resources:
		mem=6
	threads:1
	shell: "samtools faidx {input}"


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
	xmlpre = PRE + "arrow_in/" + ".".join( os.path.basename(XML).split(".")[:-2] )
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
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		threads:1
		shell:"""
mkdir -p arrow_in
dataset split --contigs --chunks {chunks} --outdir arrow_in {input.xml}
"""

	rule xml_ref:
		input:
			ref = REF,
		output:
			xml = PRE+"ref.xml",
		resources:
			mem=8
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		threads:1
		shell:"""
dataset create --type ReferenceSet {output.xml} {input.ref}
"""


	rule arrow_by_xml:
		input:
			ref = PRE+"ref.xml",
			xml = xmlpre + ".chunk{ID}." + xmlsuf,  
		output:
			fasta = PRE+"arrow_out/{ID}.fasta",
			fastq = PRE+"arrow_out/{ID}.fastq",
			vcf = PRE+"arrow_out/{ID}.vcf",
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=lambda wildcards, attempt: 8 + 4*attempt, # tested and 4 is not enough for CCS input 
		threads:6
		shell:"""
{{ variantCaller --log-level INFO -j {threads} \
	--algorithm {ALG} --noEvidenceConsensusCall lowercasereference \
	--alignmentSetRefWindows -r {input.ref} \
	-o {output.fasta} -o {output.fastq} -o {output.vcf}\
	 {input.xml} ; }} &>> {log}
"""

	rule merge_out:	
		input:
			fasta = expand(PRE+"arrow_out/{ID}.fasta", ID=IDs),
		output:
			final = PRE+"corrected.fasta",
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=16
		threads:1
		shell:"""
cat {input.fasta}  | seqtk seq -l 80 > {output.final} && samtools faidx {output.final}
"""











###############################
### Racon specifc part      ###
###############################

if(ALG in ["racon", "pilon"]):
	
	# paried end reads must be given in a apired fofn, where paired fastqs are on the same line sperated by a comma	
	reads = [ line.strip() for line in open(config["fofn"]).readlines() ]
	fais = [read + ".fai" for read in reads]
	readIDs = list(range(len(reads)))
		
	def readByID(wildcards):
		idx = int(str(wildcards.READID))
		# split them in case they are paired
		rtn = reads[idx]
		rtn = [ read.strip() for read in rtn.split(",") ]
		return(rtn)	
	
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

	def print_batch(out, df):
		for key in out:
			tmp = df[ df["contig"].isin( out[key] ) ]
			print(tmp['len'].sum(), len(out[key]), file=sys.stderr)
		print(len(out), file=sys.stderr)
		

	chunks = 100
	if("chunks" in config):
		chunks = config["chunks"]
	chunks = min(chunks, len(tigs["contig"]))
	IDs = list(range(chunks))

	batches = chunk_by_size(tigs, chunks)
	#print_batch(batches, tigs)

	rule split_ref:
		input:
			ref=REF,
			fai=FAI,
		output: 
			splitref=temp(expand(PRE+"split_ref/{ID}.fasta", ID=IDs)),
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=6
		threads:1
		run:
			d_ref = SeqIO.to_dict(SeqIO.parse(input["ref"], "fasta"))
			for ID in IDs:
				tigs = batches[ID]
				recs = []
				for tig in tigs:
					recs.append(d_ref[tig])
				SeqIO.write(recs, PRE+"split_ref/{}.fasta".format(ID), "fasta")

	#
	# CREATE ALIGNMENTS
	#
	if(ALG == "racon"):
		rule minimap:
			input:
				ref = REF,
				reads = readByID,
			output:
				bam = temp(PRE+"{READID}_alignments.bam"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=4,
				smem=2,
			threads:6
			shell:"""
minimap2  -a -x map-pb -m 5000 -t {threads} --secondary=no {input.ref}  {input.reads} \
	| samtools view -u -F 1796 - | \
	samtools sort -m {resources.smem}G -@ {threads} - > {output.bam} 
"""

	elif(ALG == "pilon"):
		rule bwa_index:
			input:
				ref = REF,
			output:
				sa  = PRE+"bwa_index/ref.sa",
				amb = PRE+"bwa_index/ref.amb",
				ann = PRE+"bwa_index/ref.ann",
				pac = PRE+"bwa_index/ref.pac",
				bwt = PRE+"bwa_index/ref.bwt",
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=16
			threads:1
			shell:"""
mkdir -p {PRE}bwa_index
bwa index -b 500000000 {input.ref} -p {PRE}bwa_index/ref
"""

		rule bwa:
			input:
				sa  = PRE+"bwa_index/ref.sa",
				amb = PRE+"bwa_index/ref.amb",
				ann = PRE+"bwa_index/ref.ann",
				pac = PRE+"bwa_index/ref.pac",
				bwt = PRE+"bwa_index/ref.bwt",
				reads = readByID,
			output:
				bam = temp(PRE+"{READID}_alignments.bam"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=8,
				sortmem=4,
			threads:8
			shell:"""
samtools --version 
bwa mem -t {threads} {PRE}bwa_index/ref {input.reads} \
	| samtools view -u -F 1796 - | \
	samtools sort -T tmp{wildcards.READID} -m {resources.sortmem}G -@ {threads} - > {output.bam} 
"""

	
	rule merge_align:
		input:
			ref = REF,
			bams = expand(PRE+"{READID}_alignments.bam", READID=readIDs),
			#bais = expand(PRE+"{READID}_alignments.bam.bai", READID=readIDs),
		output:
			bam = PRE+"alignments.bam",
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=2
		threads:8
		shell:"""
samtools merge -@ {threads} {output.bam} {input.bams}
"""
	
	rule index_merge_align:
		input:
			bam = PRE+"alignments.bam",
		output:
			bai = PRE+"alignments.bam.bai",
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=8
		threads:1
		shell:"""
samtools index {input.bam}
"""


	#
	# RUN CORRECTION 
	#
	if(ALG == "racon"):
		rule split_bam:
			input:
				bam = PRE+"alignments.bam",
				bai = PRE+"alignments.bam.bai",
			output:
				splitsam = temp(PRE+"alignments/{ID}.sam"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=4
			threads:1
			run:
				ID = int(str(wildcards["ID"]))
				tigs = sorted(batches[ID])
				out = PRE + "alignments/{}.sam".format(ID)
				shell("samtools view -h {input.bam} {tigs} > " + out)

		
		rule split_fastq:
			input:
				splitsam = PRE+"alignments/{ID}.sam",
			output:
				splitfastq = temp(PRE+"split_fasta/{ID}.fastq"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=4
			threads:1
			shell:"""
	samtools fastq {input.splitsam} > {output.splitfastq}
	"""	

		rule run_racon:
			input:
				splitfastq = PRE+"split_fasta/{ID}.fastq",
				splitref=  PRE+"split_ref/{ID}.fasta",
				splitsam = PRE+"alignments/{ID}.sam",
				#allsplitfastq = expand(PRE+"split_fasta/{ID}.fastq",ID=IDs), # forces all the fastqs to be build first 
				#allsplitref = expand(PRE+"split_ref/{ID}.fasta",ID=IDs), # forces all the refs to be build first 
			output:
				splitracon= temp(PRE+"split_cor/{ID}.fasta"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=15
			threads:1
			shell:"""
	racon {input.splitfastq} {input.splitsam} {input.splitref}  -u -t {threads} > {output.splitracon}
	"""


	elif(ALG == "pilon"):
		rule split_bam:
			input:
				bam = PRE+"alignments.bam",
				bai = PRE+"alignments.bam.bai",
			output:
				splitbam = temp(PRE+"alignments/{ID}.bam"),
				splitbai = temp(PRE+"alignments/{ID}.bam.bai"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=4
			threads:1
			run:
				ID = int(str(wildcards["ID"]))
				tigs = sorted(batches[ID])
				out = PRE + "alignments/{}.bam".format(ID)
				shell("samtools view -b {input.bam} {tigs} > " + out)
				shell("samtools index {} ".format( out) ) 

		rule run_pilon:
			input:
				splitref=  PRE+"split_ref/{ID}.fasta",
				splitbam = PRE+"alignments/{ID}.bam",
				splitbai = PRE+"alignments/{ID}.bam.bai",
				#allsplitfastq = expand(PRE+"split_fasta/{ID}.fastq",ID=IDs), # forces all the fastqs to be build first 
				#allsplitref = expand(PRE+"split_ref/{ID}.fasta",ID=IDs), # forces all the refs to be build first 
			output:
				splitracon= temp(PRE+"split_cor/{ID}.fasta"),
			log:
				"logs/{rule}_{ID}.o",
				"logs/{rule}_{ID}.e",
			benchmark:
				"logs/{rule}_{ID}.b"
			resources:
				mem=16
			threads:1
			shell:"""
java -Xmx16G -jar {PILON_JAR} --genome {input.splitref} --bam {input.splitbam} --output {PRE}split_cor/{wildcards.ID}
"""





	rule combine_cor:
		input:
			splitcor= expand(PRE+"split_cor/{ID}.fasta", ID=IDs),
			splitref= expand(PRE+"split_ref/{ID}.fasta", ID=IDs),
			bams = expand(PRE+"{READID}_alignments.bam", READID=readIDs),
			bam = PRE+"alignments.bam",
			bai = PRE+"alignments.bam.bai",
		output:
			final = PRE+"corrected.fasta",
			fai = PRE+"corrected.fasta.fai",
		log:
			"logs/{rule}_{ID}.o",
			"logs/{rule}_{ID}.e",
		benchmark:
			"logs/{rule}_{ID}.b"
		resources:
			mem=8
		threads:1
		shell:"""
cat {input.splitcor} | seqtk seq -l 80 > {output.final}
samtools faidx {output.final}
"""






