import os
import numpy as np
import pandas as pd
import json
import random
import tempfile

#
# setup the env for each execution
# this allows you to specify a file that contains all the modules you want to have loaded
# for your pipeline.
# the file should be called modules.cfg and should exist in the same locaiton as the snakemake file
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
#shell.prefix(" set -eo pipefail; ")

#
# A little complicated to find the temp dir
# temp dir is a place programs can put temp files. 
# You will never have to change or think about anything here
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


configfile: "aln.yaml"

THREADS=4
if("threads" in config):
	THREADS = config["threads"]

SAVE=False
if("save" in config):
	SAVE = config["save"]

SPLIT=False
if("split" in config):
	SPLIT = config["split"]

# window size to calcualte coverage in 
WINDOW=1000

SAMPLES = list(config.keys())
if("threads" in SAMPLES): SAMPLES.remove("threads")
if("save" in SAMPLES): SAMPLES.remove("save")
if("split" in SAMPLES): SAMPLES.remove("split")
print("Samples: ", SAMPLES, file=sys.stderr) 

# samples that have BACs
REF_SMS = []
FOFN_SMS = []
READ_SMS = {}
CMD_SMS = []

RGNS = {}

for SM in SAMPLES:
	if("ref" in config[SM]):
		REF_SMS.append(SM)
	if("fofn" in config[SM] ):
		FOFN_SMS.append(SM)
		READ_SMS[SM] = [line.strip() for line in open(config[SM]["fofn"]).readlines()]
	if("regions" in config[SM]):
		t = [line.strip() for line in open(config[SM]["regions"]).readlines()]
		RGNS[SM] = [  "__".join(line.split())  for line in t ]
		

print("FOFNs: ", FOFN_SMS, file=sys.stderr)

def get_cmd(wildcards):
	SM = str(wildcards.SM)
	if("cmd" in config[SM]):	
		return(config[SM]["cmd"])
	return("minimap2 -t 8 --secondary=no -a --eqx -x map-pb -m 1000")

def get_flag(wildcards):
	SM = str(wildcards.SM)
	if("flag" in config[SM]):	
		return(config[SM]["flag"])
	return("2308")

def get_read(wildcards):
	SM = str(wildcards.SM)
	ID = int(str(wildcards.ID))
	return(READ_SMS[SM][ID])

def get_ref(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["ref"])

def get_bams(wildcards):
	SM = str(wildcards.SM)
	IDS = list(range(len(READ_SMS[SM])))
	bams = expand("temp/" + SM + "_{ID}.bam", ID=IDS)
	return(bams)

def get_all(wildcards):
	rtn  = []
	for SM in SAMPLES:
		IDS = list(range(len(READ_SMS[SM])))
		rtn += expand("temp/" + SM + "_{ID}.bam.bai", ID=IDS)
		rtn += expand("temp/" + SM + "_{ID}.bam.pbi", ID=IDS)
	return(rtn)
	

def my_temp(myf):
	if(SAVE):
		return(myf)
	return(temp(myf))

wildcard_constraints:
	ID="\d+",
	SM = "|".join(SAMPLES)

localrules: all

rule all:
	input:
		bam = expand("results/{SM}.bam", SM=SAMPLES),
		bai = expand("results/{SM}.bam.bai", SM=SAMPLES),

rule align:
	input:
		ref=get_ref,
		reads=get_read,
	output:
		bam=my_temp("temp/{SM}_{ID}.bam"),
	params:
		cmd = get_cmd,
		flag = get_flag,
	benchmark:
		"logs/align_{SM}_{ID}.b"
	resources:
		mem = 8,
		smem = 4
	threads: THREADS
	run:
		if( "minimap2" == params["cmd"].strip()[0:8] and ".bam" == input["reads"].strip()[-4:] ):
			shell("""
			SAM_TMP="{TMPDIR}/temp_{wildcards.SM}_{wildcards.ID}"
			rm -f $SAM_TMP*
			samtools fasta {input.reads} | \
				{params.cmd} \
				{input.ref} /dev/stdin | samtools view -F {params.flag} -u - | samtools sort -T $SAM_TMP -m {resources.smem}G -@ {threads} - > {output.bam}
			""")
		else:
			shell("""
			SAM_TMP="{TMPDIR}/temp_{wildcards.SM}_{wildcards.ID}"
			rm -f $SAM_TMP*
			{params.cmd} \
				{input.ref} {input.reads} | samtools view -F {params.flag} -u - | samtools sort -T $SAM_TMP -m {resources.smem}G -@ {threads} - > {output.bam}
			""")
	

if(not SPLIT and not SAVE):
	rule merge:
		input:
			bams = get_bams, 
		output:
			bam="results/{SM}.bam",
		benchmark:
			"logs/merge_{SM}.b"
		resources:
			mem = 4
		threads: min(16, 2*THREADS)
		shell:"""
	samtools merge -@ {threads} {output.bam} {input.bams}
	"""


	rule index:
		input:
			bam="results/{SM}.bam",
		output:
			bai="results/{SM}.bam.bai",
		benchmark:
			"logs/index_{SM}.b"
		resources:
			mem = 16
		threads: 1
		shell:"""
	samtools index {input.bam}
	"""

elif(SPLIT):
	rule index:
		input:
			bam="temp/{SM}_{ID}.bam",
		output:
			bai="temp/{SM}_{ID}.bam.bai",
		benchmark:
			"logs/index_{SM}_{ID}.b"
		resources:
			mem = 16
		threads: 1
		shell:"""
	samtools index {input.bam}
	"""
	rule pb_index:
		input:
			bam="temp/{SM}_{ID}.bam",
		output:
			pbi="temp/{SM}_{ID}.bam.pbi",
		benchmark:
			"logs/pb_index_{SM}_{ID}.b"
		resources:
			mem = 16
		threads: 1
		shell:"""
	pbindex {input.bam}
	"""
	
	rule fofn:
		input:
			bai=get_all,



#
# this rule creats a bed file that is incremented by 1000 for every contig
# these will be the feautes upon which we calculate depth wtih bedtools
#
rule fai_to_bed:
	input:
		ref=get_ref,
	output:
		regions="temp/{SM}.regions.bed",
	resources:
		mem=16
	threads: 1
	run:
		fai = open(input["ref"] + ".fai")
		window = WINDOW
		out = ""
		for line in fai:
				token = line.strip().split("\t")
				length = int(token[1])
				contig = token[0]
				for start in range(0, length, window):
						end = start + window -1
						if(end > length):
								end = length
						out += "{}\t{}\t{}\n".format(contig, start, end)
		open(output["regions"], "w+").write(out)


rule bam_to_coverage:
	input:
		bam=rules.index.input.bam,
		bai=rules.index.output.bai,
		regions=rules.fai_to_bed.output.regions,
	output:
		cov="results/{SM}.coverage.bed",
	resources:
		mem=16
	threads: 1
	shell:"""
# get coverage and then sort by contig and then pos
bedtools coverage -bed -mean -sorted -a {input.regions} -b {input.bam} > {output.cov}
#	sort -k 1,1 -k2,2n
"""


rule sm_coverage:
	input: 
		bed = "results/{SM}.coverage.bed"
	output:
		tbl = "results/{SM}.coverage.stats.txt",
		hcr = "results/{SM}.coverage.hcr.bed"
	resources:
		mem=16
	threads: 1
	run:
		df = pd.read_csv(input["bed"], sep="\t", names = ["contig", "start", "end", "cov"] ) 
		df["SM"] = str(wildcards.SM)
		#stats = pd.DataFrame(group["cov"].describe().rename(columns={'cov':sample}).squeeze() for sample, group in df.groupby('SM'))
		stats = pd.DataFrame(df["cov"].describe() ).transpose()
		print(stats)
		stats.to_csv(output["tbl"], sep = "\t", index=False)

		mean = stats["mean"][0]
		SD = np.sqrt( mean )
		top =  mean + SD*3
		bot =  mean - SD*3
		print(SD,mean, top,bot)

		df[df["cov"] >= top ].to_csv(output["hcr"], header=False, index=False, sep="\t")

rule coverage:
	input: 
		expand("results/{SM}.coverage.stats.txt", SM=SAMPLES),
	output:
		tbl = "results/coverage.stats.txt",
	resources:
		mem=16
	threads: 1
	run:
		dfs = []
		for SM, f in zip(SAMPLES, input):
			print(SM, f)
			tmp = pd.read_csv(f, sep="\t") 
			tmp["SM"] = SM
			dfs.append(tmp)
		df = pd.concat(dfs, ignore_index=True)
		print(df)
		df.to_csv(output["tbl"], sep = "\t", index=False)



	
#
#
#
def get_rgn(wildcards):
	SM = str(wildcards.SM)
	RG_ID = str(wildcards.RG_ID)
	rtn = "'{}:{}-{}'".format(*RG_ID.split("__"))
	return(rtn)

rule sm_regions:	
	input:
		bam=rules.index.input.bam,
		bai=rules.index.output.bai,
	output:
		png="results/regions/{SM}/{RG_ID}/{SM}_{RG_ID}.png",
		bam="results/regions/{SM}/{RG_ID}/{SM}_{RG_ID}.bam",
		bai="results/regions/{SM}/{RG_ID}/{SM}_{RG_ID}.bam.bai",
	params:
		rgn = get_rgn,
	resources:
		mem=16
	threads: 1
	shell:"""
samtools view -b {input.bam} {params.rgn} > {output.bam}
samtools index {output.bam}
/net/eichler/vol26/home/mvollger/projects/nucfreq/NucPlot.py {output.bam} {output.png}
"""

def get_rgns(wildcards):
	SM = str(wildcards.SM)
	fmt="results/regions/{SM}/{RG_ID}/{SM}_{RG_ID}.png"
	rgns = RGNS[SM]
	out = []
	for RG_ID in rgns:
		out.append(fmt.format(SM=SM, RG_ID=RG_ID))
	return(out)

rule merge_plots:
	input:
		get_rgns,
	output:
		"results/regions/{SM}.regions.pdf"
	resources:
		mem=16
	threads: 1
	shell:"""
convert {input} {output}
"""

	
rule plots:
	input:
		pdfs = expand("results/regions/{SM}.regions.pdf", SM=list(RGNS.keys())  )








	




