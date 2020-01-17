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


configfile: "qc.yaml"

SAMPLES = list(config.keys())
print("Samples: ", SAMPLES, file=sys.stderr) 

# samples that have BACs
BACSMS = []
SDSMS = []
for SM in SAMPLES:
	if("bacs" in config[SM]):
		BACSMS.append(SM)
	if("SD" in config[SM] and config[SM]["SD"] ):
		SDSMS.append(SM)

print("Samples that have BACs: ", BACSMS, file=sys.stderr) 


def get_asm(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["asm"])

def get_bacs(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["bacs"])


MAXT=32

REF="/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta"
SD=f"{SNAKEMAKE_DIR}/data/ucsc.merged.max.segdups.bed"
SD10KB=f"{SNAKEMAKE_DIR}/data/ucsc.merged.10k.slop.segdups.bed"

localrules: all

rule all:
	input:
		index_to_hg38 = expand("results/{SM}.to.hg38.bam.bai", SM=SAMPLES),
		bed_to_hg38 = expand("results/{SM}.to.hg38.bed", SM=SAMPLES),
		# SD results
		sd_bed=expand("results/{SM}.SD.status.tbl", SM=SDSMS),
		sd_txt=expand("results/{SM}.SD.status.txt", SM=SDSMS),
		# BAC results 
		bac_tbl=expand("results/{SM}.bacs.to.asm.tbl", SM=BACSMS),
		qv_sum = "results/qv_sum.txt",


rule align_to_ref:
	input:
		ref=REF,
		asm=get_asm,
	output:
		bam="results/{SM}.to.hg38.bam",
	resources:
		mem = 4
	threads: MAXT
	shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	-m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
	 {input.ref} {input.asm} | samtools view -F 260 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}
"""

rule index_align_to_ref:
	input:
		bam="results/{SM}.to.hg38.bam",
	output:
		bai="results/{SM}.to.hg38.bam.bai",
	resources:
		mem = 8
	threads: 1
	shell:"""
samtools index {input.bam}
"""

rule bed_align_to_ref:
	input:
		bam="results/{SM}.to.hg38.bam",
	output:
		bed="results/{SM}.to.hg38.bed",
	resources:
		mem = 8
	threads: 1
	shell:"""
bedtools bamtobed -i {input.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}
"""

rule sd_align_to_ref:
	input:
		bed="results/{SM}.to.hg38.bed",
		sd=SD
	output:
		inter=temp("temp/{SM}.inter.bed"),
		rbed=temp("temp/{SM}.SD.resolved.bed"),
		ubed=temp("temp/{SM}.SD.unresolved.bed"),
		sd_bed="results/{SM}.SD.status.tbl",
	resources:
		mem = 8
	threads: 1
	shell:"""
bedtools intersect -a {input.bed} -b {input.sd} -wa -wb > {output.inter}
{SNAKEMAKE_DIR}/scripts/PrintResolvedSegdups.py {output.inter} {input.sd} \
	--extra 50000 \
	--resolved {output.rbed} \
	--unresolved {output.ubed} \
	--allseg {output.sd_bed} 
"""

rule ave_per_resolved:
	input:
		sd_bed="results/{SM}.SD.status.tbl",
	output:
		sd_txt="results/{SM}.SD.status.txt",
		sd_res_per="results/{SM}.SD.status.per.tbl",
	resources:
		mem = 8
	threads: 1
	run:
		thresh = list(range(-10000, 50000, 100))
		df = pd.read_csv(input["sd_bed"], header=None, names=["name","start","end","perid","ext"], sep="\t")
		total = sum(df.end-df.start)
		#print(df, total)
		perRes = []
		for t in thresh:
			res = sum( df[df.ext >= t ].end - df[df.ext >= t ].start)
			#print(res, total, res/total)
			perRes.append(res/total*100)
		rtn = "Average % of SDs resolved in {}: {}\n".format( wildcards["SM"] ,np.mean(perRes))
		print(rtn)
		df = pd.DataFrame({"ext":thresh, "per":perRes})
		df["SM"] = wildcards["SM"]
		df.to_csv(output["sd_res_per"], sep="\t", header=None, index=False)
		open(output["sd_txt"], "w+").write(rtn)


rule align_bac_to_ref:
	input:
		ref=REF,
		sd=SD10KB, 
		bacs=get_bacs,
	output:
		bam="temp/{SM}.bacs.to.hg38.bam",
		names="temp/{SM}.bacs.in.sd.names",
	resources:
		mem = 6
	threads: MAXT
	shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	-m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
	 {input.ref} {input.bacs} | samtools view -F 260 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}

bedtools intersect -a {output.bam} -b {input.sd} | samtools view - | awk '{{print $1}}' > {output.names}
"""


rule align_bac_to_asm:
	input:
		asm=get_asm,
		bacs=get_bacs,
	output:
		bam="results/{SM}.bacs.to.asm.bam",
	resources:
		mem = 6
	threads: MAXT
	shell:"""
minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	-m 10000 -z 10000,50 -r 500000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
	 {input.asm} {input.bacs} | samtools view -F 2308 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}
"""

rule bac_to_asm_tbl:
	input:
		bam="results/{SM}.bacs.to.asm.bam",
		names="temp/{SM}.bacs.in.sd.names",
	output:
		tbl="results/{SM}.bacs.to.asm.tbl",
	resources:
		mem = 6
	threads: 1
	shell:"""
{SNAKEMAKE_DIR}/scripts/samIdentity.py --header --tag {wildcards.SM} --mask {input.names} {input.bam} > {output.tbl}
"""

rule make_qv_sum:
	input:
		bac_tbl=expand("results/{SM}.bacs.to.asm.tbl", SM=BACSMS),
	output:
		qv_sum = "results/qv_sum.txt",
	resources:
		mem = 6
	threads: 1
	run:
		pd.options.mode.use_inf_as_na = True
		out = ""
		desc = []
		if( len(BACSMS) > 0 ): 
			for tbl in input["bac_tbl"]:
				sys.stderr.write(tbl + "\n")
				df = pd.read_csv(tbl, sep="\t")
				val = 1 - df["perID_by_all"]/100
				df["qv"] = -10 * np.log10( val )
				for mask in [False, True, "Total"]:
					if(mask == "Total"):	
						tmp = df
						tag = mask
					else:
						tmp = df[df["mask"] == mask]
						tag = "Unique"
						if(mask):
							tag = "SegDup"
					
					perfect = tmp["qv"].isna()
					stats = tmp.qv.describe()
					stats["Perfect"] = sum(perfect)
					stats["Status"] = tag
					desc.append(stats)
					#out += "{}\nPerfect\t{}\n{}\n\n".format(tbl, sum(perfect), tmp.qv.describe()   )
			
			
			open(output["qv_sum"], "w+").write(tbl + "\n" + str(pd.DataFrame(desc))  + "\n\n" )
		else:
			shell("touch {output.qv_sum}")






