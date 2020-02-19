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

# samples that have BACs/want SD run/have SDA collapses
BACSMS = []
SDSMS = []
COLSMS = []
for SM in SAMPLES:
	if("bacs" in config[SM]):
		BACSMS.append(SM)
	if("SD" in config[SM] and config[SM]["SD"] ):
		SDSMS.append(SM)
	if("collapse" in config[SM]):
		COLSMS.append(SM)

print("Samples that have BACs: ", BACSMS, file=sys.stderr) 


def get_asm(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["asm"])

def get_bacs(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["bacs"])

wildcard_constraints:
	SM="|".join(SAMPLES)	


MAXT=64

REF="/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta"
FAI="/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta.fai"
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
minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
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
minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
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
		if( len(BACSMS) > 0 ): 
			for tbl in input["bac_tbl"]:
				desc = []
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
				out += tbl + "\n" + str(pd.DataFrame(desc))  + "\n\n" 
			open(output["qv_sum"], "w+").write(out)

		else:
			shell("touch {output.qv_sum}")



rule make_ideogram:
	input:
		bed="results/{SM}.to.hg38.bed",
	output:
		pdf = "results/{SM}.ideogram.pdf"
	conda:
		"envs/karyo.yaml"
	resources:
		mem = 8
	threads: 1
	script:"""
Rscript /net/eichler/vol26/home/mvollger/projects/ideogram/ideogram.R  --asm {input.bed} --plot {output.pdf}
"""

rule make_ideograms:
	input:
		bed=expand("results/{SM}.to.hg38.bed", SM=SAMPLES),
	resources:
		mem = 8
	threads: 1






#
# Make a liftover chain
#
rule align_ref_to_asm:
	input:
		asm=REF,
		ref=get_asm,
	output:
		bam="results/hg38.to.{SM}.bam",
		bai="results/hg38.to.{SM}.bam.bai",
	resources:
		mem = 4
	threads: MAXT
	shell:"""
minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	-m 10000 -z 10000,50 -r 50000 -O 5,56 -E 4,1 -B 5 \
	 {input.ref} {input.asm} | samtools view -F 260 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}
samtools index {output.bam}
"""


rule make_chain:
	input:
		asm = get_asm ,
		bam="results/hg38.to.{SM}.bam",
		bai="results/hg38.to.{SM}.bam.bai",
	output:
		chain="results/{SM}.to.hg38.chain",
	shell:"""
{SNAKEMAKE_DIR}/scripts/bamtochain.py -i {input.bam} --fai {input.asm}.fai -o {output.chain}
"""


def get_collapse(wildcards):
	SM = str(wildcards.SM)
	return(config[SM]["collapse"])

rule liftover:
	input:
		bed = get_collapse,
		chain = rules.make_chain.output.chain, 	
	output:
		bed="results/{SM}.collapse.hg38.liftover.bed",
		unbed="results/{SM}.collapse.hg38.unlifted.bed",
	shell:"""
{SNAKEMAKE_DIR}/scripts/liftOver -minMatch=0.25 {input.bed} {input.chain} {output.bed} {output.unbed}
"""

rule collapse_to_ref:
	input:
		ref=REF,
		asm = get_asm,
		bed = get_collapse,
	output:
		fasta = temp("temp/{SM}.collapse.fasta"),
		bam=temp("temp/{SM}.collapse.to.hg38.bam"),
	resources:
		mem = 4
	threads: MAXT
	shell:"""
bedtools getfasta -fi {input.asm} -bed {input.bed} > {output.fasta}
minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	 -z 10000,50 -r 50000 -O 5,56 -E 4,1 -B 5 \
	 {input.ref} {output.fasta} | samtools view -F 2308 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}
"""

rule collapse_to_ref_bed:
	input:
		bam="temp/{SM}.collapse.to.hg38.bam",
	output:
		bed="temp/{SM}.collapse.to.hg38.bed",
	resources:
		mem = 4
	threads: 1
	shell:"""
bedtools bamtobed -i {input.bam} > {output.bed}
"""


rule col_merge:
	input:
		lo=rules.collapse_to_ref_bed.output.bed,
		bed = get_collapse,
		fai = FAI,
	output:
		bed="results/{SM}.collapse.to.hg38.bed",
	resources:
		mem = 8
	threads: 1
	run:
		cnames = ["chr", "start", "end", "mean", "median", "cmr", "length"]
		col = pd.read_csv(input["bed"], sep="\t", header=None, names=cnames)	
		col["qname"] = col["chr"] + ":" + col["start"].astype(str) + "-" +col["end"].astype(str)
		cnames = ["chr", "start", "end", "qname", "qual", "strand"]
		lo = pd.read_csv(input["lo"], sep="\t", header=None, names = cnames)	
		
		cnames = ["chr", "start", "end", "qname", "qual", "strand"]
		fai = pd.read_csv(input["lo"], sep="\t", header=None, names = cnames)	
			
		# merge dataframes 
		df = col.merge(lo, how="left", on = "qname", suffixes=("_asm", "_ref"))
		# reorder columns 
		df = df[["chr_ref", "start_ref", "end_ref"] + list(df)[0:6] ]
	
		# valid chromosomes 
		chrs = [  "chr{}".format(x) for x in list(range(1,23)) + ["X", "Y"] ]
		unknown = ~df.chr_ref.isin(chrs)
		df.chr_ref[ unknown ]= "Unknown"
		cur = 0
		for idx, row in df[unknown].iterrows():
			df.loc[idx, "start_ref"] = cur
			cur += row["end_asm"] - row["start_asm"]
			df.loc[idx, "end_ref" ]	= cur  
			cur += 1
		df[["start_ref", "end_ref"]] = df[["start_ref", "end_ref"]].astype(int)
		df.to_csv(output["bed"], sep="\t", index=False)


rule collapses:
	input:
		lo=expand(rules.col_merge.output.bed, SM=COLSMS),
	resources:
		mem = 8
	threads: 1






