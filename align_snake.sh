#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env.cfg

#
# snakemake paramenters
#
snakefile=$DIR/align.smk
jobNum=600
waitTime=60 # this really needs to be 60 on our cluster :(
retry=0 # numer of times to retry the pipeline if it failes
# I allow a retry becuase sometimes even the really long waittime is not enough,
# and the files are actaully there

#
# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
# Allow snakemake to make directories, I think it slows things down when I done with "waitTime"
#
logDir=logs
mkdir -p $logDir

#
# run snakemake
#

# make a new pipe for tee
exec 3>&1 

# unbuffer keeps the colors for tee
unbuffer snakemake -p \
        -s $snakefile \
        --drmaa " -P eichlerlab \
                -q eichler-short.q \
                -l h_rt=48:00:00  \
                -l mfree={resources.mem}G \
				-pe serial {threads} \
                -V -cwd \
                -S /bin/bash" \
		 --drmaa-log-dir $logDir \
        --jobs $jobNum \
        --latency-wait $waitTime \
        --restart-times $retry  \
         $@ 2>&1 >&3 | tee snakemake.aln.stderr.txt

# sends an email to the user that we are done
mail -s "snakemake alignment finished" $USER@uw.edu  < snakemake.aln.stderr.txt


