#!/bin/bash
set -euo pipefail

hifi_fofn=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env.cfg

threads=8
mem=128G

canu \
        -d asm \
        -p asm \
        genomeSize=3.1g \
        correctedErrorRate=0.015 \
        ovlMerThreshold=75 \
        batOptions="-eg 0.01 -eM 0.01 -dg 6 -db 6 -dr 1 -ca 50 -cp 5" \
        gridEngineThreadsOption="-pe serial THREADS" \
        gridEngineMemoryOption="-l m_mem_free=MEMORY" \
        gridOptions=" -R y -l gpfsstate=0 " \
		maxThreads=$threads \
        maxMemory=$mem \
        -pacbio-corrected $(cat $hifi_fofn) 

# useGrid=false		\

