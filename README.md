# HiFi_Asm


## Exmaples for running ALIGNMENT

To run alignments script make a file called `aln.yaml` and then execute the script `align_snake.sh`. 

Here is an example `aln.yaml`:
```
threads: 8
ont:
    fofn: fofns/ont.fofn
    ref:  RefWithDecoy/patch_plus_decoy.fasta
    flag: 2308
    cmd: "minimap2 -t 8 --secondary=no -a --eqx -x map-ont -m 5000"
clr:
    fofn: fofns/clr.fofn
    ref:  RefWithDecoy/patch_plus_decoy.fasta
    flag: 2308
    cmd: "pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j 8 "
ccs:
    fofn: fofns/ccs.fofn
    ref:  RefWithDecoy/patch_plus_decoy.fasta
    flag: 2308
    cmd: "pbmm2 align --log-level DEBUG --preset CCS --min-length 5000 -j 8 "

```



## Exmaples for running CORRECTION

To run the correction script make a file called `cor.yaml` and then execute the script `correction_snake.sh`. 

Here is an example `cor.yaml`:

```
45xHiFi.chm13:
    alg: racon
    ref: asm/asm.contigs.fasta
    fofn: ccs.fastq.fofn
```

If you want to run correction twice you can simple use to output of one correction step as the input of the other. 

```
45xHiFi.chm13:
    alg: racon
    ref: asm/asm.contigs.fasta
    fofn: ccs.fastq.fofn
45xHiFi.chm13.2x:
    alg: racon
    ref: 45xHiFi.chm13_racon.fasta
    fofn: ccs.fastq.fofn
```


