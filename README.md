# HiFi_Asm

To run the correction script make a file called `cor.yaml` and then execute the script `correction_snake.sh`. 

Here is an example cor.yaml

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


