#!/bin/sh
#PBS -N merge
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00 
#PBS -q batch             
#PBS -S /bin/bash        

Workdir=/home/user/data/
cd $Workdir
L=$Workdir/L.fastq
R=$Workdir/R.fastq
BarcodeFile=$Workdir/barcode.csv
PrimerFile=$Workdir/primer.csv
Merged=$Workdir/merged.fastq
cd $Workdir
time flash -t 8 -M 200 $L $R
time python3 join_fastq.py out.notCombined_1.fastq out.notCombined_2.fastq
cat out.extendedFrags.fastq combine.fastq > $Merged
time python3 divide.py $Merged -b $BarcodeFile -p PrimerFile -s
