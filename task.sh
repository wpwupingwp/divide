#!/bin/bash
#PBS -N merge
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00 
#PBS -q batch             
#PBS -S /bin/bash        

Workdir=/home/zhangxianchun/lichanghao/
cd $Workdir
L=$Workdir/L.fastq
R=$Workdir/R.fastq
BarcodeFile=$Workdir/barcode.csv
PrimerFile=$Workdir/primer.csv
Output=Result
time python3 divide.py $L $R -b $BarcodeFile -p $PrimerFile -o $Output
