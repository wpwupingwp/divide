# divide.py

Divide NGS data by barcode and primer.

## Changelog

### v4.0
Use regex instead of BLAST. Faster and easier.

### v3.3
Parallel version, use BLAST.

### v2.1
Single core version. Use BLAST.

### v1.0
Deprecated.

## Sequence structure

It can handle merged pair-end sequence like this:

>barcode-adapter-primer-sequence-primer-adapter-barcode

Or just handle one direction:

>barcode-adapter-primer-sequence

Sequences will be divided by barcode according to given barcode file.
If barcode is wrong even only one base, it will be dropped.

## adapter

Some one adds sequence between barcode and primer, if you do not have it, just
set adapter length to zero by "--adapter 0". The default value is 14.

## Barcode mode

Use "-m" to set barcode mode, like "8\*1", means barcode with length 5 repeats
only 1 times. The default is "5\*2", i.e., 5-base barcode repeats twice.

## Strict option

Use "-s" or "--strict" to use strict version. If set, the program will check
barcode in head and tail is equal or not and whether barcode in tail (3') is
correct. If not, it will only check barcode in head (5') of sequence.

## Barcode file

Barcode file looks like this:

>    barcode,sample

>    ATACG,BOP00001

>    ...

To avoid potential error, _please do not use space in sample info_.

And notice that here it use **English comma** to seperate two  fields rather
than **Chinese comma**.

## Primer file
Primer file looks like this:

>    gene,forward,reverse

>    rbcL,ATCGATCGATCGA,TACGTACGTACG

>    matK,AAAATTTTCCCC,GGGGTTACCAAAA

>    ...

You can use Microsoft Excel to prepare these two files and save as CSV format,
or use any text editor you prefer.

**Make sure you don't miss the first line.**


## Parallel

From v3.0, divide.py support parallel to speed up BLAST because of extremly
slow blastn-short. In default, it use all your CPU cores minus one to run, or
you could use "-c" to change number of cores you want to use.

# join_fastq.py

When you use flash or other software to combine sequence of two directions,
usually you got a lot of data could not be assemblied directly. This program
use "NNNNNNNNNN" to connect two direction.

## Usage

> python3 join_fastq.py Forward.fastq Reverse.fastq

## Result

It will produce combined file as "combine.fastq". 

# task.sh

If you use PBS task submitting system, you can use this script to submit the
task, and you can finish the work from combine two direction sequence by flash and join_fastq.py to divide them.
