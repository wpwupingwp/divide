# divide.py

Divide NGS data by barcode and primer.

## Prerequisite

* Python 3.5 or above
* Biopython
* regex
* [vsearch](https://github.com/torognes/vsearch) (Optional)

To install Biopython and regex, run as administrator:

> pip install biopython regex

## Changelog

### v4.6
Support ambiguous base.

### v4.5
Extend vsearch options.
Improve output

### v4.2
Integrate vsearch.

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

Note that the forward and reverse barcode may be different sequence, but they
*SHOULD FOLLOW THE SAME MODE!*

## Strict option

Use "-s" or "--strict" to use strict version. If set, the program will check
barcode in head and tail is equal or not and whether barcode in tail (3') is
correct. If not, it will only check barcode in head (5') of sequence.

## Barcode file

Barcode file looks like this:

>   sample,barcode-f,barcode-r

>   S0001,ATACG,ATACG

>   S0002,ATATA,TATAC

>   S0003,ATACG

>   ...

The _barcode-f_ means barcode in 5' direction and _barcode-r_ means barcode in
3' direction. All sequences should be *forward*.

If forward and reverse barcode are same, you can omit the reverse barcode in
the table.

To avoid potential error, _please do not use space in sample info_.

And notice that here it use **English comma** to seperate two  fields rather
than **Chinese comma**.

## Primer file
Primer file looks like this:

>    gene,forward,reverse

>    rbcL,ATCGATCGATCGA,TACGTACGTACG

>    matK,AAAATTTTCCCC,GGGGTTACCAAAA

>    ...

Or:

>    gene,sequence

>    rbcL-f,ATCGATCGATCGA

>    rbcL-r,TACGTACGTACG


You can use Microsoft Excel to prepare these two files and save as CSV format,
or use any text editor you prefer.

**Make sure you don't miss the first line.**

# join_fastq.py

This program use "NNNNNNNNNN" to connect forward and reverse sequence pairs in
given two fastq files.

*Warning:*  Please do not use this tool if your sequences do have overlap but
FLASH failed to merge it. Usually it means those data has extremly poor
quality. If you merge them, you may get wrong result!

_Only used for sequences do not have overlap that FLASH cannot handle._

## Usage

> python3 join_fastq.py Forward.fastq Reverse.fastq

## Result

It will produce combined file as "combine.fastq". 

# task.sh

If you use PBS task submitting system, you can use this script to submit the
task, and you can finish the work from combine two direction sequence by flash and join_fastq.py to divide them.
