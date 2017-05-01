# divide.py

Divide NGS data by barcode and primer.

## Sequence structure

>Barcode-adapter-primer-sequence

Till now, it only accept 5-2 repeat type of barcode.

## Barcode file

Barcode file looks like this:

>    barcode,sample

>    ATACG,BOP00001

**Notice that here it use English comma to seperate two fields. Make sure you use correct input method, especially for Chinese operating system.**

## Primer file

Primer file looks like this:

>    primer,gene,sequence

>    rbcL_F,rbcL,ATCGATCGATCGA

>    rbcL_R,rbcL,TACGTACGTACG

There is no limits for gene and primer name, except that you cannot use comma in them in order to avoid complicated seperator problem.

You can use Microsoft Excel to prepare these two files and save as CSV format, or use any text editor you prefer.

**Make sure you don't miss the first line.**

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
