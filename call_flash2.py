from subprocess import call
from sys import argv
import os


def flash_exist():
    for program in ['flash2', './flash2', 'flash', './flash']:
        check = call('{} --version'.format(program))
        if check.returncode == 0:
            return program
    raise Exception('FLASH not found!')


def main():
    F = argv[1]
    R = argv[2]
    call('flash2 {0} {1} -o {2}'.format(F, R, 'flash-merge'), shell=True)
    #     F, R, arg.out_prefix, arg.out_directory) , shell=True) #  test.extendedFrags.fastq test.notCombined_1.fastq test.notCombined_2.fastq
