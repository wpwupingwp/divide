from subprocess import call

call('flash2 {0} {1} -o {2} -d {3}'.format(
    F, R, arg.out_prefix, arg.out_directory) , shell=True)
