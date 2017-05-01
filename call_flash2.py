from subprocess import call


def check_flash():
    found = 0
    for program in ('flash2', 'flash'):
        check = call('{0} --version'.format(program), shell=True)
        if check.returncode == 0:
            return program
    if found == 0:
        raise Exception('flash not found, please install it!')


call('flash2 {0} {1} -o {2} -d {3}'.format(
    F, R, arg.out_prefix, arg.out_directory) , shell=True)
