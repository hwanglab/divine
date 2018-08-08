#!/usr/bin/python
from subprocess import Popen, PIPE, STDOUT
import sys
import os
import warnings


def get_command_output(cmd):

    """
    Captures and returns the output of cmd passed in
    """
    return Popen(cmd.split(), stdout=PIPE, stderr=STDOUT).communicate()[0]


def hook():

    """
    Validates the file against the PEP8 guidelines
    See http://www.python.org/dev/peps/pep-0008/
    """
    changes = os.system("hg status > /var/tmp/out.txt")
    files = []
    fp = open("/var/tmp/out.txt", 'r')
    for line in fp.readlines():
        files.append(line.strip().split(" ")[1])
    fp.close()

    errors = []
    for pyfile in files:
        if pyfile[-2:] != 'py':
            continue
        cmd = "/usr/bin/pep8 " + pyfile.strip()
        output = get_command_output(cmd)
        if output:
            for error in output.split("\n")[:-1]:
                errors.append(error)
    if errors:
        warnings.warn("\n\nThe following files are not PEP8 compliant:\n" + \
                       "\n".join(errors) + "\n")
        return -1
    return 0

if __name__ in "__main__":
        status = hook()
        print "\npep8 compliance status ", status
        if status == -1:
            sys.exit(1)
