#!/usr/bin/env python

import os, sys

def get_header(filename):
    """
    Returns each line of the input header as a list.
    """
    header = []
    f = open(filename, 'rb')
    for l in f.readlines():
        if l[0] == '#':
            header.append(l)
        else:
            break

    return header

def get_line(header, pattern):
    """
    Finds the prereq line of the header. Returns a list of them. Removes all
    empty values, so if there are no prereqs, this should return None.
    """
    for h in header:
        h = h.strip('# \n')
        if h.startswith(pattern):

            # remove empty fields, and the prereq tag
            h = filter(lambda x: '' != x, h.split(' '))
            h.remove(pattern)
            if len(h) == 0:
                return None
            else:
                return(h)

def main(filename):
    header = get_header(filename)
    prereq = get_line(header, 'prereq:')
    output = get_line(header, 'output:')
    args   = get_line(header, filename)

    if len(output) == 0:
        print('output is a copy of input')
    if len(output) > 1:
        output = output[0]
        print('ERROR: More than one output defined for {}, using {}'.format(filename, output[0]))

    print(prereq)
    print(output)
    print(args)

if __name__ == '__main__':
    main(sys.argv[1])

