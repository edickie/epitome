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

def get_prereq(header):
    """
    Finds the prereq line of the header. Returns a list of them. Removes all
    empty values, so if there are no prereqs, this should return None.
    """
    for h in header:
        h = h.strip('# \n')
        if h.startswith('prereq:'):

            # remove empty fields, and the prereq tag
            h = filter(lambda x: '' != x, h.split(' '))
            h.remove('prereq:')
            if len(h) == 0:
                return None
            else:
                return(h)

def main(filename):
    header = get_header(filename)
    prereq = get_prereq(header)

if __name__ == '__main__':
    main(sys.argv[1])

