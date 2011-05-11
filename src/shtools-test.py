#!/usr/bin/env python

import numpy
from pyshtools import pyshtools

def main():
    lmax = 10
    z = 0.42
    print 'docstring:'
    print pyshtools.plmbar.__doc__
    plm_values = pyshtools.plmbar(lmax, z)
    print 'Plm after:', plm_values

if __name__ == '__main__':
    main()

