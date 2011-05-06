#!/usr/bin/env python

import numpy
import shtools

def main():
    lmax = 10
    z = 0.42
    print 'docstring:'
    print shtools.plmbar.__doc__
    plm_values = shtools.plmbar(lmax, z)
    print 'Plm after:', plm_values

if __name__ == '__main__':
    main()

