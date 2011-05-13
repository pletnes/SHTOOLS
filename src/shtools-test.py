#!/usr/bin/env python

import numpy
from pyshtools import pyshtools
print 'pyshtools contents:'
print dir(pyshtools)

def main():
    lmax = 10
    z = 0.42
    print 'docstring:'
    print pyshtools.plmbar.__doc__

    plm_bar_values = pyshtools.plmbar(lmax, z)
    print 'Plm bar:', plm_bar_values

    plm_schmidt_values = pyshtools.plmschmidt(lmax, z)
    print 'Plm Schmidt:', plm_schmidt_values

    plm_on_values = pyshtools.plmon(lmax, z)
    print 'Plm Schmidt:', plm_on_values

    print 'PLegendre:', pyshtools.plegendre(lmax, z)

    pl_bar_values = pyshtools.plbar(lmax, z)
    print 'Pl Schmidt:', pl_bar_values

    pl_schmidt_values = pyshtools.plschmidt(lmax, z)
    print 'Pl Schmidt:', pl_schmidt_values
    
    pl_on_values = pyshtools.plon(lmax, z)
    print 'Pl Schmidt:', pl_on_values

    theta = 0.42
    phi = 0.9
    ylm_values = pyshtools.ylm(lmax, theta, phi)
    print 'Ylm:', ylm_values

    l = 9
    m = 7
    print pyshtools.plmindex(l, m)

    i = 2
    print pyshtools.yilmindex(i, l, m)


if __name__ == '__main__':
    main()

