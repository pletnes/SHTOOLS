#!/usr/bin/env python

import numpy as np
import pySHTOOLS

def main():
    """
    Run some tests...
    """
    lmax = 7
    xval = 0.32
    assert pySHTOOLS.plmindex(0,0) == 0
    assert pySHTOOLS.plmindex(1,0) == 1
    assert pySHTOOLS.plmindex(1,1) == 2
    assert pySHTOOLS.plmindex(2,0) == 3
    assert pySHTOOLS.yilmindex(1, 0, 0) == 0
    # Legendre functions
    pl = pySHTOOLS.plegendre(10, xval)
    assert (pl - np.array([1., 0.32, -0.3464, -0.39808, 0.0368752, 0.33970412, 0.16856375, -0.19099993, -0.26209324, 0.01135691, 0.24278892]) < 1.0e-8).all()
    plbar = pySHTOOLS.plbar(lmax, xval)
    print plbar
    plschmidt = pySHTOOLS.plschmidt(lmax, xval)
    print plschmidt
    plon = pySHTOOLS.plon(lmax, xval)
    print plon
    p, dp = pySHTOOLS.plbar_d1(lmax, xval)
    print ' plbar=', p, 'd_plbar=', dp
    p, dp = pySHTOOLS.plschmidt_d1(lmax, xval)
    print ' plschmidt=', p, 'd_plschmidt=', dp
    p, dp = pySHTOOLS.plon_d1(lmax, xval)
    print ' plon=', p, 'd_plon=', dp
    # Associated Legendre functions
    plmbar = pySHTOOLS.plmbar(lmax, xval)
    print plmbar
    plmschmidt = pySHTOOLS.plmschmidt(lmax, xval)
    print plmschmidt
    plmon = pySHTOOLS.plmon(lmax, xval)
    print plmon
    p, dp = pySHTOOLS.plegendrea_d1(lmax, xval)
    print ' plegendre=', p, 'd_plegendre=', dp




if __name__ == '__main__':
    main()

