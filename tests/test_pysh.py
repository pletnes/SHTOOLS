#!/usr/bin/env python

import numpy as np
import pySHTOOLS

def main():
    """
    Run some tests...
    """
    assert pySHTOOLS.plmindex(0,0) == 0
    assert pySHTOOLS.plmindex(1,0) == 1
    assert pySHTOOLS.plmindex(1,1) == 2
    assert pySHTOOLS.plmindex(2,0) == 3
    assert pySHTOOLS.yilmindex(1, 0, 0) == 0
    pl = pySHTOOLS.plegendre(10, 0.32)
    assert (pl - np.array([1., 0.32, -0.3464, -0.39808, 0.0368752, 0.33970412, 0.16856375, -0.19099993, -0.26209324, 0.01135691, 0.24278892]) < 1.0e-8).all()

    plbar = pySHTOOLS.plbar(7, 0.32)
    print plbar
    plschmidt = pySHTOOLS.plschmidt(7, 0.32)
    print plschmidt
    plon = pySHTOOLS.plon(7, 0.32)
    print plon
    plmbar = pySHTOOLS.plmbar(7, 0.32)
    print plmbar
    plmschmidt = pySHTOOLS.plmschmidt(7, 0.32)
    print plmschmidt
    plmon = pySHTOOLS.plmon(7, 0.32)
    print plmon


if __name__ == '__main__':
    main()

