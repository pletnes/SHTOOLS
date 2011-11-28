#!/usr/bin/env python

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
    #print pySHTOOLS.plegendre(10, 0.32)

if __name__ == '__main__':
    main()

