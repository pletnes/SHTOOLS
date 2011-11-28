#!/usr/bin/env python

import pySHTOOLS

def main():
    """
    Run some tests...
    """
    print dir(pySHTOOLS)
    print pySHTOOLS.plmindex(0,0)
    print pySHTOOLS.plmindex(1,0)
    print pySHTOOLS.plmindex(1,1)
    print pySHTOOLS.plmindex(2,0)


if __name__ == '__main__':
    main()

