#!/usr/bin/env python

from distutils.core import setup

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('f2py_f90_ext',parent_package,top_path)
    config.add_extension('SHTOOLS',
                         ['src/pyshtools.f90'],
                         include_dirs=['include'],
                         f2py_options=['--include_paths',
                                       config.paths('include')[0]]
                         )
    config.add_data_dir('tests')
    return config

name = 'pySHTOOLS'
setup(name=name,
        version='3.beta',
        description='Python interface to the SHTOOLS package',
        author='Paul Anton Letnes',
        author_email='paul.anton.letnes@gmail.com',
        url='https://github.com/Dynetrekk/SHTOOLS'
        packages=[name, ],
    )




