#!/usr/bin/env python3

"""
Setup script for epystatic package.
"""

import setuptools

import epystatic


def main():
    setuptools.setup(
        version=epystatic.__version__,
    )
    # Rest of options are specified in `setup.cfg`


if __name__ == '__main__':
    main()
