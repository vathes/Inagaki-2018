#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='Inagaki2018',
    version='0.0.0',
    description='Datajoint schemas for Inagaki 2018 2019 papers',
    author='Vathes',
    author_email='support@vathes.com',
    license='MIT',
    url='https://github.com/vathes/Inagaki-2018',
    packages=find_packages(exclude=[]),
    install_requires=['datajoint'],
)
