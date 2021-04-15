#!/usr/bin/env python
from setuptools import find_packages
from setuptools import setup

setup(
    name='Macrocomplex Builder',
    version='1.0',
    description='This program builds up and optimizes a model of protein macrocomplexes that can accept nucleic acids, given an imput of binary chains of the complex in pdb format',
    author='Mateusz Brodzik, Pol Ezquerra Condominas, Albert Garcia Valiente',
    author_email='mateusz.brodzik@estudiant.upf.edu, pol.ezquerra01@estudiant.upf.edu, albert.garcia11@estudiant.upf.edu',
    long_description=open('README.md').read(),
    packages=find_packages(),
    install_requires=['Biopython >= 1.78.0', 'argparse >= 1.1'],
    url='https://github.com/Albert-is-an-undefined-variable/ANNAM',
    scripts=['MB_builder','functions', 'post_optimization.py'])
