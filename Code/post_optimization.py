from functions import *
import os

def post_optimization(pdb):
    """This function optimizes a pdb after being produced by ComplexMod"""
    print('Optimizing...')
    energies = optimize(pdb, os.path.abspath(os.getcwd()))
    print('Energy before optimizing: %s' % str(energies[0][0]))
    print('Energy after optimizing: %s' % str(energies[1][0]))
    print('Model completed')

if __name__ == '__main__':
    post_optimization(str(input('Write name pdb file: ')))
