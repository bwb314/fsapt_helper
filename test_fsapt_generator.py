from fsapt_generator import *
from pymol import cmd

def test_frag():
    cmd.extend("fsapt",fsapt_generator)

#def test_writeFragFilesInFsaptDir():
#    import os
#    from fsapt_generator import fsapt_generator  
#    fsapt_generator()
#    assert 'fA.dat' in os.listdir('fsapt/')
#    assert 'fB.dat' in os.listdir('fsapt/')

def test_bfs():
    bfs()    

def test_pm_mol_init():
    phenol_dimer = txt_molecule('examples/PhenolDimer.xyz') 
