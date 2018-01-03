from fsapt_generator import *
#from pymol import cmd
#
#def test_frag():
#    cmd.extend("fsapt",fsapt_generator)

#def test_writeFragFilesInFsaptDir():
#    import os
#    from fsapt_generator import fsapt_generator  
#    fsapt_generator()
#    assert 'fA.dat' in os.listdir('fsapt/')
#    assert 'fB.dat' in os.listdir('fsapt/')
#
#def test_bfs():
#    bfs()    
#
#def test_pm_mol_init():
#    phenol_dimer = pm_mol('examples/PhenolDimer.xyz') 


def test_frag():
    # A 12
    # B 24
    # C 2 4
    diphenyl = pm_mol('examples/diphenyl_benzene.xyz')
    diphenyl.cut(A = [12], B = [24], C = [2, 4])
    assert diphenyl.A == """
 C  -3.435   3.745  -0.197
 C  -3.818   2.469  -0.637
 C  -4.291   4.834  -0.341
 C  -5.084   2.315  -1.223
 C  -5.546   4.670  -0.931
 C  -5.940   3.406  -1.372
 H  -2.459   3.876   0.258
 H  -3.978   5.813   0.010
 H  -5.389   1.334  -1.577
 H  -6.212   5.521  -1.046
 H  -6.913   3.268  -1.836"""
    assert diphenyl.B == """
 C   0.922   3.945  -1.861
 C   0.298   2.913  -1.159
 C   0.325   4.472  -3.006
 C  -0.933   2.396  -1.587
 C  -0.898   3.960  -3.446
 C  -1.522   2.932  -2.743
 H   1.872   4.338  -1.509
 H   0.757   2.511  -0.260
 H   0.809   5.276  -3.553
 H  -1.366   4.361  -4.340
 H  -2.473   2.537  -3.085"""


