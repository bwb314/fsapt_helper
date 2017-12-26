def test_frag():
    from pymol import cmd
    from fsapt_generator import fsapt_generator  
    cmd.extend("fsapt",fsapt_generator)

def test_writeFragFilesInFsaptDir():
    import os
    from fsapt_generator import fsapt_generator  
    fsapt_generator()
    assert 'fA.dat' in os.listdir('fsapt/')
    assert 'fB.dat' in os.listdir('fsapt/')
     
