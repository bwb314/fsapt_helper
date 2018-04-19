This tool only works for Python 2.7 at the moment due to PyMOL's 2.7 dependencies.

To generate an F-ISAPT input and F-SAPT fragment files, 

0) Open an XYZ file of the geometry that is to be fragmented

1) Select border atoms of fragment. Border atoms are atoms connected
to other atoms that are not part of the fragment.

2) Name this selection 'fragmentName_X', where X is the monomer the
fragment belongs to ('A', 'B', or 'C').

3) Repeat 1 and 2 until all fragments have been selected. If there are
atoms unaccounted for (i.e. solvent) after the fragment selection has
occurred, these will automatically be grouped into a fragment named
'ISAPT_C'.

4) Type 'run path/to/fsapt_generator.py' into PyMOL console.

5) Type 'fisapt' into PyMOL console. This should produce a visual
representation of the fragmentation scheme, an input file with the same
name as the XYZ file except with a '.in' extension, and 'fA.dat' and
'fB.dat' fragment files.

WARNING: This program will overwrite any preexisting fragment files
or input files of the same name.
