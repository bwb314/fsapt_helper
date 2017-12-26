import os

def fsapt_generator():
    
    os.mkdir('fsapt/')
    with open('fsapt/fA.dat','w') as fil:
        fil.close()    
    with open('fsapt/fB.dat','w') as fil:
        fil.close()    

#import os, sys, math, commands, string
#export PYTHONPATH=/anaconda/bin/:/Users/brandonbakr/psi4/install/lib/
#from pymol import cmd
#
##deprecated atm
#def write_frags(frags,base_name):
#    for frag in range(len(frags)):
#        fil = open(os.getcwd()+'/'+base_name+'_'+string.ascii_uppercase[frag]+'.xyz','w')
#        fil.write(str(len(frags[frag]))+'\n\n')
#        for lin in frags[frag]:
#            fil.write(lin)
#        fil.close()
#
##read main geometry
#def read_original_geometry(fil_name):
#    geom = []
#    for i in open(os.getcwd()+'/'+fil_name,'r'):
#        if len(i.split()) == 4: geom.append(i)
#    return geom
#
##read fragmentation input for eventual fA and fB generation
#def read_fragment_files():
#    frags = []
#    for i in cmd.get_names("all")[1:]:
#        #pymol is bugging out atm, ignore selection "_drag"
#        if i == '_drag': continue
#        if not i.endswith('_A') and not i.endswith('_B') and not i.endswith('_a') and not i.endswith('_b'): 
#            print "End custom fragment names with '_A' or '_B'."
#            print_usage()
#            sys.exit()
#        stored.list=[]
#        cmd.iterate("("+i+")","stored.list.append((name,rank))")
#        frag = [i]
#        for atom in stored.list:
#            frag.append(atom[1])
#        frags.append(frag)
#    return frags
#
##write fragmentation input for fA and fB generation
#def write_prefrag_file(frags):
#    fragfile = open(os.getcwd()+'/'+'fraginfo.dat','w')
#    for frag in frags:
#        fragfile.write(str(frag[0])+' ')
#        for el in frag[1:]:
#            fragfile.write(str(el+1)+' ')
#        fragfile.write('\n')
#    fragfile.close()
#
##write fragment geometries and move them to check_geoms/ directory 
#def write_frag_geoms():
#    newgeomz = open(os.getcwd()+'/'+'Geomz.xyz','r')   
#    fil = ''
#    fils = []
#    for lin in newgeomz:
#        if len(lin.split()) == 1:
#            if fil != '': fil.close() 
#            num = lin
#            name = next(newgeomz)
#            fil = open(os.getcwd()+'/'+name.split()[1]+'.xyz','w')
#            fils.append(name.split()[1]+'.xyz') 
#            fil.write(num)
#            fil.write(name)
#        if len(lin.split()) == 4: 
#            fil.write(lin)
#    fil.close()
#    os.system("mkdir check_geoms")
#    for fil in fils: os.system("mv "+fil+" check_geoms/")
#
##map atoms from monomers to abs position in input and fA and fB
#def get_maps_for_separating(geom):
#    clone = [x for x in range(len(geom))]
#    fA = []
#    #frag inds starting from 0
#    Aatoms = []
#    for lin in open(os.getcwd()+'/'+'fA.dat','r'):
#        fA.append(lin.split())
#    map2A = {} 
#    for frag in fA:
#        for atom in range(len(frag[1:])):
#            frag[1+atom] = int(frag[1+atom]) - 1
#            Aatoms.append(frag[1+atom])
#            map2A[frag[1+atom]] = len(Aatoms)-1
#            clone.remove(int(frag[1+atom]))
#    fB = []
#    #frag inds starting from 0
#    Batoms = []
#    for lin in open(os.getcwd()+'/'+'fB.dat','r'):
#        fB.append(lin.split())
#    map2B = {}
#    for frag in fB:
#        for atom in range(len(frag[1:])):
#            frag[1+atom] = int(frag[1+atom]) - 1
#            Batoms.append(frag[1+atom])
#            map2B[frag[1+atom]] = len(Batoms)-1
#            clone.remove(int(frag[1+atom]))
#    Catoms = []
#    if clone != []:
#        for i in clone: 
#            Catoms.append(i)
#
#    os.system('rm fA.dat')
#    os.system('rm fB.dat')
#    if sorted(Aatoms)[0] > sorted(Batoms)[0]:
#        datoms = list(Aatoms)
#        Aatoms = list(Batoms)
#        Batoms = list(datoms)
#        fd = list(fA)
#        fA = list(fB)
#        fB = list(fd)
#        for frag in fA:
#            frag[0] = ''.join(frag[0].split('_')[:-1])+'_A'
#        for frag in fA:
#            frag[0] = ''.join(frag[0].split('_')[:-1])+'_B'
#        map2d = map2A.copy()
#        map2A = map2B.copy()
#        map2B = map2d.copy()
#    map2orig = {}
#    new_geom = []
#    for a in range(len(Aatoms)):
#        new_geom.append(geom[Aatoms[a]])
#    for b in range(len(Batoms)):
#        new_geom.append(geom[Batoms[b]])
#    for c in range(len(Catoms)):
#        new_geom.append(geom[Catoms[c]])
#    map2orig = {}
#    for a in range(len(geom)):
#        for b in range(len(new_geom)): 
#            if geom[a] == new_geom[b]:
#                map2orig[b] = a 
#    return Aatoms, map2A, fA, Batoms, map2B, fB, Catoms, map2orig
#
##Write input file
#def write_input_file(Aatoms, Batoms, Catoms, geom, fil_name, fsapt_path, fsapt_plot_path, inp_fil, template):
#    if inp_fil == '': inp_name = os.getcwd()+'/'+fil_name.split('.xyz')[0]+'.in'
#    else: inp_name = os.getcwd() + '/' + inp_fil
#    inp = open(inp_name,'w')
#    if template == '':
#        inp.write('memory 4 Gb\n\n')
#        inp.write('molecule mol {\n')
#        #inp.write('0 1 #make sure this is accurate for your system\n')
#        inp.write('CHARGE MULTIPLICITY\n')
#        for a in sorted(Aatoms):
#            inp.write(geom[a])
#        #inp.write('--\n0 1 #make sure this is accurate for your system\n')
#        inp.write('--\nCHARGE MULTIPLICITY\n')
#        for b in sorted(Batoms):
#            inp.write(geom[b])
#        if len(Catoms) > 0:
#            inp.write('--\nCHARGE MULTIPLICITY\n')
#            for c in sorted(Catoms):
#                inp.write(geom[c])
#        inp.write('no_com\n')
#        inp.write('no_reorient\n')
#        inp.write('symmetry c1\n')
#        inp.write('}\n\n')
#        inp.write('set basis jun-cc-pvdz\n')
#        inp.write('set df_basis_scf  jun-cc-pvdz-jkfit\n')
#        inp.write('set df_basis_sapt jun-cc-pvdz-ri\n')
#        inp.write('set scf_type df\n')
#        inp.write('set guess sad\n')
#        inp.write('set d_convergence 9\n')
#        inp.write('set freeze_core true\n')
#        if fsapt_plot_path != '': 
#            inp.write('set FISAPT_DO_PLOT TRUE\n')
#            inp.write('set FISAPT_PLOT_FILEPATH '+fsapt_plot_path+'\n')
#        if fsapt_path != 'fsapt': inp.write('set FISAPT_FSAPT_FILEPATH '+fsapt_path+'\n\n')
#        inp.write("energy('fisapt0')\n")
#        inp.close()
#    elif template not in os.listdir(os.getcwd()):
#        print "Please move "+template+" to the current working" 
#        print "directory in order to use as template."
#        sys.exit()
#    else:
#        inpfil = []
#        template_fil = open(os.getcwd()+'/'+template,'r')
#        for lin in template_fil: 
#            if '--' in lin:
#                for a in sorted(Aatoms):
#                    inpfil.append(geom[a])
#                inpfil.append(lin)
#                inpfil.append(next(template_fil))
#                for b in sorted(Batoms):
#                    inpfil.append(geom[b])
#                if len(Catoms) > 0:
#                    inpfil.append(next(template_fil))
#                    inpfil.append(next(template_fil))
#                    for c in sorted(Catoms):
#                        inpfil.append(geom[c])
#                continue    
#            inpfil.append(lin)    
#    
#        template_fil.close()
#        for lin in inpfil:
#            if 'FISAPT_FSAPT_FILEPATH' in lin.upper(): fsapt_path = lin.split()[-1]
#            inp.write(lin)
#        inp.close()
#    return fsapt_path
#
##rewrite fA and fB with new geometry
#def write_fA_and_fB(map2A, fA, map2B, fB, Aatoms, fsapt_path):
#    os.system('mkdir '+os.getcwd()+'/'+fsapt_path)
#    fAfil = open(os.getcwd()+'/'+fsapt_path+'/fA.dat','w')
#    for frag in fA:
#        fAfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fAfil.write(str(map2A[frag[1+atom]]+1)+' ')
#        fAfil.write('\n')
#    fAfil.close()
#    fBfil = open(os.getcwd()+'/'+fsapt_path+'/fB.dat','w')
#    for frag in fB:
#        fBfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fBfil.write(str(map2B[frag[1+atom]]+len(Aatoms)+1)+' ')
#        fBfil.write('\n')
#    fBfil.close()
#
##rewrite the original fA and fB with new frag scheme chosen post-calc
#def try_rewrite_fA_and_fB(map2orig, map2A, fA, map2B, fB, Aatoms, fsapt_path):
#    os.system('mkdir '+os.getcwd()+'/'+fsapt_path)
#    fAfil = open(os.getcwd()+'/'+fsapt_path+'/fA.dat','w')
#    for frag in fA:
#        fAfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fAfil.write(str(map2orig[map2A[frag[1+atom]]]+1)+' ')
#        fAfil.write('\n')
#    fAfil.close()
#    fBfil = open(os.getcwd()+'/'+fsapt_path+'/fB.dat','w')
#    for frag in fB:
#        fBfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fBfil.write(str(map2orig[map2B[frag[1+atom]]+len(Aatoms)]+1)+' ')
#        fBfil.write('\n')
#    fBfil.close()
#
##rewrite the original fA and fB with new frag scheme chosen post-calc
#def rewrite_fA_and_fB(map2orig, map2A, fA, map2B, fB, Aatoms):
#    fAfil = open(os.getcwd()+'/'+'fA.dat','w')
#    for frag in fA:
#        fAfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fAfil.write(str(map2orig[map2A[frag[1+atom]]]+1)+' ')
#        fAfil.write('\n')
#    fAfil.close()
#    fBfil = open(os.getcwd()+'/'+'fB.dat','w')
#    for frag in fB:
#        fBfil.write(frag[0]+' ')
#        for atom in range(len(frag[1:])):
#            fBfil.write(str(map2orig[map2B[frag[1+atom]]+len(Aatoms)]+1)+' ')
#        fBfil.write('\n')
#    fBfil.close()
#    
#
##remove files kept for reference
#def clean_up():
#    os.system('rm Geomz.xyz')        
#    os.system('rm fraginfo.dat')       
#
##'distance' between colors
#def dist(color1, color2):
#    r1,g1,b1 = color1
#    r2,g2,b2 = color2
#    dist = ((r1-r2)**2+(g1-g2)**2+(b1-b2)**2)**0.5
#    return dist
#
##check for distinctness of color
#def distinct(color,colorlist):
#    for c in colorlist:
#        if dist(color,c) < 0.5/len(colorlist): return False
#    return True
#
##color fragments based on fragmentation from bfs
#def color_frags(map2orig, map2A, fA, map2B, fB, Aatoms, Catoms):
#    warm = []
#    cool = []
#    x = math.ceil((len(fA)/8.)**(1./3.))
#    for li in range(int(2*x)):
#        i = 1.0 - li*(1.0-0.6)/(2*x-1)
#        for lj in range(int(4*x)):
#            j = lj*(1.0-0.0)/(4*x-1)
#            for lk in range(int(x)):
#                safe = x
#                if x != 1: safe = x - 1
#                k = lk*(0.25-0.0)/(safe)
#                warm.append([i,j,k])
#    
#    x = math.ceil((len(fB)/8.)**(1./3.))
#    for lk in range(int(2*x)):
#        k = 1.0 - lk*(1.0-0.6)/(2*x-1)
#        for lj in range(int(4*x)):
#            j = lj*(1.0-0.0)/(4*x-1)
#            for li in range(int(x)):
#                safe = x
#                if x != 1: safe = x - 1
#                i = li*(0.25-0.0)/(safe)
#                cool.append([i,j,k])
#
#    for frag in range(len(fA)):
#        selection = ''
#        for atom in range(len(fA[frag][1:])):
#            key = str(map2orig[map2A[fA[frag][1+atom]]])
#            selection += 'rank '+key+' '                
#        cmd.select("("+fA[frag][0]+")", selection)
#        #linspace through color vector
#        safefA = len(fA)
#        if len(fA) != 1: safefA = len(fA) - 1
#        fragind = int(math.floor(frag*((len(warm))-1)/safefA))
#        color = warm[fragind]
#        cmd.set_color("A_"+fA[frag][0], color)
#        cmd.color("A_"+fA[frag][0],"("+fA[frag][0]+")")
#    for frag in range(len(fB)):
#        selection = ''
#        for atom in range(len(fB[frag][1:])):
#            key = str(map2orig[map2B[fB[frag][1+atom]]+len(Aatoms)])
#            selection += 'rank '+key+' '           
#        cmd.select("("+fB[frag][0]+")", selection)
#        #linspace through color vector
#        safefB = len(fB)
#        if len(fB) != 1: safefB = len(fB) - 1
#        fragind = int(math.floor(frag*(len(cool)-1)/safefB))
#        color = cool[fragind]
#        cmd.set_color("B_"+str(fB[frag]), color)
#        cmd.color("B_"+str(fB[frag]),"("+fB[frag][0]+")")
#    if len(Catoms) > 0:
#        selection = ''
#        for atom in Catoms:
#            selection += 'rank '+str(atom)+' '
#        cmd.select("(ISAPT_C)", selection)
#        cmd.color("grey", "(ISAPT_C)")
#        
#
#def print_usage():
#    print "WARNING: Any 'fA.dat' or 'fB.dat' files within"
#    print "the current working directory will be removed"
#    print "upon successful execution of this script." 
#    print "Please move them to a safe place before"
#    print "continuing."
#    print ""
#    print "This program is meant to be used within Pymol."
#    print ""
#    print "Directions (if making F-SAPT input file):"
#    print "1. Change to the directory containing the"
#    print "   desired geometry with .xyz file extension"
#    print "   and open the geometry with Pymol."
#    print "  -Example: 'pymol H2Odimer.xyz'"
#    print "2. In Pymol, make a selection of all of the" 
#    print "   border atoms for each F-SAPT fragment."
#    print "  -Name each selection with an appropriate"
#    print "   fragment name plus '_A' or '_B' for the"
#    print "   corresponding monomer (A or B). This"
#    print "   designation is case insensitive."
#    print "  -Example: Methyl_A, phenyl_B"
#    print "3. In the Pymol command line, type"
#    print "   'run /FULL/PATH/TO/PROGRAM/geom_separator.py'"
#    print "   -Example; 'run /Users/john_doe/Desktop/frag_gen/geom_separator.py"
#    print ""
#    print "That is it! Now an input file has been"
#    print "generated, along with 'fA.dat' and 'fB.dat'"
#    print "for F-SAPT analysis. Move the .dat files"
#    print "into fsapt/ and/or s-fsapt/ when ready."
#    print ""
#    print "The input name will have the same name"
#    print "as the '.xyz' file but with a '.in' extension."
##    print ""
##    print "Directions (if just making fA.dat and fB.dat)"
##    print "1. Open the geometry used for F-SAPT input"
##    print "   by typing 'pymol geom.xyz -- just_fragment'."
##    print "2. In Pymol, make separate selections of border" 
##    print "   atoms for F-SAPT fragmentation."
##    print "  -Name each selection with an appropriate"
##    print "   fragment name plus '_A' or '_B' for the"
##    print "   corresponding monomer (A or B)."
##    print "  -Example: Methyl_A, phenyl_B"
##    print "3. In the Pymol command line, type"
##    print "   'run fsapt_filemaker.py'"
##    print ""
##    print "Now you will have fA.dat and fB.dat for a"
##    print "new round of F-SAPT analysis. If these files"
##    print "already exist in the cwd, they will be"
##    print "replaced by the new files."
##    print "" 
##    print "WARNING: Be sure that the ranks of atoms for"
##    print "monomer A and B used for F-SAPT analysis are"
##    print "nonoverlapping, meaning rank of atoms in"
##    print "fragments for monomer A are in range [0,A]"
##    print "and the rank of atoms in fragments for"
##    print "monomer B are in range [A+1,Total Atoms - 1]."
##    print "In short, use geom.xyz from fsapt/ or s-fsapt/"
##    print "as stated earlier because this preserves the"
##    print "nonoverlapping condition as well as having the"
##    print "atom rank correspond to what was used in the" 
##    print "original F-SAPT job."
##    print "" 
##    print "" 
##    print "" 
#
#def sanity_check(fil_name):
#    if not fil_name.endswith('.xyz'): 
#        print "Use .xyz extension for geometry!"
#        print_usage() 
#        sys.exit()
#    if len(sys.argv) == 0:
#        print_usage() 
#        sys.exit()
#
####The bfs function was mostly written by Trent M. Parker and modified by Brandon W. Bakr.
#def bfs(fil_name,fraginfo):
#
#    ## CONSTANTS ##
#    
#    # factor beyond average of covalent radii to determine bond cutoff
#    bond_thresh = 1.20
#    
#    # covalent (or ionic) radii by atomic element in angstroms from
#    # "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
#    cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
#      'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
#      'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
#      'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
#      'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
#      'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
#      'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}
#    
#    ## FUNCTIONS ##
#    
#    # create 2d array of strings from input file name
#    def get_file_string_array(file_name):
#      try:
#        file = open(os.getcwd()+'/'+file_name, "r")
#      except IOError:
#        print 'Error: file (%s) does not exist!' % (file_name)
#        sys.exit()
#      file_data = file.readlines()
#      file.close()
#      file_array = []
#      for line in file_data:
#        file_array.append(line.split())
#      return file_array
#    
#    # return key string from x, y, z values and block resolution
#    def get_key(x, y, z, b):
#      return '%i,%i,%i' % (x - x%b, y - y%b, z - z%b)
#    
#    # square distance between two 3-d cartesian coordinates
#    def get_r2ij(coords_i, coords_j):
#      r2ij = 0.0
#      for p in range(3):
#        r2ij += (coords_j[p] - coords_i[p])**2
#      return r2ij
#    
#    ## CLASSES ##
#    
#    # geometry information for a molecule
#    class Geometry:
#      # constructor
#      def __init__(self, xyz_file_name):
#        self.source = xyz_file_name
#        self.read_xyz(self.source)
#    
#      # read in data from file
#      def read_xyz(self, xyz_file):
#        xyz_array = get_file_string_array(xyz_file)
#        self.n_atoms = int(xyz_array[0][0])
#        self.comment = ' '.join(xyz_array[1])
#        self.at_types = ['' for i in range(self.n_atoms)]
#        self.coords = [[0.0 for j in range(3)] for i in range(self.n_atoms)]
#        for i in range(self.n_atoms):
#          self.at_types[i] = xyz_array[i+2][0].capitalize()
#          for j in range(3):
#            self.coords[i][j] = float(xyz_array[i+2][j+1])
#        self.get_covradii()
#        self.bond_tree = []
#        self.fragments = []
#        self.fragment_names = []
#    
#      # print geometry to output file
#      def print_xyz(self, out_file, print_n_atoms):
#        if (print_n_atoms):
#          out_file.write('%i\n%s\n' % (self.n_atoms, self.comment))
#        for i in range(self.n_atoms):
#          out_file.write(' %-2s' % (self.at_types[i]))
#          for j in range(3):
#            out_file.write(' %12.6f' % (self.coords[i][j]))
#          out_file.write('\n')
#    
#      # print fragment geometry to screen
#      def print_xyz_frag(self, print_n_atoms, frag_num, geomfile):
#        fil = open(geomfile,'w')
#        n_frag_atoms = len(self.fragments[frag_num])
#        fil.write('%i\n  fragment' % (n_frag_atoms))
#        if (frag_num < len(self.fragment_names)):
#          fil.write('%s\n' % (self.fragment_names[frag_num]))
#        else:
#          fil.write('%i\n' % (frag_num+1))
#        for q in range(n_frag_atoms):
#          i = self.fragments[frag_num][q]
#          fil.write(' %-2s' % (self.at_types[i]))
#          for j in range(3):
#            fil.write('%12.6f' % (self.coords[i][j]))
#        fil.close()    
#
#      # print geometry of all fragments to screen
#      def print_xyz_allfrags(self, print_n_atoms):
#        for p in range(len(self.fragments)):
#          self.print_xyz_frag(print_n_atoms, p, 'Geomz.xyz')
#    
#      # determine molecular fragment structure from bond tree
#      def get_fragments(self, intrafrag_file):
#        self.get_bond_tree()
#        self.get_intrafrags(intrafrag_file)
#        break_list = []
#        new_list = []
#        unfound_list = [i for i in range(self.n_atoms)]
#        for p in range(len(self.fragments)):
#          new_list.append([])
#          for q in range(len(self.fragments[p])):
#            at_num = self.fragments[p][q]
#            new_list[p].append(at_num)
#            break_list.append(at_num)
#            unfound_list.remove(at_num)
#        while (len(unfound_list) > 0):
#          for p in range(len(new_list)):
#            while (len(new_list[p]) > 0):
#              for q in range(len(new_list[p])-1, -1, -1):
#                at1 = new_list[p][q]
#                for r in range(len(self.bond_tree[at1])):
#                  at2 = self.bond_tree[at1][r]
#                  if (at2 in unfound_list and not at2 in break_list):
#                    self.fragments[p].append(at2)
#                    new_list[p].append(at2)
#                    unfound_list.remove(at2)
#                new_list[p].remove(at1)
#          if (len(unfound_list) > 0):
#            at_new = unfound_list[0]
#            self.fragments.append([at_new])
#            new_list.append([at_new])
#            unfound_list.remove(at_new)
#        for a in range(len(self.fragments)):
#          self.fragments[a] = sorted(self.fragments[a])
#    
#      def print_fsapt_frags(self):
#        fA = open(os.getcwd()+'/'+'fA.dat','w')
#        fB = open(os.getcwd()+'/'+'fB.dat','w') 
#        for frag_num in range(len(self.fragment_names)):
#            if self.fragment_names[frag_num].endswith('_A') or self.fragment_names[frag_num].endswith('_a'): fil = fA
#            elif self.fragment_names[frag_num].endswith('_B') or self.fragment_names[frag_num].endswith('_b'): fil = fB
#            else: 
#                print 'If F-SAPT fragmentation files are desired,\nadd _A and/or _B to corresponding fragments.'
#                sys.exit()
#            fil.write(self.fragment_names[frag_num]+' ')
#            for atom in range(len(self.fragments[frag_num])):
#                fil.write(str(self.fragments[frag_num][atom]+1)+' ')
#            fil.write('\n')
#        fA.close()
#        fB.close()
#    
#      # read in list of intramolecular fragment border atoms from file
#      def get_intrafrags(self, intrafrag_file):
#        if (os.path.isfile(intrafrag_file)):
#          intrafrag_array = get_file_string_array(intrafrag_file)
#          for p in range(len(intrafrag_array)):
#            if (len(intrafrag_array[p]) >= 2):
#              self.fragment_names.append(intrafrag_array[p][0])
#              self.fragments.append([])
#              for q in range(1, len(intrafrag_array[p])):
#                self.fragments[p].append(int(intrafrag_array[p][q])-1)
#    
#      # create bond tree from atomic coordinates
#      def get_bond_tree(self):
#        self.bond_tree = [[] for i in range(self.n_atoms)]
#        self.get_blocks()
#        for block in self.blocks:
#          neighbor_blocks = self.get_neighbor_blocks(block)
#          atom_list = self.get_atoms_from_blocks(neighbor_blocks)
#          for p in range(len(self.blocks[block])):
#            i = self.blocks[block][p]
#            for q in range(len(atom_list)):
#              j = atom_list[q]
#              r2_ij = get_r2ij(self.coords[i], self.coords[j])
#              r2_thresh = bond_thresh * (self.covrad[i] + self.covrad[j])**2
#              if (r2_ij <= r2_thresh and not i == j):
#                if (not j in self.bond_tree[i]):
#                  self.bond_tree[i].append(j)
#                if (not i in self.bond_tree[j]):
#                  self.bond_tree[j].append(i)
#    
#      # parition atoms into spatial blocks
#      def get_blocks(self):
#        self.blocks = {}
#        self.block_res = int(math.ceil(2.0 * bond_thresh * self.maxcovrad))
#        for i in range(self.n_atoms):
#          x, y, z = (int(math.floor(self.coords[i][j])) for j in range(3))
#          xyz_key = get_key(x, y, z, self.block_res)
#          if (not xyz_key in self.blocks):
#            self.blocks[xyz_key] = []
#          self.blocks[xyz_key].append(i)
#    
#      # find occupied blocks which neighbor given block, including self
#      def get_neighbor_blocks(self, block):
#        neighbor_blocks = []
#        x1, y1, z1 = (int(block.split(',')[j]) for j in range(3))
#        br = self.block_res
#        for i in range(3):
#          x = x1 + br*(i-1)
#          for j in range(3):
#            y = y1 + br*(j-1)
#            for k in range(3):
#              z = z1 + br*(k-1)
#              key = get_key(x, y, z, self.block_res)
#              if (key in self.blocks):
#                neighbor_blocks.append(key)
#        return neighbor_blocks
#    
#      # get list of atoms in a set of blocks
#      def get_atoms_from_blocks(self, block_list):
#        atom_list = []
#        for block in block_list:
#          block_atom_list = self.blocks[block]
#          for p in range(len(block_atom_list)):
#            atom_list.append(block_atom_list[p])
#        return atom_list
#    
#      # return covalent radii for all atoms
#      def get_covradii(self):
#        self.maxcovrad = 0.0
#        self.covrad = [0.0 for i in range(self.n_atoms)]
#        for i in range(self.n_atoms):
#          self.covrad[i] = cov_rad[self.at_types[i]]
#          self.maxcovrad = max(self.maxcovrad, self.covrad[i])
#    
#    # check for proper input syntax and assign input variables
#    def get_inputs():
#        input_file, fragment_file = fil_name, fraginfo
#        fsapt = True
#        return input_file, fragment_file, fsapt
#    
#    ## MAIN BLOCK ##
#    
#    input_file, fragment_file, fsapt = get_inputs()
#    geom = Geometry(input_file)
#    geom.get_fragments(fragment_file)
#    geom.print_xyz_allfrags(True)
#    if (fsapt): geom.print_fsapt_frags()
#
#
# 
#### Code calls start here
#try: import pymol
#except: 
#    print_usage()
#    sys.exit()
#
#
##This is for older versions of pymol (< 1.8)
##For some reason it doesn't like drawing bonds with boron.
#B_N_bonds = cmd.find_pairs('n. b', 'n. n', cutoff=2.5)
#for pair in B_N_bonds:
#    cmd.bond('index %s' % pair[0][1], 'index %s' % pair[1][1])
#B_O_bonds = cmd.find_pairs('n. b', 'n. o', cutoff=2.5)
#for pair in B_O_bonds:
#    cmd.bond('index %s' % pair[0][1], 'index %s' % pair[1][1])
#B_C_bonds = cmd.find_pairs('n. b', 'n. c', cutoff=2.0)
#for pair in B_C_bonds:
#    cmd.bond('index %s' % pair[0][1], 'index %s' % pair[1][1])
#
#def fisapt(just_fragment = False, fsapt_path = 'fsapt', fsapt_plot_path = '', input_file = '', template = ''):
#    #remove addition of 'xyz', make read in file name of command line
#    fil_name = cmd.get_names("all")[0]+'.xyz'
#    sanity_check(fil_name)
#    geom = read_original_geometry(fil_name)
#    frags = read_fragment_files()
#    if frags != []:
#        write_prefrag_file(frags)
#        bfs(fil_name, 'fraginfo.dat')
#        #way to doublecheck, coloring should completely eradicate this function
#        #write_frag_geoms()
#        clean_up()
#    cmd.show("sticks", "all")
#    Aatoms, map2A, fA, Batoms, map2B, fB, Catoms, map2orig = get_maps_for_separating(geom)
#    if (just_fragment): 
#        rewrite_fA_and_fB(map2orig, map2A, fA, map2B, fB, Aatoms)
#        color_frags(map2orig, map2A, fA, map2B, fB, Aatoms, Catoms)
#    else:
#        fsapt_path = write_input_file(Aatoms, Batoms, Catoms, geom, fil_name, fsapt_path, fsapt_plot_path, input_file, template)
#        #write_fA_and_fB(map2A, fA, map2B, fB, Aatoms, fsapt_path)
#        try_rewrite_fA_and_fB(map2orig, map2A, fA, map2B, fB, Aatoms,fsapt_path)
#        color_frags(map2orig, map2A, fA, map2B, fB, Aatoms, Catoms)
#        #workflow adaptation
#        #fsapt_path = write_input_file(Aatoms, Batoms, Catoms, geom, fil_name, fsapt_path, fsapt_plot_path, input_file, template)
#        #write_fA_and_fB(map2A, fA, map2B, fB, Aatoms, fsapt_path)
#        #color_frags(map2orig, map2A, fA, map2B, fB, Aatoms, Catoms)
#
#cmd.extend("fisapt",fisapt)
