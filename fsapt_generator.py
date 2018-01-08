import os
import numpy as np
import copy
from pymol import cmd 
import random
import sys

def fsapt_generator():
    pass 
    # Should take border atoms (only!) from selections for monomers A & B (and C if C)
    # Should return fA and fB 

# covalent (or ionic) radii by atomic element in angstroms from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}

def bound(A_el, A_coord, B_el, B_coord):
    dist = np.linalg.norm(A_coord - B_coord)
    limit = cov_rad[A_el] + cov_rad[B_el]
    return dist <= limit

class pm_mol():

    def __init__(self, mol): 
    
        if type(mol) is str and mol.endswith('.xyz'):
            self.fil = mol
            with open(mol,'r') as ofil:
                n = int(next(ofil).split('\n')[0])
                next(ofil)
                els = []
                geom = np.array([]).reshape(0,3)
                for l in range(n):
                    lin = next(ofil)
                    ls = lin.split()
                    els.append(ls[0])
                    coords = [float(x) for x in ls[1:]]
                    geom = np.vstack([geom,coords])
                mol = [els,geom]
        self.mol = mol
        self.frags = []

    def print_out(self):
    
        els = self.mol[0]
        print(len(els))
        print("In Angstrom")
        coords = self.mol[1]
        for i in range(len(els)): 
            el = els[i]
            c = coords[i]
            print("%2s %7.3f %7.3f %7.3f" % (el,c[0],c[1],c[2]))
        print("")
    
    def copy(self):
        return copy.deepcopy(self)

    def write_out(self,fil_name):
    
        with open(fil_name,'w') as fil:
            els = self.mol[0]
            fil.write("%d\n" % (len(els)))
            fil.write("In Angstrom\n")
            coords = self.mol[1]
            for i in range(len(els)): 
                el = els[i]
                c = coords[i]
                fil.write("%2s %7.3f %7.3f %7.3f\n" % (el,c[0],c[1],c[2]))

    def bfs(self):
    
        els = self.mol[0]
        coords = self.mol[1]
        lc = len(coords)
        bonds = []
        #find all bonded atoms
        for i in range(1,lc):
            for j in range(i):
                if not bound(els[i],coords[i],els[j],coords[j]): continue
                bonds.append(set([i,j]))
        #Construct subsets such that each pair has a null intersection (BFS)
        minds = []
        while bonds != []:
            mol = bonds[0]
            bonds.remove(bonds[0])
            intersect = True
            while intersect: 
                intersect = False
                remove = []
                for i in bonds:
                    if i & mol == set([]): continue
                    for j in i: mol.add(j)
                    intersect = True
                    remove.append(i)
                for i in remove: bonds.remove(i)
            minds.append(mol)
        #build mol out of separated inds
        mols = []
        for mol in minds:
            mels = []
            mcoords = []
            for i in mol:
                mels.append(els[i])
                mcoords.append(coords[i])
            mols.append(pm_mol([mels,mcoords]))
        self.frags = mols

    def cut(self, frags):
        # remove atoms of other frags
        # bfs
        # find frag with desired atoms
        # extract and set equal to frag
        
        for frag in frags.keys():
            border_atoms = frags[frag]
            
            other_borders = []
            for frag2 in frags.keys():
                if frag == frag2: continue
                for border_atom in frags[frag2]: 
                    other_borders.append(border_atom)           
 
            copy_mol = self.copy() 
            
            for border in other_borders:
                copy_mol.mol[0][border] = 'XXX'
            
            new_copy_els = []
            new_copy_coords = []
            for i in range(len(copy_mol.mol[0])):
                if copy_mol.mol[0][i] == 'XXX': continue
                new_copy_els.append(copy_mol.mol[0][i]) 
                new_copy_coords.append(copy_mol.mol[1][i]) 
            
            new_copy = pm_mol([new_copy_els, new_copy_coords])
            new_copy.bfs()
            for border in border_atoms:
                target_coords = self.mol[1][border]
                found_frag = False
                for new_frag in new_copy.frags:
                    for coords in new_frag.mol[1]:
                        if not np.array_equal(target_coords, coords): continue
                        found_frag = new_frag
                        break
                    if found_frag != False: break
                if found_frag != False: break
            for atom_coords in found_frag.mol[1]:
                 for ind in range(len(self.mol[1])):
                    if not np.array_equal(self.mol[1][ind], atom_coords): continue
                    if ind not in frags[frag]: frags[frag].append(ind)
                    break         
       
        self.color_frags(frags)
        self.write_frags(frags)


    # Need to resolve naming of `frags`, already exists as deprecated member data from crystal class
    def write_frags(self, frags):

        with open('fA.dat', 'w') as fA:
            for frag in frags:
                classification = frag.split('_')[-1].upper()
                if classification != 'A': continue
                fA.write(frag+' ')
                for atom in frags[frag]: 
                    fA.write(str(atom+1)+' ')
                fA.write('\n')

        with open('fB.dat', 'w') as fB:
            for frag in frags:
                classification = frag.split('_')[-1].upper()
                if classification != 'B': continue
                fB.write(frag+' ')
                for atom in frags[frag]: 
                    fB.write(str(atom+1)+' ')
                fB.write('\n')
                

    #color fragments based on fragmentation from bfs
    def color_frags(self, frags):

        for frag in frags:
            selection = ''
            for atom in frags[frag]:
                selection += 'rank '+str(atom)+' '
            cmd.select("("+frag+")", selection)
            g = random.uniform(0, 1)           
            classification = frag.split('_')[-1].upper()
            if classification == 'C':
                cmd.color("grey", "("+frag+")")
                continue
            elif classification == 'A':
                r = random.uniform(0.5, 1)           
                b = random.uniform(0, 0.5)           
            elif classification == 'B':
                r = random.uniform(0, 0.5)           
                b = random.uniform(0.5, 1) 
            color = [r, g, b]
            cmd.set_color("color_"+frag, color) 
            cmd.color("color_"+frag, "("+frag+")")
    

# Read main geometry
def read_original_geometry():

    fil_name = cmd.get_names("all")[0]+'.xyz'
    geometry = pm_mol(fil_name)
    return geometry


def fisapt():
    
    # Make it pretty 
    cmd.show("sticks", "all")

    # Initialize pm molecule object
    total_molecule = read_original_geometry()

    allowed_classifications = ['A', 'B', 'C']
    # Take in user border atoms
    frag_names = cmd.get_names("all")[1:]
    frags = {}
    for name in frag_names:
        classification = name.split('_')[-1].upper()
        if classification not in allowed_classifications: 
            print("USAGE!")
            sys.exit()
        frags[name] = []
        stored.list=[]
        cmd.iterate("("+name+")","stored.list.append((name,rank))")
        for atom in stored.list: 
            frags[name].append(atom[1])

    # Fill up with fragments and color the fragments 
    fragments = total_molecule.cut(frags)

    # Should write input, fA, and fB outside of class, need info from above

cmd.extend("fisapt",fisapt)


#diphenyl = pm_mol('examples/diphenyl_benzene.xyz')
#diphenyl.cut(A = [12], B = [24], C = [2, 4])


#### OLD FRAG GEN ###
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
#
#
#os.mkdir('fsapt/')
#with open('fsapt/fA.dat','w') as fil:
#    fil.close()    
#with open('fsapt/fB.dat','w') as fil:
#    fil.close()    
#
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
#
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

