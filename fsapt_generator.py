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
cov_rad = { 'H' : 0.37,
            'C' : 0.77,
            'O' : 0.73,
            'N' : 0.75,
            'F' : 0.71,
            'P' : 1.10,
            'S' : 1.03,
            'Cl': 0.99,
            'Br': 1.14,
            'I' : 1.33,
            'He': 0.30,
            'Ne': 0.84,
            'Ar': 1.00,
            'Li': 1.02,
            'Be': 0.27,
            'B' : 0.88,
            'Na': 1.02,
            'Mg': 0.72,
            'Al': 1.30,
            'Si': 1.18, 
            'K' : 1.38,
            'Ca': 1.00,
            'Sc': 0.75,
            'Ti': 0.86,
            'V' : 0.79,
            'Cr': 0.73,
            'Mn': 0.67,
            'Fe': 0.61,
            'Co': 0.64,
            'Ni': 0.55,
            'Cu': 0.46,
            'Zn': 0.60,
            'Ga': 1.22,
            'Ge': 1.22,
            'As': 1.22,
            'Se': 1.17,
            'Br': 1.14,
            'Kr': 1.03  }

def bound(A_el, A_coord, B_el, B_coord):
    dist = np.linalg.norm(A_coord - B_coord)
    limit = cov_rad[A_el] + cov_rad[B_el]
    return dist <= 1.1 * limit

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
    
        with open(fil_name, 'w') as fil:
            els = self.mol[0]
            fil.write("%d\n" % (len(els)))
            fil.write("In Angstrom\n")
            coords = self.mol[1]
            for i in range(len(els)): 
                el = els[i]
                c = coords[i]
                fil.write("%2s %7.3f %7.3f %7.3f\n" % (el,c[0],c[1],c[2]))
    
    def write_input(self, fragments, fil_name = "input.dat"):

        template = """memory 60 Gb

molecule mol {
***GEOM***
}

set {
basis         jun-cc-pvdz
df_basis_scf  jun-cc-pvdz-jkfit
df_basis_sapt jun-cc-pvdz-ri

scf_type df
guess sad
freeze_core true

sSAPT0_SCALE TRUE

d_convergence 9

}

energy('fisapt0')"""

        els = self.mol[0]
        coords = self.mol[1]
        A = [key for key in fragments.keys() if '_A' in key]
        A_geom = []
        for a in A:
            inds = fragments[a]
            for ind in inds:
                coord = coords[ind]
                coord = [str(x) for x in coord]
                A_geom.append(els[ind] + ' ' + ' '.join(coord) + '\n')
        B = [key for key in fragments.keys() if '_B' in key]
        B_geom = []
        for b in B:
            inds = fragments[b]
            for ind in inds:
                coord = coords[ind]
                coord = [str(x) for x in coord]
                B_geom.append(els[ind] + ' ' + ' '.join(coord) + '\n')
        C = [key for key in fragments.keys() if '_C' in key]
        C_geom = []
        for c in C:
            inds = fragments[c]
            for ind in inds:
                coord = coords[ind]
                coord = [str(x) for x in coord]
                C_geom.append(els[ind] + ' ' + ' '.join(coord) + '\n')

        geom_to_write = '0 1\n'
        geom_to_write += ''.join(A_geom)
        geom_to_write += '--\n0 1\n'
        geom_to_write += ''.join(B_geom)
        if len(C_geom) > 0: 
            geom_to_write += '--\n0 1\n'
            geom_to_write += ''.join(C_geom)

        input_file = template.replace("***GEOM***", geom_to_write)

        with open(fil_name, 'w') as fil:
            fil.write(input_file)

    def union_find(self):
    
        els = self.mol[0]
        coords = self.mol[1]
        lc = len(coords)
        bonds = []
        #find all bonded atoms
        for i in range(1,lc):
            for j in range(i):
                if not bound(els[i],coords[i],els[j],coords[j]): continue
                bonds.append(set([i,j]))

        #Construct subsets such that each pair has a null 
        #intersection (union find)
        minds = []
        while bonds != []:
            mol = bonds[0]
            bonds.remove(bonds[0])
            intersect = True
            while intersect: 
                intersect = False
                remove = []
                for i in bonds:
                    if i & mol == set([]): 
                        continue
                    for j in i: 
                        mol.add(j)
                    intersect = True
                    remove.append(i)
                for i in remove: 
                    bonds.remove(i)
            minds.append(mol)

        # catch atoms
        for i in range(lc):
            ipresent = False
            for mol in minds:
                if i in mol:
                    ipresent = True
                    break
            if not ipresent:
                minds.append(set([i]))
        
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
        # union find
        # find frag with desired atoms
        # extract and set equal to frag
       
        for frag in frags.keys():
            border_atoms = frags[frag]
            
            other_borders = []
            for frag2 in frags.keys():
                if frag == frag2: 
                    continue
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
            new_copy.union_find()
            for border in border_atoms:
                target_coords = self.mol[1][border]
                found_frag = False
                for new_frag in new_copy.frags:
                    for coords in new_frag.mol[1]:
                        if not np.array_equal(target_coords, coords): 
                            continue
                        found_frag = new_frag
                        break
                    if found_frag != False: break
                if found_frag != False: break
            for atom_coords in found_frag.mol[1]:
                 for ind in range(len(self.mol[1])):
                    if not np.array_equal(self.mol[1][ind], atom_coords):
                        continue
                    if ind not in frags[frag]: 
                        frags[frag].append(ind)
                    break         

        # Clean and combine all C
        copy_mol = self.copy() 
            
        for frag in frags.keys():
            classification = frag.split('_')[-1].upper()
            if classification != 'C':
                for atom in frags[frag]:
                    copy_mol.mol[0][atom] = 'XXX'
                continue
            del frags[frag]
    
        frags['ISAPT_C'] = []
        for ind in range(len(self.mol[0])):
            atom = copy_mol.mol[0][ind]
            if atom == 'XXX': 
                continue
            frags['ISAPT_C'].append(ind)
        if frags['ISAPT_C'] == []:
            del frags['ISAPT_C']    
      
        self.color_frags(frags)
        self.write_frags(frags)
       
        return frags

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
                

    #color fragments based on fragmentation from union find 
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
        
        # All C fragments are in "ISAPT_C"
        if classification == 'C':
            cmd.delete(name)

    # Fill up with fragments and color the fragments 
    fragments = total_molecule.cut(frags)

    # Should write input, fA, and fB outside of class, need info from above
    total_molecule.write_input(fragments, fil_name = cmd.get_names("all")[0]+".in")

cmd.extend("fisapt",fisapt)
cmd.load("/Users/brandonbakr/Desktop/demo/hexaphenylbenzene.xyz")
cmd.show("sticks", "all")
cmd.hide("spheres", "all")

#diphenyl = pm_mol('examples/diphenyl_benzene.xyz')
#diphenyl.cut(A = [12], B = [24], C = [2, 4])


