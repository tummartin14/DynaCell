#! /usr/bin/python
###################################################
##                                               ##
##   Author :: Antoine MARION                    ##
##     Date :: 2017.05.03                        ##
## Comments :: Converts Amber topology and       ##
##             coordinate files to Dynacell      ##
##             input format.                     ##
##             This script uses ParmEd classes   ##
##             written by J. Swail and available ##
##             in AmberTools15-16                ##
##                                               ##
###################################################
# Edit August 2020 by Martin Zachmann: works now
# without the need to provide a ligand (apo), and
# pdb-output for atom-names is now correctly aligned

import os
import sys
import math
import re
import numpy as np

# update this list if you have new solvent residue
Known_solvent_residue = [
'WAT',
'HOH',
'DMS'
]

List_of_chains_prot = [
'A',
'B',
'C',
'D',
'E',
'F',
'G',
'H',
'I',
'J',
'K',
'M',
'N',
'O',
'P',
'Q',
'R',
'T',
'U',
'V',
'W',
'X',
'Y',
'Z'
]

ligand_chain = 'L'

solvent_chain = 'S'

Cal2Joule = 4.184
# For vdW, Amber uses the "radius" of the atoms.
# For the two identical atoms Rmij = 2*radius_i
# Rmij = 2**(1/6) * sigma_ij
# thus radius_i = 1/2  2**(1/6) * sigma_i
# thus sigma_i = 2*(5/6) * radius_i
Radius2Sigma = (2.**(5./6.)) * 0.1
rad2deg = 180./math.pi

## Functions

def PrintDefault(file):
    file.write('[ defaults ]\n')
    file.write(";combination-rule  fudgeLJ  fudgeQQ  nrExcls\n")
    # combination-rule=2 (LJ params treated as sigma, epsilon and 
    #                    combined as sigma=1/2*(sigma1+sigma2)
    #                    epsilon=sqrt(epsilon1*epsilon2))
    comb_rule = 2
    # fudgeLJ and fudgeQQ are the factors by which LJ and coulomb
    # 14-interactions will be scaled in the pairs section if 
    # these params are not set specifically
    fudgeLJ = 1.0/2.0
    fudgeQQ = 1.0/1.2 
    # nrExcls is the number of neighbors up to which all nb 
    # interactions are automatically excluded amber very likely
    # uses 3 as all other forcefields, however you can set any number here
    nrExcls = 3
    file.write('%17d%20.8f%20.8f%9d\n' % (comb_rule,fudgeLJ,fudgeQQ,nrExcls))

def PrintAtoms(parm,file):
    atnum = 0
    file.write('[ atoms ]\n')

    lj_funct_form = 4

    file.write("\
;     nr    type    name    chargeGrNr          charge            mass           sigma         epsilon    atomicNr\n")

    for at in parm.atoms:
        atnum += 1
        at.Atnum = atnum # need to set it explicitely to use it later
        attype = at.type
        atname = at.name
        chg_grp = atnum # not used so set equal to atom number
        atcharge = at.charge
        atmass = at.mass 
        lj_sigma = at.rmin * Radius2Sigma
        lj_epsilon = at.epsilon * Cal2Joule
        atomic_number = at.atomic_number

        file.write("%8d%8s%8s%14d%16.6e%16.6e%16.6e%16.6e%12d\n" % (
        atnum,
        attype,
        atname,
        chg_grp,
        atcharge,
        atmass,
        lj_sigma,
        lj_epsilon,
        atomic_number
        ))

def PrintBonds(parm,file):
    file.write('[ bonds ]\n')
    
    file.write(";  ati_idx   atj_idx   bond_func              b0             kb0    rotatable? (R)\n")

    for bond in parm.bonds:
        ati = bond.atom1.Atnum
        atj = bond.atom2.Atnum
        bond_funct = 1
        req = bond.type.req * 0.1
        k = bond.type.k * 2. * 100. * Cal2Joule

        if bond.Flexible:
            file.write("%10d%10d%12d%16.4e%16.6e%18s\n" % (
            ati,
            atj,
            bond_funct,
            req,
            k,
            'R'
            ))
        else:
            file.write("%10d%10d%12d%16.4e%16.6e\n" % (
            ati,
            atj,
            bond_funct,
            req,
            k
            ))

def PrintPairs(parm,file):
    ## 1-4 pairs are defined from the dihedral partners of each atoms.
    ## The list of dihedral partners contains only the pairs that are NOT
    ## in the list of bond and angle partners. Thus only the 1-4 partners.
    file.write('[ pairs ]\n')

    file.write(";  ati_idx   atj_idx\n")
    # Use the following line if explicit setting of LJ parameters for each pair
    #file.write("; ati_idx   atj_idx           sigma         epsilon\n")

    pair_funct = 4

    already_printed_pairs = []
    for atomi in parm.atoms:
        ati = atomi.Atnum
        for atomj in atomi.dihedral_partners:
            atj = atomj.Atnum
            if (not (ati,atj) in already_printed_pairs) and (not (atj,ati) in already_printed_pairs):
                already_printed_pairs.append((ati,atj))
                file.write("%10d%10d\n" % (
                ati,
                atj
                ))
                # the following can be used to set explicitely the LJ parameters for a given pair
                # need to check before using that combinatin rule and scaling is correct
                #dih = FindDihedral(atomi,atomj)
                #scnb = 1. / dih.type.scnb
                #sigma = (atomi.rmin + atomj.rmin) * 0.5 * Radius2Sigma * scnb
                #epsilon = math.sqrt(atomi.epsilon*atomj.epsilon) * Cal2Joule * scnb

                #file.write("%9d%10d%16.6e%16.6e\n" % (
                #ati,
                #atj,
                #sigma,
                #epsilon
                #))


def FindDihedral(atom1,atom2):
    # Finds which dihedral contains the two atoms
    for dih1 in atom1.dihedrals:
        for dih2 in atom2.dihedrals:
            if dih1 == dih2:
                return dih1
    # If we arrive here, it means that we didn't find two equivalent dihedrals
    sys.stderr.write("Error in pairs definition...\n")
    exit(2)

def PrintAngles(parm,file):
    file.write('[ angles ]\n')

    file.write(";  ati_idx   atj_idx   atk_idx   func-form           theta          ktheta\n")
    
    angle_funct = 1

    for ang in parm.angles:
        ati = ang.atom1.Atnum
        atj = ang.atom2.Atnum
        atk = ang.atom3.Atnum
        theta_eq = ang.type.theteq
        k = ang.type.k * 2. * Cal2Joule
         
        file.write("%10d%10d%10d%12d%16.6e%16.6e\n" % (
        ati,
        atj,
        atk,
        angle_funct,
        theta_eq,
        k
        ))

def GetDih(parm,old_amber=False):
    if old_amber:
        fact = rad2deg 
    else:
        fact = 1.0

    all_dih = {}
    all_imp = {}
    for dih in parm.dihedrals:
        ati = dih.atom1.Atnum
        atj = dih.atom2.Atnum
        atk = dih.atom3.Atnum
        atl = dih.atom4.Atnum

        phase = dih.type.phase * fact
        k = dih.type.phi_k * Cal2Joule
        mult = dih.type.per

        #if dih.signs[0] < 0 and dih.signs[1] < 0: ## improper
        if dih.improper: ## improper
            if (ati,atj,atk,atl) in all_imp.keys():
                all_imp[(ati,atj,atk,atl)].append((phase,k,mult))
            elif (atl,atk,atj,ati) in all_imp.keys():
                all_imp[(atl,atk,atj,ati)].append((phase,k,mult))
            else:
                all_imp[(ati,atj,atk,atl)] = [(phase,k,mult)]
        else: # dihedral
            if (ati,atj,atk,atl) in all_dih.keys():
                all_dih[(ati,atj,atk,atl)].append((phase,k,mult))
            elif (atl,atk,atj,ati) in all_dih.keys():
                all_dih[(atl,atk,atj,ati)].append((phase,k,mult))
            else:
                all_dih[(ati,atj,atk,atl)] = [(phase,k,mult)]
    return all_dih, all_imp

def PrintDih(all_dih,file,section='dihedrals'):
    func_form = 1
    file.write('[ %s ]\n' % section)
    file.write("; ati_idx  atj_idx  atk_idx  atl_idx")
    if section == 'impropers': file.write("  func_form")
    file.write("           phase            kphi  mult (last 3 repeated if multiple)\n")
    for key in sorted(all_dih.keys()):
        ati, atj, atk, atl = key
        line = "%9d%9d%9d%9d" % (
        ati,
        atj,
        atk,
        atl
        )
        if section == 'impropers': line += "%11d" % func_form
        for phase,k,mult in all_dih[key]:
            line += "%16.4e%16.6e%6d" % (phase,k,mult)
        line += '\n'
        file.write(line)

def ParseInputFile(parm,input_file):
    ligand, flexible = [], []
    m = None
    temp_flex = []

    if input_file != None:
        try:
            f = open(input_file,'r')
        except:
            sys.stderr.write("Error while opening input file %s\n" % input_file)
            exit(2)

        line = f.readline()
        while line != '':
            if re.search("&lig",line) != None:
                m = line.split('=')[1].strip()
            elif re.search("&flex",line):
                line = f.readline()
                while re.search("&endflex",line) == None:
                    if not line.startswith('#'):
                        temp = line.split('#')[0].split(',') # remove comment if there is
                        if len(temp) != 2:
                            sys.stderr.write("Error while reading flexible bonds in %s\n" % input_file)
                            sys.stderr.write("Please verify the following line in your input file :\n%s" % line)
                            exit(2)
                        else:
                            temp_flex.append((temp[0].strip(),temp[1].strip()))
                    else: # skip lines starting with "#"
                        pass
                    line = f.readline()
            else:
                pass
            line = f.readline()

    if m == '' or m == None:
        sys.stdout.write("Note: No ligand was defined\n")
        # return value unused
        ligand = SetAtomChain(parm)
    else:
        try:
            mask = AmberMask(parm,m)
            sel = mask.Selection()
            ligand = SetAtomChain(parm,sel)
            if len(ligand) == 0:
                sys.stderr.write("Error: The ligand mask did not match any atoms\n")
                exit(2)
        except:
            sys.stderr.write("Error while parsing ambermask %s\n" % m)
            exit(2)

    if len(temp_flex) != 0: flexible = SetFlexibleBonds(parm,temp_flex)
    if len(flexible) == 0:
        sys.stdout.write("Note: No flexible bonds were defined\n")

    return

def SetFlexibleBonds(parm,flex):
    # Check if each flexible bond is part of the ligand and set it as flexible
    all_flex_bonds = []
    for m1, m2 in flex:
        at1 = GetAtom(parm,m1)
        at2 = GetAtom(parm,m2)
        if at1.Chain != 'L' or at2.Chain != 'L':
            sys.stderr.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            sys.stderr.write("!!! WARNING !!! : One or both of the atom selections   !!!\n")
            sys.stderr.write("!!! %s,%s defining a flexible bond are not           !!!\n" % (m1,m2))
            sys.stderr.write("!!! part of the ligand.                                !!!\n")
            sys.stderr.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

        bond = FindBond(parm,at1,at2)
        if bond == None:
            sys.stderr.write("Error while setting flexible bonds.\n")
            sys.stderr.write("The selection %s,%s does not correspond to any bond.\n" % (m1,m2))
            exit(2)
        else:
            bond.Flexible = True
            all_flex_bonds.append(bond)
    return all_flex_bonds

def FindBond(parm,a1,a2):
    # Find the bond defined by two atoms
    this_bond = None
    #for b in parm.bonds_inc_h + parm.bonds_without_h: # Only in AmberTools14
    for b in parm.bonds:
        if (a1 in b) and (a2 in b):
            this_bond = b
            break
    return this_bond

def GetAtom(parm,m):
    # Get the atom that corresponds to a atom mask
    mask = AmberMask(parm,m)
    sel = np.array(mask.Selection())
    idx = np.where(sel == 1)[0]
    if len(idx) == 0:
        sys.stderr.write("Error while parsing mask %s for flexible bond definition.\n" % m)
        sys.stderr.write("The mask did not return any atom.\n")
        exit(2)
    elif len(idx) > 1:
        sys.stderr.write("Error while parsing mask %s for flexible bond definition.\n" % m)
        sys.stderr.write("The mask returned more than one atom.\n")
        exit(2)
    else:
        #at = parm.atom_list[idx[0]] # Only in AmberTools14
        at = parm.atoms[idx[0]]
    return at

def IsTer(res):
    ister = False
    all_res = []
    for atom in res.atoms:
        for bond in atom.bonds:
            all_res += [bond.atom1.residue.idx,bond.atom2.residue.idx]
    if not res.idx+1 in all_res:
        ister = True
    return ister

def SetAtomChain(parm,mask=[]):
    # Assing a chain label to each atom according to the defined mask
    lig = []

    i = -1
    Nchain_prot = 0
    prot = False
    for res in parm.residues:
        res.ter = IsTer(res)
        for atom in res.atoms:
            i += 1
            if len(mask) > 0 and mask[i] == 1:
                atom.Chain = ligand_chain 
                lig.append(atom)
                prot = False
            elif res.name in Known_solvent_residue:
                atom.Chain = solvent_chain
                prot = False
            else:
                try:
                    atom.Chain = List_of_chains_prot[Nchain_prot] 
                except:
                    Nchain_prot = 0
                    atom.Chain = List_of_chains_prot[Nchain_prot] 
                prot = True
        if res.ter and prot: Nchain_prot += 1

    return lig

def PrintCoordinates(parm,file,pf,cf):
    file.write("HEADER  Dynacell coordinate file generated using Amber2Dynacell from %s and %s\n" % (pf,cf))
    file.write("MODEL     1\n")
    for atom in parm.atoms:
        line = MakePDBline(atom)
        file.write("%s\n" % line)
    file.write("ENDMDL\n")
    return

def MakePDBline(atom):
    pdbline = ''
    record = "ATOM"
    iat = atom.Atnum
    #print dir(atom)
    #print  dir(atom.residue)
    #print  atom.residue.ter
    #exit(0)
    
    atname = atom.name
    locid = ""
    resname = atom.residue.name
    chain = atom.Chain
    ires = atom.residue.idx + 1 # seems like it starts at 0
    codeinsres = ""
    # Coordinates are in Ang
    x = atom.xx
    y = atom.xy
    z = atom.xz
    occ = 0.0
    #temp = 0.0
    #segid = ""
    #elem = ""
    #chg = ""
    isLongElem = True
    an = atom.atomic_number
    if (an == 1) or (an >=5 and an <= 9) or (an == 15) or (an == 16) or (an == 19) or (an == 23) or (an == 39) or (an == 53) or (an == 74) or (an == 92):
        isLongElem = False

    if (len(atname) >= 4) or (len(atname) == 0) or (isLongElem): # just to align the atom names properly
        #fmt = "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s"
        fmt = "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f"
    elif atname[0].isalpha():
        #fmt = "%-6s%5d  %-3s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s"
        fmt = "%-6s%5d  %-3s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f"
    else:
        #fmt = "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s"
        fmt = "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f"

    pdbline = fmt % (\
    record,
    iat,
    atname,
    locid,
    resname,
    chain,
    ires,
    codeinsres,
    x,
    y,
    z,
    occ,
    #temp,
    #segid,
    #elem,
    #chg
    )

    return pdbline

#############################


## AmberTools16-17
try:
    from parmed.amber.readparm import AmberParm
    from parmed.amber.mask import AmberMask
    from parmed.topologyobjects import Atom
    from parmed.topologyobjects import Bond
    oldAmber = False
except:
## AmberTools15
    try:
        from chemistry.amber.readparm import AmberParm
        from chemistry.amber.mask import AmberMask
        from chemistry.topologyobjects import Atom
        from chemistry.topologyobjects import Bond
        oldAmber = True
    except:
        sys.stderr.write("Problem: couldn't load parmed packages.\n")
        sys.stderr.write("Please check the version of AmberTools that you are using.\n")
        sys.stderr.write("This version only works with AmberTools15-16.\n")
        sys.stderr.write("Make sure that you sourced Amber's configuration file $AMBERHOME/amber.sh\n")
        exit(1)

from argparse import ArgumentParser

## Add atributes to Atom and Bond classes
Atom.Chain = ''
Atom.Atnum = None
Bond.Flexible = False

## Parse command line options
parser = ArgumentParser()
parser.add_argument('-O','--Overwrite',dest='overwrite',default=False,
                    help='Allow Amber2Dynacell to overwrite existing files', action='store_true')
parser.add_argument('-p','--prmtop',dest='prmtop_file',default=None, required=True,
                    help='Amber format topology file to be converted to Dynacell format', action='store')
parser.add_argument('-c','--coordinate',dest='crd_file',default=None,
                    help='Amber format restart coordinate file to be converted to Dynacell pdb format. Equivalent to .inpcrd', action='store')
parser.add_argument('-i','--input',dest='input_file',default=None,
                    help='Input file containing the selection of atoms that define the ligand in ambermask format (optional, \"&lig = myMask\") and a list of flexible bonds (optional, between \"&flex\" and \"&endflex\" in the form \"myMask1, myMask2\" each on a new line)', action='store')
parser.add_argument('-o','--output',dest='output_suffix',default=None, required=True,
                    help='Output suffix to be used for naming the Dynacell .top and .pdb files', action='store')
options = parser.parse_args() 

## Setup

out_top_name = options.output_suffix + ".top"
if os.path.exists(out_top_name) and not options.overwrite:
    sys.stderr.write("Error : file %s exists.\n" % out_top_name)
    sys.stderr.write("        Use the \'-O\' option to overwrite files.\n")
    exit(1)
else:
    out_top_file = open(out_top_name,'w')

if options.crd_file != None:
    out_crd_name = options.output_suffix + ".pdb"
    if os.path.exists(out_crd_name) and not options.overwrite:
        sys.stderr.write("Error : file %s exists.\n" % out_crd_name)
        sys.stderr.write("        Use the \'-O\' option to overwrite files.\n")
        exit(1)
    else:
        out_crd_file = open(out_crd_name,'w')
        do_crd = True
else:
    sys.stdout.write("Note: No coordinate file will be written, please make sure\n")
    sys.stdout.write("      that the coordinate file you will be using matches the\n")
    sys.stdout.write("      topology written here, and contains correct chain-IDs\n")
    do_crd = False

if do_crd:
    prmtop = AmberParm(options.prmtop_file,options.crd_file)
else:
    prmtop = AmberParm(options.prmtop_file)

ParseInputFile(prmtop,options.input_file)


out_top_file.write('; AMBER\n')
out_top_file.write('; Dynacell topology file generated using Amber2Dynacell')
out_top_file.write('; from %s and %s\n' % (options.prmtop_file,options.crd_file))
PrintDefault(out_top_file)
out_top_file.write('\n')
PrintAtoms(prmtop,out_top_file)
out_top_file.write('\n')
PrintBonds(prmtop,out_top_file)
out_top_file.write('\n')
PrintPairs(prmtop,out_top_file)
out_top_file.write('\n')
PrintAngles(prmtop,out_top_file)
out_top_file.write('\n')

dihedrals, impropers = GetDih(prmtop,oldAmber)
PrintDih(dihedrals,out_top_file,section='dihedrals')
out_top_file.write('\n')
PrintDih(impropers,out_top_file,section='impropers')

if do_crd: PrintCoordinates(prmtop,out_crd_file,options.prmtop_file,options.crd_file)

## End
out_top_file.close()
if do_crd: out_crd_file.close()


