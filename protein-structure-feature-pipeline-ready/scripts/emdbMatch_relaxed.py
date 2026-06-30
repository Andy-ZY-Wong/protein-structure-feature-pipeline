#!/usr/bin/env python3
import argparse, sys, collections, mmcif.io.PdbxReader, Bio.SeqIO, numpy as np
# from IPython import embed  # removed: not used

three2one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'} 

if __name__ != '__main__': sys.exit()
parser = argparse.ArgumentParser(description="match CryoNet.Fold model to input sequence, in-situ template chains and PDB chains")
parser.add_argument('-i', '--input', required=True, help = 'path to predicted .cryofold file')
parser.add_argument('-s', '--seq', required=True, help = 'path to input .manAlign file, with coupled .faa file')
parser.add_argument('-t', '--templ', help = 'path to template .cif file')
parser.add_argument('-o', '--output', help = 'path to output .autoDel file')

args = parser.parse_args()
if args.output is None: args.output = args.input.replace('.cryofold', '.autoDel')

# Gen.py already matched .faa to .cif; .manAlign has description and gapped
tnamesD, t2c, c2gap, descD= collections.defaultdict(list), {}, {}, {}
for x in Bio.SeqIO.parse(args.seq.replace('manAlign', 'faa'), 'fasta'): 
    cname, tname = x.id.rsplit('_', 1)
    t2c.update({tname: cname}), tnamesD[ str(x.seq) ].append( tname ) 
for x in Bio.SeqIO.parse(args.seq, 'fasta'): 
    if x.id.count('_')!=2: # assume protein in front of chains
        gquery, seq = str(x.seq), str(x.seq).replace('X', 'G').replace('-', '')
        descD[seq] = x.description.split(' ', 3)[3]
    if x.id.count('_')==2: 
        c2gap[x.id] = ''.join(y for x,y in zip(gquery, str(x.seq)) if x!='-')  
t2gap = {t: c2gap[c] for t,c in t2c.items()}  # template name to gapped sequence

# CryoNet.Fold predicts full-length
predD, pnamesD, lsD = {}, collections.defaultdict(list), collections.defaultdict(list)
for l in open(args.input).readlines(): 
    if l.startswith('ATOM') and ' CA ' in l: lsD[l[20:22].strip()].append(l)
for n, ls in lsD.items(): 
    seq = ''.join([ three2one[l[17:20]] for l in ls])
    coordA = np.array([(l[30:38], l[38:46], l[46:54]) for l in ls]).astype(float)
    pnamesD[seq].append(n), predD.update({n: coordA})

atts = 'auth_asym_id label_seq_id label_alt_id label_comp_id label_atom_id Cartn_x Cartn_y Cartn_z'
if args.templ is None: templD = {}
else:
    cif, templD, attL = [], collections.defaultdict(list), atts.split()
    mmcif.io.PdbxReader.PdbxReader(open(args.templ)).read(cif)
    atom_site, pos2coord = cif[0].getObj('atom_site'), collections.defaultdict(list)
    for i in range(len(atom_site)): 
        auth, pos, alt, res, atom, x, y, z = [atom_site.getValue(k, i) for k in attL]
        if atom_site.getValue('pdbx_PDB_model_num', i)!='1' or pos=='.': continue
        pos2coord[ (auth, pos) ].append( (atom, x, y, z) ) 
    for (auth, pos), axyzL in pos2coord.items(): 
        xyz = np.array([v[1:] for v in axyzL]).astype(float).mean(axis=0)
        for v in axyzL: xyz = np.array(v[1:]).astype(float) if v[0]=='CA' else xyz 
        templD[auth].append(xyz) # use CA if available, otherwise center of mass
    templD = {k: np.array(v) for k,v in templD.items()}

# find nearset templ chain, and transfer gaps to pred chain
def get_diffA(predA, templA, gap): return predA[ np.array(list(gap))!='-' ]-templA
def get_rmsd(*args): return np.sqrt(np.mean( get_diffA(*args)**2 ))*np.sqrt(3)  
n2seq = {n:seq for seq, ns in pnamesD.items() for n in ns}
rcsbL, posLD = [], collections.defaultdict(list)
for n, predA in predD.items(): 
    tS = set(tnamesD[ n2seq[n] ]) & templD.keys()
    rmsdD = { get_rmsd(predA, templD[t], t2gap[t] ):t for t in tS }
    if rmsdD=={}: continue
    t = rmsdD[ np.min(list(rmsdD)) ]
    for i,x in enumerate(t2gap[t]): posLD[n].append(i+1) if x=='-' else None
    rcsbL.append(f">{t2c[t]}_{n}\n{n2seq[n]}\n")
    templD.pop(t)

# write FASTA sequence file and deletion file
open(args.seq.replace('manAlign', 'toRCSB'), 'w').write( ''.join(rcsbL) ) 
fastaL, delL = [], []
for i, (seq, pnames) in enumerate(pnamesD.items()): 
    desc = descD.get(seq, "unknown_sequence_not_in_manAlign")
    fastaL.append(f">{i+1}|Chains {', '.join(pnames)}|{desc}\n{seq}\n")
    for n in pnames: 
        oneL = ''.join(('-' if i+1 in posLD[n] else c) for i,c in enumerate(seq))
        delL.append(f">{i+1}_{n}\n{seq}\n>{i+1}_{n}\n{oneL}\n")
open(args.seq.replace('manAlign', 'fasta'), 'w').write( ''.join(fastaL) )  
open(args.output, 'w').write( ''.join(delL) ) 

