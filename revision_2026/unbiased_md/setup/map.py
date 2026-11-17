import collections, pickle, numpy as np, parmed as pmd
from rdkit import Chem
from rdkit.Chem import rdFMCS
T="/media/rahul/EXTERNAL_USB/binding_affinity/TI/ti_rep1/bundle/ti_rep1/topologies/"
C="/home/rahul/Downloads/CPD_DNA_DAMAGE_data_files/NAR_paper_data_files/INITIAL_STRUCTURE/3_cpd_pb.wb"
pb=pmd.load_file(T+"system.prmtop",T+"system.inpcrd")
dk=pmd.load_file(T+"system_docked_cleaned.prmtop",T+"system_docked_cleaned.inpcrd")
cp=pmd.load_file(C+".prmtop",C+".rst7")
for n,s in (("pb",pb),("docked",dk)):
    c=collections.Counter(r.name for r in s.residues); print(n,"atoms",len(s.atoms),"res",len(s.residues),{k:v for k,v in c.items() if k in("WAT","Na+","T64")},"solute",sum(len(r.atoms) for r in s.residues if r.name not in("WAT","Na+")),"box",s.box[:3] if s.box is not None else None, "res602:",s.residues[601].name,len(s.residues[601].atoms))
def mol(res):
    m=Chem.RWMol(); idx={}
    for a in res.atoms:
        el=a.element_name if hasattr(a,'element_name') else 'C'
        at=Chem.Atom(el); idx[a.idx]=m.AddAtom(at)
    for a in res.atoms:
        for p in a.bond_partners:
            if p.idx in idx and a.idx<p.idx: m.AddBond(idx[a.idx],idx[p.idx],Chem.BondType.SINGLE)
    m=m.GetMol(); Chem.SanitizeMol(m,Chem.SANITIZE_NONE); m.UpdatePropertyCache(strict=False); Chem.FastFindRings(m); return m
rt=pb.residues[601]; rc=cp.residues[601]
mt,mc=mol(rt),mol(rc)
mh=lambda m:Chem.RemoveHs(m,sanitize=False)
res=rdFMCS.FindMCS([mt,mc],atomCompare=rdFMCS.AtomCompare.CompareElements,bondCompare=rdFMCS.BondCompare.CompareAny,ringMatchesRingOnly=True,completeRingsOnly=False,timeout=60)
print("MCS atoms",res.numAtoms,"bonds",res.numBonds)
q=Chem.MolFromSmarts(res.smartsString); 
mt_m=mt.GetSubstructMatch(q); mc_m=mc.GetSubstructMatch(q)
pairs=[(rt.atoms[i].name,rc.atoms[j].name) for i,j in zip(mt_m,mc_m)]
print(len(pairs),pairs)
print("T64 unmatched:",[a.name for i,a in enumerate(rt.atoms) if i not in mt_m]); print("MOL unmatched:",[a.name for j,a in enumerate(rc.atoms) if j not in mc_m])
pickle.dump(pairs,open("pairs.pkl","wb"))
