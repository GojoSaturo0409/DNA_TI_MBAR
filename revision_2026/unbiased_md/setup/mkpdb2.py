import pickle, numpy as np, parmed as pmd
T="/media/rahul/EXTERNAL_USB/binding_affinity/TI/ti_rep1/bundle/ti_rep1/topologies/"
C="/home/rahul/Downloads/CPD_DNA_DAMAGE_data_files/NAR_paper_data_files/INITIAL_STRUCTURE/3_cpd_pb.wb"
pairs=dict((m,t) for t,m in pickle.load(open("pairs.pkl","rb")))      # MOL atom -> T64 atom
pairs.update({"O14":"O8","H20":"H12"})                                 # the one unmatched hydroxyl O/H pair
cp=pmd.load_file(C+".prmtop",C+".rst7"); res=cp.residues[601]
for tag,top in (("pb","system"),("docked","system_docked_cleaned")):
    s=pmd.load_file(T+top+".prmtop",T+top+".inpcrd"); t64={a.name:a for a in s.residues[601].atoms}
    assert len(pairs)==64 and set(pairs.values())==set(t64), (len(pairs),set(t64)-set(pairs.values()))
    mol=cp[":602"]
    for a in mol.atoms:
        t=t64[pairs[a.name]]; a.xx,a.xy,a.xz=t.xx,t.xy,t.xz
    # bond lengths of the CPD topology on the T64-derived coordinates
    bad=[]
    for a in res.atoms:
        for p in a.bond_partners:
            if p.residue.idx==601 and a.idx<p.idx:
                d=np.linalg.norm([mol.atoms[a.idx-res.atoms[0].idx].xx-mol.atoms[p.idx-res.atoms[0].idx].xx,
                                  mol.atoms[a.idx-res.atoms[0].idx].xy-mol.atoms[p.idx-res.atoms[0].idx].xy,
                                  mol.atoms[a.idx-res.atoms[0].idx].xz-mol.atoms[p.idx-res.atoms[0].idx].xz])
                bad.append((round(float(d),2),a.name,p.name))
    bad.sort(reverse=True); print(tag,"MOL bonds:",len(bad),"| >1.8 A:",[b for b in bad if b[0]>1.8],"| >1.7 A count:",sum(1 for b in bad if b[0]>1.7),"| min",bad[-1])
    mol.save(f"{tag}_mol_t64map.pdb",overwrite=True)
    s.save(f"{tag}_all.pdb",overwrite=True)
    L=[l for l in open(f"{tag}_all.pdb") if l.startswith(("ATOM","HETATM"))]
    start=sum(len(r.atoms) for r in s.residues[:601]); n=len(s.residues[601].atoms); cols=L[start][21:27]
    M=[l[:21]+cols+l[27:] for l in open(f"{tag}_mol_t64map.pdb") if l.startswith(("ATOM","HETATM"))]
    L=L[:start]+M+L[start+n:]
    names=[r.name for r in s.residues]; ters={565}|{i+1 for i,x in enumerate(names[:620]) if x.endswith("3") and x.startswith("D")}
    final=[];ri=0;prev=None
    for l in L:
        key=(l[17:20],l[21:27])
        if prev is not None and key!=prev:
            if ri in ters: final.append("TER\n")
            ri+=1
        elif prev is None: ri=1
        prev=key; final.append(l)
    final+=["TER\n","END\n"]; open(f"{tag}_t64map.pdb","w").writelines(final); print(tag,"written; residues",ri,"TER after",sorted(ters))
