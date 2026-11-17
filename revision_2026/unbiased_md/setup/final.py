import pickle, numpy as np, parmed as pmd
T="/media/rahul/EXTERNAL_USB/binding_affinity/TI/ti_rep1/bundle/ti_rep1/topologies/"
pairs=dict((m,t) for t,m in pickle.load(open("pairs.pkl","rb"))); pairs.update({"O14":"O8","H20":"H12"})
for tag,top in (("pb","system"),("docked","system_docked_cleaned")):
    ref=pmd.load_file(T+top+".prmtop",T+top+".inpcrd"); new=pmd.load_file(f"cpd_{tag}.prmtop",f"cpd_{tag}.leap.inpcrd")
    ls=sum(len(r.atoms) for r in ref.residues[:601]); le=ls+64; idx=np.array([i for i in range(len(ref.atoms)) if not(ls<=i<le)])
    X=np.array([[a.xx,a.xy,a.xz] for a in new.atoms]); R=np.array([[a.xx,a.xy,a.xz] for a in ref.atoms])
    sh=X[idx]-R[idx]; print("==",tag,"leap output offset vs 64pp (non-lesion): mean",sh.mean(0).round(4),"std",sh.std(0).round(5),"max dev from mean %.5f"%np.abs(sh-sh.mean(0)).max())
    s0=sh.mean(0)
    X[ls:le]-=s0                      # bring the lesion into the 64pp frame
    X[idx]=R[idx]                     # non-lesion atoms: exact 64pp coordinates
    new.coordinates=X; new.box=ref.box; new.save(f"cpd_{tag}_final.rst7",format="rst7",overwrite=True)
    chk=pmd.load_file(f"cpd_{tag}.prmtop",f"cpd_{tag}_final.rst7"); Y=np.array([[a.xx,a.xy,a.xz] for a in chk.atoms])
    t64={a.name:a for a in ref.residues[601].atoms}
    d=[np.linalg.norm(Y[ls+i]-np.array([t64[pairs[a.name]].xx,t64[pairs[a.name]].xy,t64[pairs[a.name]].xz])) for i,a in enumerate(chk.residues[601].atoms)]
    print("   lesion vs mapped T64 positions: max %.4f A | non-lesion max |dxyz| vs 64pp %.1e | box %s"%(max(d),np.abs(Y[idx]-R[idx]).max(),np.round(chk.box[:3],6)))
    # clashes of lesion with surroundings (waters/DNA/protein), minimum-image
    box=np.array(ref.box[:3]); owner=np.array([a.residue.idx for a in chk.atoms]); worst=[]
    for i in range(ls,le):
        dd=Y-Y[i]; dd-=box*np.round(dd/box); r=np.linalg.norm(dd,axis=1); r[ls:le]=99; r[[j for j in range(len(r)) if abs(owner[j]-601)==1 and chk.atoms[j].name in("O3'","P","O5'","C3'","C5'","OP1","OP2","O1P","O2P") and np.linalg.norm(dd[j])<2.2]]=99
        j=int(np.argmin(r)); worst.append((round(float(r[j]),2),chk.atoms[i].name,chk.atoms[j].residue.name+str(chk.atoms[j].residue.idx+1)+"."+chk.atoms[j].name))
    worst.sort(); print("   closest lesion contacts:",worst[:5],"| <1.5 A:",sum(1 for w in worst if w[0]<1.5))
    # lesion bonded parameters vs original CPD topology
    C=pmd.load_file("/home/rahul/Downloads/CPD_DNA_DAMAGE_data_files/NAR_paper_data_files/INITIAL_STRUCTURE/3_cpd_pb.wb.prmtop"); o0=C.residues[601].atoms[0].idx; n0=chk.residues[601].atoms[0].idx
    kb=lambda s,b0:{tuple(sorted((x.atom1.name,x.atom2.name))):(round(x.type.k,2),round(x.type.req,3)) for x in s.bonds if all(a.residue.idx==601 for a in (x.atom1,x.atom2))}
    bn,bo=kb(chk,0),kb(C,0); print("   lesion bond params vs original CPD topology: %d bonds, mismatches %d"%(len(bn),sum(1 for k in bn if bo.get(k)!=bn[k])))
    ka=lambda s:{tuple(sorted((x.atom1.name,x.atom3.name)))+(x.atom2.name,):(round(x.type.k,2),round(x.type.theteq,2)) for x in s.angles if all(a.residue.idx==601 for a in (x.atom1,x.atom2,x.atom3))}
    an,ao=ka(chk),ka(C); print("   lesion angle params: %d angles, mismatches %d"%(len(an),sum(1 for k in an if ao.get(k)!=an[k])))
