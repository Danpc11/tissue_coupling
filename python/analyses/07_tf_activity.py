"""Actividad de factores de transcripcion con regulones DoRothEA A-C.

La actividad es la media ponderada por el signo de la interaccion de los z de las dianas. Los
genes del reloj se excluyen de cada regulon antes de correlacionar con la salida del reloj, porque
si no el resultado es circular.
Salidas: results/new/AAGMEx_TF_activity_vs_clock.tsv, NR3C1_activity_vs_mRNA.tsv
"""
import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import os
import pandas as pd, numpy as np, sys, gzip, warnings
warnings.filterwarnings("ignore")
from scipy import stats
from gene_sets import SETS, resolve
RAW=os.environ.get("T2D_RAW","data/raw")
d=pd.read_csv(f"{RAW}/dorothea_ABC.tsv",sep="\t",low_memory=False,usecols=["source_genesymbol","target_genesymbol","is_stimulation","is_inhibition","dorothea_level"])
d["mor"]=np.where(d.is_stimulation&~d.is_inhibition,1,np.where(d.is_inhibition&~d.is_stimulation,-1,0))
d=d[d.mor!=0].drop_duplicates(["source_genesymbol","target_genesymbol"])
CL=set(SETS["Clock output"]["genes"])|set(SETS["Core oscillator"]["genes"])|{"BMAL1"}   # BMAL1: simbolo vigente de ARNTL
A=pd.read_pickle(f"{_OUT}/new/AAGMEx_adipose_log2.pkl"); A=A.sub(A.mean(0),axis=1)   # centrado por muestra (Methods)
with gzip.open(f"{_RAW}/GSE95674_series_matrix.txt.gz","rt") as f:
    for l in f:
        if l.startswith("!Sample_title"): t=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
cov=pd.DataFrame({"sex":sx,"age":ag},index=t).reindex(A.columns)
base=[(cov.sex=="M").astype(float).values,cov.age.values,A.loc["LEP"].values]
def pres(y,cv):
    D=np.column_stack([np.ones(len(y))]+cv); b,*_=np.linalg.lstsq(D,np.asarray(y,float),rcond=None); return np.asarray(y,float)-D@b
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
co=pres(sc(A,SETS["Clock output"]["genes"]),base)
def act(tf,excl_clock=True):
    r=d[d.source_genesymbol==tf]
    g=[(resolve(x,A.index),m) for x,m in zip(r.target_genesymbol,r.mor) if not (excl_clock and x in CL)]
    g=[(a,b) for a,b in g if a is not None]
    if len(g)<8: return None,0
    Z=A.loc[[a for a,_ in g]]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    w=np.array([b for _,b in g]); return (Z.values*w[:,None]).sum(0)/np.sqrt(len(w)), len(w)
print("TF activity vs salida del reloj, AAGMEx adiposo, SIN genes del reloj en el regulón y ajustado por sexo/edad/LEP:")
res=[]
for tf in d.source_genesymbol.unique():
    a,n=act(tf)
    if a is not None: res.append((tf,n,stats.spearmanr(pres(a,base),co).statistic))
R=pd.DataFrame(res,columns=["TF","n_dianas","rho"]).sort_values("rho",ascending=False)
# Fondo: la actividad de un regulon es una suma de genes, y una suma de genes cualesquiera ya
# correlaciona algo con la puntuacion. Para cada tamano de regulon se sortean conjuntos de dianas al
# azar (mismos signos, genes fuera del reloj) y se da la posicion de la rho observada en ese nulo.
rng=np.random.default_rng(0); pool=np.array([g for g in A.index if g not in CL]); NULLN=500; cache={}
Zall=((A.T-A.mean(1))/(A.std(1)+1e-9)).T
def null_for(n):
    if n not in cache:
        v=[]
        for _ in range(NULLN):
            g=rng.choice(pool,n,replace=False); w=rng.choice([1.0,-1.0],n)
            v.append(stats.spearmanr(pres((Zall.loc[g].values*w[:,None]).sum(0)/np.sqrt(n),base),co).statistic)
        cache[n]=np.array(v)
    return cache[n]
R["null_sd"]=[null_for(int(n)).std() for n in R.n_dianas]
R["z_vs_random_targets"]=[(r-null_for(int(n)).mean())/null_for(int(n)).std() for r,n in zip(R.rho,R.n_dianas)]
R["p_vs_random_targets"]=[(np.sum(np.abs(null_for(int(n)))>=abs(r))+1)/(NULLN+1) for r,n in zip(R.rho,R.n_dianas)]
R.to_csv(f"{_OUT}/new/AAGMEx_TF_activity_vs_clock.tsv",sep="\t",index=False)
print(R.head(12).to_string(index=False)); print("..."); print(R.tail(10).to_string(index=False))
print(f"\nfrente a regulones al azar del mismo tamano: {int((R.p_vs_random_targets<0.05).sum())} de {len(R)} factores con P < 0.05 (se esperan {0.05*len(R):.0f} por azar); "
      f"FOXO3 z = {R[R.TF=='FOXO3'].z_vs_random_targets.iloc[0]:+.1f}, P = {R[R.TF=='FOXO3'].p_vs_random_targets.iloc[0]:.3f}; "
      f"NR3C1 z = {R[R.TF=='NR3C1'].z_vs_random_targets.iloc[0]:+.1f}, P = {R[R.TF=='NR3C1'].p_vs_random_targets.iloc[0]:.3f}")
print(f"\nNR3C1: puesto {list(R.TF).index('NR3C1')+1} de {len(R)}, ρ {R[R.TF=='NR3C1'].rho.iloc[0]:+.3f} con {R[R.TF=='NR3C1'].n_dianas.iloc[0]} dianas")
print("regulón NR3C1 de DoRothEA A-C:",sorted(d[d.source_genesymbol=='NR3C1'].target_genesymbol))
print("¿incluye alguna diana GC clásica?",[g for g in SETS["Glucocorticoid"]["genes"]+["PER1"] if g in set(d[d.source_genesymbol=='NR3C1'].target_genesymbol)])
