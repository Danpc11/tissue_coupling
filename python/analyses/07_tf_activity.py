"""Actividad de factores de transcripcion con regulones DoRothEA A-C.

La actividad es la media ponderada por el signo de la interaccion de los z de las dianas. Los
genes del reloj se excluyen de cada regulon antes de correlacionar con la salida del reloj, porque
si no el resultado es circular.
Salidas: results/new/AAGMEx_TF_activity_vs_clock.tsv, NR3C1_activity_vs_mRNA.tsv
"""
import pandas as pd, numpy as np, sys, gzip, warnings
warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
d=pd.read_csv("data/raw/dorothea_ABC.tsv",sep="\t",low_memory=False,usecols=["source_genesymbol","target_genesymbol","is_stimulation","is_inhibition","dorothea_level"])
d["mor"]=np.where(d.is_stimulation&~d.is_inhibition,1,np.where(d.is_inhibition&~d.is_stimulation,-1,0))
d=d[d.mor!=0].drop_duplicates(["source_genesymbol","target_genesymbol"])
CL=set(SETS["Clock output"]["genes"])|set(SETS["Core oscillator"]["genes"])|{"PER1","PER2","PER3","CIART","ARNTL","BMAL1"}
A=pd.read_pickle("results/new/AAGMEx_adipose_log2.pkl")
with gzip.open("data/raw/GSE95674_series_matrix.txt.gz","rt") as f:
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
R.to_csv("results/new/AAGMEx_TF_activity_vs_clock.tsv",sep="\t",index=False)
print(R.head(12).to_string(index=False)); print("..."); print(R.tail(10).to_string(index=False))
print(f"\nNR3C1: puesto {list(R.TF).index('NR3C1')+1} de {len(R)}, ρ {R[R.TF=='NR3C1'].rho.iloc[0]:+.3f} con {R[R.TF=='NR3C1'].n_dianas.iloc[0]} dianas")
print("regulón NR3C1 de DoRothEA A-C:",sorted(d[d.source_genesymbol=='NR3C1'].target_genesymbol))
print("¿incluye alguna diana GC clásica?",[g for g in ["FKBP5","TSC22D3","ZBTB16","KLF15","DDIT4","SESN1","GLUL","SGK1","CEBPD","PER1"] if g in set(d[d.source_genesymbol=='NR3C1'].target_genesymbol)])
