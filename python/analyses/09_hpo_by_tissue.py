import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import pandas as pd, numpy as np, sys, gzip, warnings; warnings.filterwarnings("ignore")
from scipy import stats
from gene_sets import SETS, resolve
ann=pd.read_csv(f"{_RAW}/Human_GRCh38_p13_annot_tsv.gz",sep="\t",usecols=["GeneID","Symbol"],dtype=str)
sym=dict(zip(ann.GeneID,ann.Symbol))
import islets as _isl
def islets():
    """Islotes de donantes de GSE159984, log2 CPM. Que GSM son donantes se lee de la matriz de serie."""
    don,_=_isl.groups(_RAW)
    return _isl.log2cpm_symbols(f"{_RAW}/GSE159984_raw_counts.tsv.gz",don,f"{_RAW}/Human_GRCh38_p13_annot_tsv.gz")
def rank_vs(X,cv=None):
    Xc=X.sub(X.mean(0),axis=1)                    # centrado por muestra (§30)
    D=np.column_stack([np.ones(Xc.shape[1])]+(cv or []))
    Q,_=np.linalg.qr(D); V=Xc.values.T; Vr=V-Q@(Q.T@V)
    rk=np.apply_along_axis(stats.rankdata,0,Vr)
    nr=rk[:,list(Xc.index).index("NR3C1")]
    z=(rk-rk.mean(0))/rk.std(0); nz=(nr-nr.mean())/nr.std()
    return pd.Series((z*nz[:,None]).mean(0),index=Xc.index).dropna()
A=pd.read_pickle(f"{_OUT}/new/AAGMEx_adipose_log2.pkl"); M=pd.read_pickle(f"{_OUT}/new/AAGMEx_muscle_log2.pkl")
with gzip.open(f"{_RAW}/GSE95674_series_matrix.txt.gz","rt") as f:
    for l in f:
        if l.startswith("!Sample_title"): t=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
cov=pd.DataFrame({"sex":sx,"age":ag},index=t)
RK={}
RK["adipose"]=rank_vs(A,[(cov.sex.reindex(A.columns)=="M").astype(float).values,cov.age.reindex(A.columns).values,A.loc["LEP"].values,A.loc["CD68"].values])
RK["muscle"]=rank_vs(M,[(cov.sex.reindex(M.columns)=="M").astype(float).values,cov.age.reindex(M.columns).values])
RK["islets"]=rank_vs(islets())
syn=pd.read_csv(f"{_OUT}/new/ncbi_synonym_map.tsv",sep="\t",index_col=0)["current"].to_dict()
rng=np.random.default_rng(0); rows=[]
gmt={}
_hpo=f"{_RAW}/msigdb/c5.hpo.v2026.2.Hs.symbols.gmt"      # donde lo deja download_public.sh
for l in open(_hpo if _os.path.exists(_hpo) else "data/msigdb/c5_hpo_v2026_2_Hs_symbols.gmt"):
    p=l.rstrip("\n").split("\t"); gmt[p[0]]=set(p[2:])
WANT=["HP_TRUNCAL_OBESITY","HP_INCREASED_URINARY_CORTISOL_LEVEL","HP_INSULIN_RESISTANCE",
      "HP_HYPERGLYCEMIA","HP_TYPE_II_DIABETES_MELLITUS","HP_OBESITY","HP_ABNORMAL_CIRCULATING_CORTISOL_CONCENTRATION",
      "HP_HYPERCORTISOLISM","HP_PRIMARY_HYPERCORTISOLISM","HP_DIABETES_MELLITUS","HP_HYPERINSULINEMIA",
      "HP_ABNORMAL_GLUCOSE_HOMEOSTASIS","HP_INCREASED_BODY_WEIGHT","HP_MUSCLE_WEAKNESS","HP_HYPERTENSION"]
# El ranking es la correlacion de cada gen con NR3C1, y la de NR3C1 consigo mismo vale 1. Varios de
# estos terminos lo contienen (entre ellos el de cortisol urinario, de 16 genes): contarlo es
# circular. `z` lo excluye; `z_including_NR3C1` es el calculo anterior, tal como se publico.
absent=[w for w in WANT if w not in gmt]
if absent: print("terminos que no existen en esta version de MSigDB y no se puntuan:", ", ".join(absent))
rng2=np.random.default_rng(1)
for tis,R in RK.items():
    R.index=[syn.get(g,g) for g in R.index]; R=R.groupby(level=0).mean()
    idx={g:i for i,g in enumerate(R.index)}; vals=R.values; cache={}; cache2={}
    ia=idx.get("NR3C1"); pool=np.array([i for i in range(len(vals)) if i!=ia])
    for name in WANT:
        if name not in gmt: continue
        g=[idx[x] for x in gmt[name] if x in idx]; n=len(g)
        if n<8: continue
        if n not in cache:
            m=np.array([vals[rng.integers(0,len(vals),n)].mean() for _ in range(3000)]); cache[n]=(m.mean(),m.std())
        mu,sd=cache[n]; z0=(vals[g].mean()-mu)/sd; z=z0; has=ia in g
        if has:
            g2=[i for i in g if i!=ia]; k=len(g2)
            if k not in cache2:
                m=np.array([vals[rng2.choice(pool,k)].mean() for _ in range(3000)]); cache2[k]=(m.mean(),m.std())
            z=(vals[g2].mean()-cache2[k][0])/cache2[k][1]
        rows.append(dict(tissue=tis,term=name,n=n,z=z,z_including_NR3C1=z0,contains_NR3C1=has))
T=pd.DataFrame(rows); T.to_csv(f"{_OUT}/new/HPO_by_tissue.tsv",sep="\t",index=False)
P=T.pivot_table(index="term",columns="tissue",values="z").reindex(columns=["adipose","muscle","islets"])
print(P.round(2).to_string())
