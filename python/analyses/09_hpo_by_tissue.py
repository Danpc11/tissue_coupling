import pandas as pd, numpy as np, sys, gzip, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
ann=pd.read_csv("data/raw/Human_GRCh38_p13_annot_tsv.gz",sep="\t",usecols=["GeneID","Symbol"],dtype=str)
sym=dict(zip(ann.GeneID,ann.Symbol))
def islets():
    C=pd.read_csv("data/raw/GSE159984_raw_counts.tsv.gz",sep="\t",index_col=0)
    C=C[[g for g in C.columns if int(g[3:])<=4852118]]
    C.index=[sym.get(str(i),str(i)) for i in C.index]; C=C[~C.index.str.isdigit()]
    C=C.groupby(level=0).sum(); C=C[(C>5).sum(1)>=10]
    return np.log2(C/C.sum()*1e6+1)
def rank_vs(X,cv=None):
    Xc=X.sub(X.mean(0),axis=1)                    # centrado por muestra (§30)
    D=np.column_stack([np.ones(Xc.shape[1])]+(cv or []))
    Q,_=np.linalg.qr(D); V=Xc.values.T; Vr=V-Q@(Q.T@V)
    rk=np.apply_along_axis(stats.rankdata,0,Vr)
    nr=rk[:,list(Xc.index).index("NR3C1")]
    z=(rk-rk.mean(0))/rk.std(0); nz=(nr-nr.mean())/nr.std()
    return pd.Series((z*nz[:,None]).mean(0),index=Xc.index).dropna()
A=pd.read_pickle("results/new/AAGMEx_adipose_log2.pkl"); M=pd.read_pickle("results/new/AAGMEx_muscle_log2.pkl")
with gzip.open("data/raw/GSE95674_series_matrix.txt.gz","rt") as f:
    for l in f:
        if l.startswith("!Sample_title"): t=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
cov=pd.DataFrame({"sex":sx,"age":ag},index=t)
RK={}
RK["adipose"]=rank_vs(A,[(cov.sex.reindex(A.columns)=="M").astype(float).values,cov.age.reindex(A.columns).values,A.loc["LEP"].values,A.loc["CD68"].values])
RK["muscle"]=rank_vs(M,[(cov.sex.reindex(M.columns)=="M").astype(float).values,cov.age.reindex(M.columns).values])
RK["islets"]=rank_vs(islets())
syn=pd.read_csv("results/new/ncbi_synonym_map.tsv",sep="\t",index_col=0)["current"].to_dict()
rng=np.random.default_rng(0); rows=[]
gmt={}
for l in open("data/msigdb/c5_hpo_v2026_2_Hs_symbols.gmt"):
    p=l.rstrip("\n").split("\t"); gmt[p[0]]=set(p[2:])
WANT=["HP_TRUNCAL_OBESITY","HP_INCREASED_URINARY_CORTISOL_LEVEL","HP_INSULIN_RESISTANCE",
      "HP_HYPERGLYCEMIA","HP_TYPE_II_DIABETES_MELLITUS","HP_OBESITY","HP_ABNORMAL_CIRCULATING_CORTISOL_CONCENTRATION",
      "HP_HYPERCORTISOLISM","HP_PRIMARY_HYPERCORTISOLISM","HP_DIABETES_MELLITUS","HP_HYPERINSULINEMIA",
      "HP_ABNORMAL_GLUCOSE_HOMEOSTASIS","HP_INCREASED_BODY_WEIGHT","HP_MUSCLE_WEAKNESS","HP_HYPERTENSION"]
for tis,R in RK.items():
    R.index=[syn.get(g,g) for g in R.index]; R=R.groupby(level=0).mean()
    idx={g:i for i,g in enumerate(R.index)}; vals=R.values; cache={}
    for name in WANT:
        if name not in gmt: continue
        g=[idx[x] for x in gmt[name] if x in idx]; n=len(g)
        if n<8: continue
        if n not in cache:
            m=np.array([vals[rng.integers(0,len(vals),n)].mean() for _ in range(3000)]); cache[n]=(m.mean(),m.std())
        mu,sd=cache[n]; rows.append(dict(tissue=tis,term=name,n=n,z=(vals[g].mean()-mu)/sd))
T=pd.DataFrame(rows); T.to_csv("results/new/HPO_by_tissue.tsv",sep="\t",index=False)
P=T.pivot_table(index="term",columns="tissue",values="z").reindex(columns=["adipose","muscle","islets"])
print(P.round(2).to_string())
