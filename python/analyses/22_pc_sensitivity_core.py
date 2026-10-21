#!/usr/bin/env python3
"""Sensibilidad del eje receptor-reloj al ajuste por componentes principales, en las cohortes
del argumento central.

Quitar los primeros componentes de la matriz elimina la heterogeneidad entre muestras, pero una
senal sistemica que afecta a muchos genes a la vez vive tambien ahi: el ajuste es a la vez control
contra composicion y posible sobrecorreccion. Se reporta como analisis de sensibilidad, no como
resultado principal.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results"); RAW=os.environ.get("T2D_RAW","data/raw")

def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]
    Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)

def islets():
    ann=pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot_tsv.gz",sep="\t",usecols=["GeneID","Symbol"],dtype=str)
    sym=dict(zip(ann.GeneID,ann.Symbol))
    C=pd.read_csv(f"{RAW}/GSE159984_raw_counts.tsv.gz",sep="\t",index_col=0)
    C=C[[g for g in C.columns if int(g[3:])<=4852118]]
    C.index=[sym.get(str(i),str(i)) for i in C.index]; C=C[~C.index.str.isdigit()]
    C=C.groupby(level=0).sum(); C=C[(C>5).sum(1)>=10]
    L=np.log2(C/C.sum()*1e6+1); return L.sub(L.mean(0),axis=1)

def run(X, cov, label, npc=(0,1,2,3,5)):
    Z=X.sub(X.mean(1),axis=0).div(X.std(1)+1e-9,axis=0)
    U,S,Vt=np.linalg.svd(Z.values-Z.values.mean(1,keepdims=True),full_matrices=False)
    var=S**2/(S**2).sum()
    out=[]
    for k in npc:
        cv=list(cov)+[Vt[j] for j in range(k)]
        D=np.column_stack([np.ones(X.shape[1])]+cv) if cv else np.ones((X.shape[1],1))
        Q,_=np.linalg.qr(D); r=lambda y: np.asarray(y,float)-Q@(Q.T@np.asarray(y,float))
        co=r(sc(X,SETS["Clock output"]["genes"]))
        rho=stats.spearmanr(r(X.loc["NR3C1"]),co).statistic
        rng=np.random.default_rng(0)
        null=np.array([stats.spearmanr(r(X.loc[g]),co).statistic
                       for g in rng.choice(X.index,800,replace=False)])
        out.append(dict(cohort=label,n=X.shape[1],n_pc=k,pc1_var=var[0]*100,
                        rho=rho,null_sd=null.std(),pct=np.mean(null<rho)*100))
        print(f"  {label:22s} PC quitados {k}: ρ {rho:+.3f}  nulo sd {null.std():.3f}  percentil {np.mean(null<rho)*100:3.0f}%")
    return out

if __name__=="__main__":
    import gzip
    rows=[]
    with gzip.open(f"{RAW}/GSE95674_series_matrix.txt.gz","rt",errors="ignore") as f:
        for l in f:
            if l.startswith("!Sample_geo_accession"): ids=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
    cov0=pd.DataFrame({"sex":sx,"age":ag},index=ids)
    for pk,lab in (("AAGMEx_adipose_log2","AAGMEx adipose"),("AAGMEx_muscle_log2","AAGMEx muscle")):
        X=pd.read_pickle(f"{OUT}/new/{pk}.pkl"); X=X.sub(X.mean(0),axis=1)
        c=cov0.reindex(X.columns)
        # las columnas de AAGMEx no siempre son GSM: si no casan, se corre sin covariables
        if c.sex.isna().all():
            cv=[]; print(f"  [{lab}] sin covariables: los identificadores no casan con la serie matriz")
        else:
            cv=[(c.sex=="M").astype(float).values, c.age.values.astype(float)]
            cv=[np.nan_to_num(v, nan=np.nanmean(v)) for v in cv]
        rows+=run(X,cv,lab)
    for pk,lab in (("GSE41168_adipose","GSE41168 adipose"),("GSE41168_muscle","GSE41168 muscle")):
        if os.path.exists(f"{OUT}/new/{pk}.pkl"):
            X=pd.read_pickle(f"{OUT}/new/{pk}.pkl"); X=X.sub(X.mean(0),axis=1)
            rows+=run(X,[],lab)
    rows+=run(islets(),[],"GSE159984 islets")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/pc_sensitivity_core.tsv",sep="\t",index=False)
