"""Paneles de distribucion por persona de la figura 1: score de reloj y score glucocorticoide,
caso frente a control en las cinco cohortes, con los valores de P explicitos.
"""
import pandas as pd, numpy as np, sys, os, warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
OUT="results"; RAW="data/raw"
COLA,COLM,COLI="#2a9d8f","#d95f02","#6a3d9a"
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
ann=pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot_tsv.gz",sep="\t",usecols=["GeneID","Symbol"],dtype=str)
symmap=dict(zip(ann.GeneID,ann.Symbol))
def islets():
    C=pd.read_csv(f"{RAW}/GSE159984_raw_counts.tsv.gz",sep="\t",index_col=0)
    C=C[[g for g in C.columns if int(g[3:])<=4852118]]
    C.index=[symmap.get(str(i),str(i)) for i in C.index]; C=C[~C.index.str.isdigit()]; C=C.groupby(level=0).sum(); C=C[(C>5).sum(1)>=10]
    g=pd.Series(["T2D" if int(c[3:])<=4852060 else "ND" for c in C.columns],index=C.columns)
    return np.log2(C/C.sum()*1e6+1),g
COH=[]
X=pd.read_pickle(f"{OUT}/new/GSE25462_entrez.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE25462_meta.pkl").set_index("gsm").reindex(X.columns)
k=md.index[md.fh.isin(["DM","Family history negative"])]
COH.append(("GSE25462\nmúsculo",X[list(k)],(md.fh.reindex(k)=="DM").values,COLM))
X=pd.read_pickle(f"{OUT}/new/GSE182120_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE182120_meta.pkl").set_index("gsm").reindex(X.columns)
COH.append(("GSE182120\nmúsculo",X,(md.disease=="T2D").values,COLM))
X=pd.read_pickle(f"{OUT}/new/GSE26637_entrez.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE26637_meta.pkl").set_index("gsm").reindex(X.columns)
k=md.index[md.stim=="fasting"]; COH.append(("GSE26637\nadiposo",X[list(k)],(md.status.reindex(k)=="resistant").values,COLA))
X=pd.read_pickle(f"{OUT}/new/GSE78721_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE78721_meta.pkl").set_index("gsm").reindex(X.columns)
COH.append(("GSE78721\nadiposo",X,(md.dx!="Normal").values,COLA))
Xi,gi=islets(); COH.append(("GSE159984\nislotes",Xi,(gi=="T2D").values,COLI))
fig,ax=plt.subplots(1,2,figsize=(13,4.4))
rng=np.random.default_rng(0)
for j,(setname,a) in enumerate((("Clock output",ax[0]),("Glucocorticoid",ax[1]))):
    for i,(name,X,grp,c) in enumerate(COH):
        s=np.asarray(sc(X,SETS[setname]["genes"]),float)
        for k2,(m,lab) in enumerate(((~grp,"ctrl"),(grp,"caso"))):
            v=s[m]; xx=i*2.6+k2*0.9
            a.scatter(np.full(v.shape,xx)+rng.uniform(-.2,.2,len(v)),v,s=10,alpha=.55,
                      color=c,marker="o" if k2 else "s",edgecolors="none")
            a.hlines(np.median(v),xx-.3,xx+.3,color="k",lw=1.6)
        p=stats.mannwhitneyu(s[grp],s[~grp]).pvalue
        a.text(i*2.6+.45,a.get_ylim()[1] if False else 2.35,f"P = {p:.3f}",fontsize=7,ha="center")
    a.set_xticks([i*2.6+.45 for i in range(len(COH))]); a.set_xticklabels([c[0] for c in COH],fontsize=7.5)
    a.set(ylabel=f"score {'de salida del reloj' if j==0 else 'glucocorticoide'} (z)",ylim=(-2.6,2.7),
          title=("c  Distribución por persona: control (cuadrado) frente a caso (círculo)" if j==0 else "d  El eje glucocorticoide no separa en ninguna cohorte"))
    a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9)
fig.tight_layout(); fig.savefig(f"{OUT}/figs/F1_distribuciones.png",dpi=160,bbox_inches="tight")
print("ok")
