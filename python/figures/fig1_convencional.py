#!/usr/bin/env python3
"""Figura 1. Lo que la comparación de grupos ve, y lo que no puede ver."""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results")
BLUE,ORANGE,GREEN,GREY="#1f77b4","#d95f02","#2a9d8f","#999999"
# Valores publicados de las cohortes que no están en este entorno (correr 01b en el clúster para
# regenerarlos). Se marcan en la figura con borde punteado para que nunca se confundan con lo corrido.
PUB=pd.DataFrame([dict(cohorte="GSE159984 · músculo",tejido="músculo",casos=np.nan,controles=np.nan,genes_FDR01=0,fuente="publicado")])
PUBZ=pd.DataFrame(columns=["cohorte","conjunto","z"])
D=pd.read_csv(f"{OUT}/fig1_de_per_cohort.tsv",sep="\t"); D["fuente"]="corrido"
D=pd.concat([D,PUB[~PUB.cohorte.isin(D.cohorte)]],ignore_index=True)
S=pd.read_csv(f"{OUT}/fig1_signature_genes.tsv",sep="\t")
Z=pd.read_csv(f"{OUT}/fig1_set_enrichment.tsv",sep="\t")
PAIR=pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv",sep="\t")
fig,ax=plt.subplots(1,4,figsize=(17,4.3))
# a) tamaño de la lista diferencial
d=D.sort_values("genes_FDR01")
c=[GREEN if t=="adiposo" else ORANGE for t in d.tejido]
ax[0].barh(range(len(d)),d.genes_FDR01,color=c,
           edgecolor=["k" if f=="publicado" else "none" for f in d.fuente],
           linestyle="--",linewidth=1.0)
for i,(n,ca,co) in enumerate(zip(d.genes_FDR01,d.casos,d.controles)):
    tag="" if pd.isna(ca) else f"  ({int(ca)} vs {int(co)})"
    ax[0].text(max(n,2)*1.4+0.6,i,f"{int(n)}{tag}",va="center",fontsize=8)
ax[0].set_yticks(range(len(d))); ax[0].set_yticklabels(d.cohorte,fontsize=8)
ax[0].set(xlim=(0,22),xlabel="genes a FDR < 0,1",
          title="a  Cinco cohortes, casi ningún gen\n(borde punteado: valor publicado)")
# b) los 18 genes de las firmas
S=S[S.cohorte.isin(D.cohorte)]
piv=S.pivot_table(index="gen",columns="cohorte",values="t")
ordr=[g for g in SETS["Clock output"]["genes"]+SETS["Glucocorticoid"]["genes"] if g in piv.index]
piv=piv.loc[ordr]
x=np.arange(len(piv))
for k,c_ in enumerate(piv.columns):
    col=GREEN if "adiposo" in c_ else (ORANGE if "13070" in c_ else GREY)
    ax[1].scatter(x+(k-1)*0.22,piv[c_],s=26,color=col,label=c_,zorder=3)
ax[1].axhline(0,color="k",lw=.8)
for v in (-2,2): ax[1].axhline(v,color="k",lw=.6,ls=":")
ax[1].axvline(len(SETS["Clock output"]["genes"])-0.5,color="k",lw=.8,ls="--")
ax[1].set_xticks(x); ax[1].set_xticklabels(piv.index,rotation=90,fontsize=7)
ax[1].set(ylabel="t moderada (caso − control)",
          title="b  Gen a gen no se obtiene la firma:\nningún gen pasa FDR 0,1")
ax[1].legend(fontsize=6.5,frameon=False,loc="lower left")
ax[1].text(1.5,ax[1].get_ylim()[1]*.92,"salida del reloj",fontsize=7,ha="center")
ax[1].text(len(ordr)-4,ax[1].get_ylim()[1]*.92,"glucocorticoide",fontsize=7,ha="center")
# c) enriquecimiento de conjunto
p=Z.pivot_table(index="cohorte",columns="conjunto",values="z").reindex(columns=["Clock output","Glucocorticoid","Insulin core","Oxidative"]).astype(float)
im=ax[2].imshow(p.values,cmap="RdBu_r",vmin=-8,vmax=8,aspect="auto")
ax[2].set_xticks(range(p.shape[1])); ax[2].set_xticklabels(["salida\nreloj","gluco-\ncorticoide","insulina","oxidativo"],fontsize=8)
ax[2].set_yticks(range(p.shape[0])); ax[2].set_yticklabels(p.index,fontsize=8)
for i in range(p.shape[0]):
    for j in range(p.shape[1]):
        v=p.values[i,j]
        ax[2].text(j,i,"—" if np.isnan(v) else f"{v:.1f}",ha="center",va="center",fontsize=8,
                   color="w" if (not np.isnan(v) and abs(v)>4) else "k")
ax[2].set(title="c  A nivel de conjunto el reloj YA sale alterado\nen músculo y en adiposo; el eje GC no")
plt.colorbar(im,ax=ax[2],label="z",shrink=.8)
# d) lo que sí hay, y que la comparación de grupos no puede ver
sub=PAIR[(PAIR.contraste=="NR3C1~reloj")&(PAIR.cohorte.isin(["GSE13070","GSE182120"]))].copy()
if "GSE182120" not in set(sub.cohorte):
    import pickle
    X=pd.read_pickle(f"{OUT}/new/GSE182120_log2.pkl")
    def _sc(X,g):
        gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
    r=stats.spearmanr(X.loc["NR3C1"],_sc(X,SETS["Clock output"]["genes"]))
    n=X.shape[1]; z=np.arctanh(r.statistic); se=1/np.sqrt(n-3)
    sub=pd.concat([sub,pd.DataFrame([dict(contraste="NR3C1~reloj",cohorte="GSE182120",tejido="músculo",n=n,rho=r.statistic,z=z,se=se)])],ignore_index=True)
sub=sub.sort_values("rho",ascending=False)
lab=[f"{r.cohorte}\n{r.tejido}\n(n={int(r.n)})" for _,r in sub.iterrows()]
lo=np.tanh(sub.z-1.96*sub.se); hi=np.tanh(sub.z+1.96*sub.se)
cc=[GREEN if "adiposo" in t else ORANGE for t in sub.tejido]
ax[3].bar(range(len(sub)),sub.rho,color=cc,yerr=[sub.rho-lo,hi-sub.rho],capsize=4)
ax[3].axhline(0,color="k",lw=.8); ax[3].set_xticks(range(len(sub))); ax[3].set_xticklabels(lab,fontsize=8)
ax[3].set(ylabel="ρ receptor ~ salida del reloj",ylim=(0,.72),
          title="d  La señal que falta está entre personas,\nno entre grupos")
ax[3].text(len(sub)/2-.5,.62,"el conjunto del reloj sale alterado (c),\npero no qué señal lo explica",fontsize=7.5,ha="center",style="italic")
for a in fig.axes: a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9)
fig.tight_layout(); fig.savefig(f"{OUT}/figs/F1_convencional.png",dpi=160,bbox_inches="tight")
fig.savefig(f"{OUT}/figs/F1_convencional.pdf",bbox_inches="tight"); print("ok")
