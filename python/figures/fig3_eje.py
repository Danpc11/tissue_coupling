#!/usr/bin/env python3
"""Figura 3. El reloj periférico sigue al receptor de glucocorticoides.

a  forest de 13 estimaciones cohorte-tejido, agrupadas por tejido, con los combinados
b  distribución de la salida del reloj por quintiles de NR3C1, cada punto una persona
c  señalización glucocorticoide frente a reloj, estratificado vivos / post mortem
d  diseño dentro del donante: cuánto del efecto es del individuo y cuánto del tejido
"""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results")
COL={"adiposo":"#2a9d8f","músculo":"#d95f02","islotes":"#6a3d9a"}
def fam(t): return "islotes" if "islote" in t else ("adiposo" if "adiposo" in t else "músculo")
def dl(sub):
    w=1/sub.se**2; m=(w*sub.z).sum()/w.sum(); Q=(w*(sub.z-m)**2).sum(); k=len(sub)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(sub.se**2+tau2)
    m2=(w2*sub.z).sum()/w2.sum(); se2=np.sqrt(1/w2.sum()); I2=max(0,(Q-(k-1))/Q*100) if Q>0 else 0
    return np.tanh(m2),np.tanh(m2-1.96*se2),np.tanh(m2+1.96*se2),I2,2*(1-stats.norm.cdf(abs(m2/se2))),k
R=pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv",sep="\t")
R["fam"]=R.tejido.apply(fam); R["lo"]=np.tanh(R.z-1.96*R.se); R["hi"]=np.tanh(R.z+1.96*R.se)
fig=plt.figure(figsize=(16,8.6)); gs=fig.add_gridspec(2,3,height_ratios=[1.25,1],hspace=.42,wspace=.3)
# ---- a
ax=fig.add_subplot(gs[0,0:2])
sub=R[R.contraste=="NR3C1~reloj"].copy()
sub["k"]=sub.fam.map({"islotes":0,"adiposo":1,"músculo":2}); sub=sub.sort_values(["k","rho"],ascending=[True,False]).reset_index(drop=True)
y=0; ticks=[]; labs=[]
for f in ("islotes","adiposo","músculo"):
    q=sub[sub.fam==f]
    for _,r in q.iterrows():
        ax.plot([r.lo,r.hi],[y,y],color=COL[f],lw=1.3); ax.scatter(r.rho,y,s=14+r.n/14,color=COL[f],zorder=3)
        ticks.append(y); labs.append(f"{r.cohorte} · {r.tejido} (n={int(r.n)})"); y+=1
    v=dl(q); ax.plot([v[1],v[2]],[y,y],color=COL[f],lw=2.6); ax.scatter(v[0],y,marker="D",s=60,color=COL[f],zorder=3)
    ticks.append(y); labs.append(f"{f}: ρ {v[0]:.2f}  (k={v[5]}, I²={v[3]:.0f}%)"); y+=1.4
v=dl(sub); ax.plot([v[1],v[2]],[y,y],color="k",lw=3); ax.scatter(v[0],y,marker="D",s=75,color="k",zorder=3)
ticks.append(y); labs.append(f"global: ρ {v[0]:.2f} (IC95 {v[1]:.2f}–{v[2]:.2f}), P = {v[4]:.0e}")
ax.axvline(0,color="k",lw=.8); ax.axvline(v[0],color="k",ls=":",lw=.7)
ax.set_yticks(ticks); ax.set_yticklabels(labs,fontsize=8); ax.invert_yaxis()
ax.set(xlabel="ρ receptor (NR3C1) ~ salida del reloj",xlim=(-0.2,1.0),
       title="a  Trece estimaciones, cuatro tejidos, seis cohortes")
# ---- b distribuciones por quintil, incluyendo islotes
ax=fig.add_subplot(gs[0,2]); rng=np.random.default_rng(0)
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
ann=pd.read_csv(f"{os.environ.get('T2D_RAW','data/raw')}/Human_GRCh38_p13_annot_tsv.gz",sep="\t",usecols=["GeneID","Symbol"],dtype=str)
symmap=dict(zip(ann.GeneID,ann.Symbol))
C=pd.read_csv(f"{os.environ.get('T2D_RAW','data/raw')}/GSE159984_raw_counts.tsv.gz",sep="\t",index_col=0)
C=C[[g for g in C.columns if int(g[3:])<=4852118]]
C.index=[symmap.get(str(i),str(i)) for i in C.index]; C=C[~C.index.str.isdigit()]; C=C.groupby(level=0).sum(); C=C[(C>5).sum(1)>=10]
I=np.log2(C/C.sum()*1e6+1)
panels=[("islotes (n=85)",I.loc["NR3C1"],sc(I,SETS["Clock output"]["genes"]),COL["islotes"]),
        ("GTEx adiposo visceral (n=587)",*[None,None],COL["adiposo"])]
G=pd.read_pickle(f"{OUT}/new/GTEx_adipose_visc.pkl")
panels[1]=("GTEx visceral (n=587)",G.loc["NR3C1"],sc(G,SETS["Clock output"]["genes"]),COL["adiposo"])
off=0
for name,x,yv,c in panels:
    x=np.asarray(x,float); yv=np.asarray(yv,float); q=pd.qcut(x,5,labels=False)
    for k in range(5):
        v=yv[q==k]-np.median(yv)
        ax.scatter(np.full(v.shape,k+off)+rng.uniform(-.22,.22,len(v)),v,s=5,alpha=.4,color=c)
        ax.hlines(np.median(v),k+off-.33,k+off+.33,color="k",lw=1.8)
    ax.text(2+off,2.55,name,fontsize=7.5,ha="center",color=c)
    off+=6
ax.set_xticks(list(range(5))+list(range(6,11))); ax.set_xticklabels(["1","2","3","4","5"]*2,fontsize=7)
ax.set(xlabel="quintil de NR3C1 (bajo → alto)",ylabel="salida del reloj (centrada)",
       ylim=(-2.1,2.95),title="b  Cada punto es una persona:\nla mediana sube de forma monótona")
# ---- c GC vivos vs post mortem
ax=fig.add_subplot(gs[1,0])
s2=R[R.contraste=="GC~reloj"].copy(); s2["pm"]=(s2.cohorte=="GTEx").astype(int)
s2=s2.sort_values(["pm","rho"],ascending=[True,False]).reset_index(drop=True)
for i,r in s2.iterrows():
    c="#6a3d9a" if r.pm==0 else "#999999"
    ax.plot([r.lo,r.hi],[i,i],color=c,lw=1.3); ax.scatter(r.rho,i,s=14+r.n/14,color=c,zorder=3)
vv=dl(s2[s2.pm==0]); pm=dl(s2[s2.pm==1]); yy=len(s2)
for lab,v,c in ((f"vivos (k={vv[5]})",vv,"#6a3d9a"),(f"post mortem (k={pm[5]})",pm,"#999999")):
    yy+=.6; ax.plot([v[1],v[2]],[yy,yy],color=c,lw=2.6); ax.scatter(v[0],yy,marker="D",s=60,color=c,zorder=3)
ax.set_yticks(list(range(len(s2)))+[len(s2)+.6,len(s2)+1.2])
ax.set_yticklabels([f"{r.cohorte} · {r.tejido}" for _,r in s2.iterrows()]+
                   [f"vivos: {vv[0]:.2f}, P {vv[4]:.0e}",f"post mortem: {pm[0]:.2f}, P {pm[4]:.2f}"],fontsize=7)
ax.invert_yaxis(); ax.axvline(0,color="k",lw=.8)
ax.set(xlabel="ρ señalización glucocorticoide ~ reloj",title="c  La señalización se mide en vivos")
# ---- d dentro del donante
ax=fig.add_subplot(gs[1,1])
W=pd.read_csv(f"{OUT}/new/within_donor_decomposition.tsv",sep="\t")
dat=[(f"{r.cohorte}\n{r.tejido}",r.rho,r.rho_sin_nivel_donante) for _,r in W.iterrows()]
x=np.arange(len(dat))
ax.bar(x-.2,[d[1] for d in dat],.4,color="#999999",label="crudo")
ax.bar(x+.2,[d[2] for d in dat],.4,color="#2a9d8f",label="sin el nivel global del donante")
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels([d[0] for d in dat],fontsize=7)
ax.set(ylabel="ρ NR3C1 ~ reloj",title="d  La mayor parte del efecto es del individuo,\nno del tejido")
ax.legend(fontsize=7,frameon=False)
# ---- e islotes cultivados
ax=fig.add_subplot(gs[1,2])
C2=pd.read_csv(f"{os.environ.get('T2D_RAW','data/raw')}/GSE159984_raw_counts.tsv.gz",sep="\t",index_col=0)
C2=C2[[g for g in C2.columns if int(g[3:])>=4852119]]
C2.index=[symmap.get(str(i),str(i)) for i in C2.index]; C2=C2[~C2.index.str.isdigit()]; C2=C2.groupby(level=0).sum(); C2=C2[(C2>5).sum(1)>=10]
X2=np.log2(C2/C2.sum()*1e6+1); co2=sc(X2,SETS["Clock output"]["genes"])
r=stats.spearmanr(X2.loc["NR3C1"],co2)
ax.scatter(X2.loc["NR3C1"],co2,s=22,color=COL["islotes"],alpha=.75)
ax.set(xlabel="NR3C1 (log2 CPM)",ylabel="salida del reloj (z)",
       title=f"e  Islotes cultivados, sin eje HPA ni cortisol\ncirculante: ρ = {r.statistic:.2f}, P = {r.pvalue:.1e} (n={X2.shape[1]})")
for a in fig.axes: a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9)
fig.savefig(f"{OUT}/figs/F3_eje.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F3_eje.pdf",bbox_inches="tight")
print("ok")
