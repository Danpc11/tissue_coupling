#!/usr/bin/env python3
"""Extended Data Fig. 4: regulatory architecture of glucocorticoid signalling in adipose tissue.

El eje usado aqui es el score de nueve genes de senalizacion glucocorticoide, el mismo que en el
resto del manuscrito, y NO el eje de concordancia entre tejidos del envio anterior. Con este score
las correlaciones son menores, de modo que el umbral de dianas es 0.30 en lugar de 0.45: con 0.45
sobreviven tres dianas, con 0.35 diez y con 0.30 dieciseis. El umbral se declara en el pie.

A  chord diagram: transcription factors whose inferred activity correlates with the glucocorticoid
   axis score, and the targets they share.
B  correlation of each FOXO3 target with the glucocorticoid axis score.
The ranking of TF activity against clock output is Fig. 4A and is not repeated here.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results"); RAW=os.environ.get("T2D_RAW","data/raw")
POS,NEG="#c0392b","#2980b9"; ORA,BLU="#d95f02","#1f77b4"

A=pd.read_pickle(f"{OUT}/new/AAGMEx_adipose_log2.pkl"); Ac=A.sub(A.mean(0),axis=1)
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
gc=sc(Ac,SETS["Glucocorticoid"]["genes"])
d=pd.read_csv(f"{RAW}/dorothea_ABC.tsv",sep="\t",low_memory=False,
              usecols=["source_genesymbol","target_genesymbol","is_stimulation","is_inhibition"])
d["mor"]=np.where(d.is_stimulation&~d.is_inhibition,1,np.where(d.is_inhibition&~d.is_stimulation,-1,0))
d=d[d.mor!=0].drop_duplicates(["source_genesymbol","target_genesymbol"])
def activity(tf):
    r=d[d.source_genesymbol==tf]
    g=[(resolve(t,Ac.index),m) for t,m in zip(r.target_genesymbol,r.mor)]
    g=[(a,b) for a,b in g if a]
    if len(g)<8: return None
    Z=Ac.loc[[a for a,_ in g]]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    w=np.array([b for _,b in g]); return (Z.values*w[:,None]).sum(0)/np.sqrt(len(w))
act={}
for tf in d.source_genesymbol.unique():
    a=activity(tf)
    if a is not None: act[tf]=stats.spearmanr(a,gc).statistic
S=pd.Series(act).sort_values()
# Los cinco factores con mayor correlacion de actividad inferida con el eje glucocorticoide
TFS=list(S.reindex(S.abs().sort_values(ascending=False).index).head(5).index)
rho_t={g:stats.spearmanr(Ac.loc[g],gc).statistic for g in Ac.index if g in set(d.target_genesymbol)}
# Dianas curadas de esos factores que correlacionan por encima de 0.45 en valor absoluto
THR=float(os.environ.get("ED4_THR","0.30"))
conn=set(d[d.source_genesymbol.isin(TFS)].target_genesymbol)
TG=[g for g,v in sorted(rho_t.items(),key=lambda x:-abs(x[1])) if g in conn and abs(v)>THR]
print(f"{len(TFS)} TF, {len(TG)} dianas con |rho| > {THR}")

fig=plt.figure(figsize=(16,10)); gs=fig.add_gridspec(1,2,width_ratios=[2.35,1],wspace=.18)
# ---- A chord
ax=fig.add_subplot(gs[0,0]); ax.set_aspect("equal"); ax.axis("off")
nT,nG=len(TFS),len(TG)
aT=np.linspace(162,18,nT); aG=np.linspace(-14,-166,nG)
pos={}
for i,t in enumerate(TFS):
    a0,a1=aT[i]-11,aT[i]+11
    ax.add_patch(Wedge((0,0),1.0,a0,a1,width=.085,color=POS if S[t]>0 else NEG,zorder=3))
    pos[t]=np.deg2rad(aT[i]); aa=np.deg2rad(aT[i]); r=1.17
    ax.text(r*np.cos(aa),r*np.sin(aa),f"{t}\n{S[t]:+.2f}",fontsize=11,ha="center",va="center",fontweight="bold")
for i,g in enumerate(TG):
    a0,a1=aG[i]-3.2,aG[i]+3.2
    ax.add_patch(Wedge((0,0),1.0,a0,a1,width=.065,color=ORA if rho_t[g]>0 else BLU,zorder=3))
    pos[g]=np.deg2rad(aG[i]); aa=np.deg2rad(aG[i]); r=1.05
    ax.text(r*np.cos(aa),r*np.sin(aa),g,fontsize=8.5,rotation=aG[i]-180 if aG[i]<-90 else aG[i],
            rotation_mode="anchor",ha="right" if aG[i]<-90 else "left",va="center")
for _,e in d[(d.source_genesymbol.isin(TFS))&(d.target_genesymbol.isin(TG))].iterrows():
    a,b=pos[e.source_genesymbol],pos[e.target_genesymbol]
    p0=(0.91*np.cos(a),0.91*np.sin(a)); p1=(0.93*np.cos(b),0.93*np.sin(b))
    ax.add_patch(PathPatch(Path([p0,(0,0),p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3]),fill=False,
                 lw=1.0,edgecolor=POS if e.mor>0 else NEG,alpha=.5,zorder=1))
ax.set_xlim(-1.42,1.42); ax.set_ylim(-1.52,1.34)
ax.text(-1.40,1.30,"A",fontsize=13,fontweight="bold")
h=[plt.Line2D([],[],color=POS,lw=2.4),plt.Line2D([],[],color=NEG,lw=2.4),
   plt.Line2D([],[],marker="s",ls="",color=ORA,ms=9),plt.Line2D([],[],marker="s",ls="",color=BLU,ms=9)]
ax.legend(h,["factor activity, positive; activating edge","factor activity, negative; repressing edge",
             "target, positive correlation","target, negative correlation"],
          frameon=False,fontsize=8.5,loc="lower center",bbox_to_anchor=(0.5,-0.055),ncol=2)
ax.text(0,-1.49,f"{nT} transcription factors; {nG} curated targets with |ρ| > {THR:.2f}",
        fontsize=8,ha="center",color="#666666",style="italic")
# ---- B  todas las dianas de FOXO3 en un solo panel
R=pd.read_csv(f"{OUT}/new/FOXO3_regulon_ED4.tsv",sep="\t")
R["rho"]=[stats.spearmanr(Ac.loc[g],gc).statistic if g in Ac.index else np.nan for g in R.target]
R=R.dropna(subset=["rho"]).sort_values("rho",ascending=False)
HI={"TXNIP","FBXO32"}
a=fig.add_subplot(gs[0,1])
a.barh(range(len(R)),R.rho,color=[ORA if v>0 else BLU for v in R.rho])
a.set_yticks(range(len(R))); a.set_yticklabels(R.target,fontsize=7.2)
for i,g in enumerate(R.target):
    if g in HI: a.get_yticklabels()[i].set_fontweight("bold")
a.invert_yaxis(); a.axvline(0,color="k",lw=.8); a.set_xlim(-0.45,0.78)
a.set_xlabel("ρ with the glucocorticoid signalling score",fontsize=8.5)
a.spines[["top","right"]].set_visible(False)
a.text(-0.42,1.02,"B",transform=a.transAxes,fontsize=13,fontweight="bold")
a.text(0,1.02,"FOXO3 regulon, 38 targets",transform=a.transAxes,fontsize=8,color="#666666",style="italic")
for g,lab in (("TXNIP","insulin core"),("FBXO32","atrophy")):
    if g in list(R.target):
        k=list(R.target).index(g)
        a.annotate(lab,(R.rho.iloc[k],k),fontsize=7.5,xytext=(8,0),textcoords="offset points",
                   va="center",color="#555555")
fig.savefig(f"{OUT}/figs/ED4.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/figs/ED4.png",dpi=400,bbox_inches="tight")
print("ok", len(TFS), len(TG))
