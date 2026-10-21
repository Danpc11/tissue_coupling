import pandas as pd, numpy as np, warnings, sys; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
R=pd.read_csv("results/new/FOXO3_regulon_ED4.tsv",sep="\t")
A=pd.read_pickle("results/new/AAGMEx_adipose_log2.pkl"); Ac=A.sub(A.mean(0),axis=1)
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
co=sc(Ac,SETS["Clock output"]["genes"])
R["rho"]=[stats.spearmanr(Ac.loc[g],co).statistic if g in Ac.index else np.nan for g in R.target]
R=R.sort_values(["sign","rho"],ascending=[False,False]).reset_index(drop=True)
n=len(R); ang=np.linspace(0,2*np.pi,n,endpoint=False)+np.pi/2
rad=np.where(np.arange(n)%2==0,1.02,1.30)
x,y=rad*np.cos(ang),rad*np.sin(ang)
fig,ax=plt.subplots(figsize=(9.2,9.2)); ax.set_aspect("equal"); ax.axis("off")
AD,MU,IS="#2a9d8f","#d95f02","#6a3d9a"
HI={"TXNIP","FBXO32","CDKN1A","GADD45A"}
for i,r in R.iterrows():
    ax.plot([0,x[i]],[0,y[i]],color=AD if r.sign>0 else MU,lw=1.0,alpha=.5,zorder=1)
for i,r in R.iterrows():
    face=plt.cm.RdBu_r((np.clip(r.rho,-.5,.5)+.5)) if np.isfinite(r.rho) else "#eeeeee"
    hi=r.target in HI
    ax.scatter(x[i],y[i],s=760 if hi else 620,color=face,edgecolors=IS if hi else "k",
               linewidths=2.2 if hi else .6,zorder=3)
    lum=0.299*face[0]+0.587*face[1]+0.114*face[2] if isinstance(face,tuple) else 1
    ax.text(x[i],y[i],r.target,fontsize=6.4,ha="center",va="center",zorder=4,
            color="w" if lum<0.45 else "k",fontweight="bold" if hi else "normal")
ax.scatter([0],[0],s=5200,color="#f7f7f7",edgecolors="k",linewidths=1.5,zorder=3)
ax.text(0,0.05,"FOXO3",fontsize=15,ha="center",va="center",fontweight="bold",zorder=4)
ax.text(0,-0.10,"38 targets",fontsize=8,ha="center",va="center",color="#666",zorder=4)
h=[plt.Line2D([],[],color=AD,lw=2.2),plt.Line2D([],[],color=MU,lw=2.2),
   plt.Line2D([],[],marker="o",ls="",mfc="#eeeeee",mec=IS,mew=2.2,ms=11)]
ax.legend(h,["activation","repression","gene referred to in the text"],frameon=False,
          loc="lower center",bbox_to_anchor=(0.5,-0.10),ncol=3,fontsize=8)
sm=plt.cm.ScalarMappable(cmap="RdBu_r",norm=plt.Normalize(-.5,.5))
cb=fig.colorbar(sm,ax=ax,shrink=.42,pad=0.01,aspect=16)
cb.set_label("ρ with clock-output score (AAGMEx adipose, n = 232)",fontsize=7.5); cb.ax.tick_params(labelsize=7)
ax.text(-1.52,1.46,"A",fontsize=11,fontweight="bold")
ax.set_xlim(-1.55,1.55); ax.set_ylim(-1.55,1.5)
fig.savefig("results/figs/ED4.pdf",bbox_inches="tight"); fig.savefig("results/figs/ED4.png",dpi=170,bbox_inches="tight")
print("ok")
